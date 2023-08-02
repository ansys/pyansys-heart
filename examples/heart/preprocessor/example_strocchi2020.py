# example to process one of the examples of Strocchi2020.

# 1. extract surfaces. [pyvista/vtk]
# 2. label surfaces. [pyvista/vtk]
# 3. set input model.
# 4. call mesher
# 5. create heart model.
import os
import pathlib

from ansys.heart.misc.downloader import download_case, unpack_case
import pyvista as pv
import numpy as np

from ansys.heart.preprocessor.mesh.connectivity import (
    face_tetra_connectivity,
)


run_extraction = True
if run_extraction:
    # download case from remote repository
    case_num = 1  # patient number 1
    database = "Strocchi2020"

    download_folder: pathlib.Path = os.path.join(pathlib.Path(__file__).parents[3], "downloads")
    case_path: pathlib.Path = download_case(
        database=database, case_number=case_num, download_folder=download_folder, overwrite=False
    )
    mesh_path = os.path.join(
        pathlib.Path(case_path).parents[0], "%02d" % (case_num,), "%02d.case" % (case_num,)
    )

    if not os.path.isfile(mesh_path):
        unpack_case(case_path)

    # read full unstructured grid.
    mesh: pv.UnstructuredGrid = pv.read(mesh_path)
    if isinstance(mesh, pv.MultiBlock):
        mesh = mesh[0]

    if not isinstance(mesh, pv.UnstructuredGrid):
        raise TypeError("Expecting unstructured grid. Check inputs.")
    else:
        mesh: pv.UnstructuredGrid = mesh

    # read database labels
    from ansys.heart.preprocessor.database_labels import Strocchi2020

    # remove caps and spaces in keys
    import copy

    strocchi2020_labels = {"-".join(k.lower().split()): v for k, v in Strocchi2020.items()}
    strocchi2020_labels_original = copy.deepcopy(strocchi2020_labels)

    tag_to_label = {v: "-".join(k.lower().split(" ")) for k, v in strocchi2020_labels.items()}

    # split the surfaces by connectivity to give them anatomical meaning:
    # 1. compute intersections between the parts
    # 2. make up new id's and store in dictionary.
    # 3. prep the input: first using individual polydata's, but combine them later.

    # use face-tetra connectivity to find interface faces, and tag these as individual boundaries.
    tetras = np.reshape(mesh.cells, (mesh.n_cells, 5))[:, 1:]
    faces, c0, c1 = face_tetra_connectivity(tetras)
    interface_mask = mesh.cell_data["tags"][c0] != mesh.cell_data["tags"][c1]
    t0 = np.array(mesh.cell_data["tags"][c0], dtype=int)
    t1 = np.array(mesh.cell_data["tags"][c1], dtype=int)

    # all available tag pairs (interface pairs)
    pairs = np.unique(
        np.sort(
            np.array([mesh.cell_data["tags"][c0], mesh.cell_data["tags"][c1]], dtype=int).T, axis=1
        ),
        axis=0,
    )
    pairs = np.array([p for p in pairs if p[0] != p[1]])
    polydatas: pv.PolyData = []
    # extract pairs:
    for pair in pairs:
        mask1 = np.all(np.array([t0 == pair[0], t1 == pair[1]]), axis=0)
        mask2 = np.all(np.array([t0 == pair[1], t1 == pair[0]]), axis=0)
        mask = np.any(np.array([mask1, mask2]), axis=0)

        name = "interface_" + tag_to_label[pair[0]].lower() + "_" + tag_to_label[pair[1]].lower()
        surface_id = np.max(list(strocchi2020_labels.values())) + 1
        strocchi2020_labels[name] = surface_id

        faces_interface = np.hstack(
            [np.ones((np.sum(mask), 1), dtype=int) * 3, faces[mask, :]]
        ).flatten()

        pd = pv.PolyData(mesh.points, faces=faces_interface)
        pd.cell_data.set_scalars(
            name="surface-id", scalars=np.ones(pd.n_cells, dtype=int) * surface_id
        )

        polydatas += [pd]

    all_p = pv.PolyData()
    for p in polydatas:
        all_p += p

    # re-compute tag to label dict
    tag_to_label = {v: "-".join(k.split(" ")) for k, v in strocchi2020_labels.items()}

    # mesh as multiple parts
    part_definitions = {}
    geoms = pv.PolyData()
    geom_all = mesh.extract_geometry()

    tag_offset = max(tag_to_label.keys()) + 1
    geom_all.cell_data["orig_ids"] = np.arange(0, geom_all.n_cells)

    new_tag_to_label = copy.deepcopy(tag_to_label)
    for tag, label in tag_to_label.items():
        # split myocardial surfaces
        if "myocardium" in label and not "interface" in label:
            mask = geom_all.cell_data["tags"] == tag
            sub_geom = geom_all.extract_cells(mask).extract_geometry()
            sub_geom = sub_geom.connectivity()

            # get connected regions and sort by bounding box volume
            sub_sub_geoms = []
            for region_id in np.unique(sub_geom.cell_data["RegionId"]):
                sub_sub_geom = sub_geom.extract_cells(
                    sub_geom.cell_data["RegionId"] == region_id
                ).extract_geometry()
                sub_sub_geoms += [sub_sub_geom]

            sub_sub_geoms.sort(
                key=lambda x: (x.bounds[1] - x.bounds[0])
                * (x.bounds[3] - x.bounds[2])
                * (x.bounds[5] - x.bounds[4]),
                reverse=False,
            )

            if len(sub_sub_geoms) == 3:
                names = ["septum", "endocardium", "epicardium"]
            elif len(sub_sub_geoms) == 2:
                names = ["endocardium", "epicardium"]

            # update dictionary and geometry cell data
            for ii, sub in enumerate(sub_sub_geoms):
                geom_all.cell_data["tags"][sub.cell_data["orig_ids"]] = tag_offset
                new_tag_to_label[tag_offset] = (
                    tag_to_label[tag].replace("-myocardium", "") + "-" + names[ii]
                )
                tag_offset += 1
            # remove tag from dict.
            del new_tag_to_label[tag]

    label_to_tag = {v: k for k, v in new_tag_to_label.items()}

    # geom_all.set_scalars(name="surface-id", geom_all.cell_data[])
    geom_all.rename_array("tags", "surface-id")
    geom_all.cell_data["surface-id"] = np.array(geom_all.cell_data["surface-id"], dtype=np.int32)
    geom_all_int = geom_all + all_p

    # extract tag-id and also those of interface. These make up the part.
    # form part definitions:
    part_definitions = {}
    for original_label, original_tag in strocchi2020_labels_original.items():
        # boundary_names = [original_label] + interface_keys
        if "myocardium" in original_label:
            part_label = original_label.replace("-myocardium", "")
        else:
            part_label = original_label

        enclosed_by_boundaries = {
            label: int(label_to_tag[label]) for label in label_to_tag if part_label in label
        }

        part_definitions[original_label] = {
            "id": original_tag,
            "enclosed_by_boundaries": enclosed_by_boundaries,
        }

    from ansys.heart.preprocessor.models_new import HeartModel, ModelInfo, FullHeart, BiVentricle

    # truncate part definitions to just left, right ventricle/atrium and arteries.
    part_definitions1 = {
        k: v
        for k, v in part_definitions.items()
        if "myocardium" in k or "aorta" in k or "ventricle" in k or "pulmonary-artery" in k
    }

    part_definitions1["Left ventricle"] = part_definitions1.pop("left-ventricle-myocardium")
    part_definitions1["Right ventricle"] = part_definitions1.pop("right-ventricle-myocardium")
    part_definitions1["Left atrium"] = part_definitions1.pop("left-atrium-myocardium")
    part_definitions1["Right atrium"] = part_definitions1.pop("right-atrium-myocardium")
    part_definitions1["Aorta"] = part_definitions1.pop("aorta-wall")
    part_definitions1["Pulmonary artery"] = part_definitions1.pop("pulmonary-artery-wall")
    # add pulmonary artery.

    part_definitions1["Left ventricle"]["enclosed_by_boundaries"][
        "right-ventricle-septum"
    ] = part_definitions1["Left ventricle"]["enclosed_by_boundaries"].pop("left-ventricle-septum")

    del part_definitions1["Left atrium"]
    del part_definitions1["Right atrium"]
    del part_definitions1["Aorta"]
    del part_definitions1["Pulmonary artery"]

    info = ModelInfo(
        geom_all_int,
        scalar="surface-id",
        part_definitions=part_definitions1,
        work_directory=r"D:\development\pyheart-lib\pyheart-lib\downloads\Strocchi2020\01\test_new_model",
    )

    # model = BiVentricle(info)
    # model.load_input()
    # model._input.plot()

    from ansys.heart.simulator.support import preprocess_model

    model = preprocess_model(info, "BiVentricle", clean_workdir=False)


import os
import pathlib

import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.simulator import MechanicsSimulator
from ansys.heart.simulator.support import _deprecated_run_preprocessor

model = models.HeartModel.load_model(
    r"D:\development\pyheart-lib\pyheart-lib\downloads\Strocchi2020\01\test_new_model\heart_model.pickle"
)

lsdyna_path = (
    r"D:\development\dyna-versions\ls-dyna_smp_d_Dev_97584-g1b99fd817b_winx64_ifort190.exe"
)
simulator = MechanicsSimulator(
    model=model,
    lsdynapath=lsdyna_path,
    dynatype="smp",
    num_cpus=4,
    simulation_directory=os.path.join(model.info.workdir, "simulation-mechanics"),
)

# load default settings.
simulator.settings.load_defaults()
# compute the fiber orientation
simulator.compute_fibers()
# visualize computed fibers
simulator.model.plot_fibers(n_seed_points=2000)
# compute the stress free configuration
# simulator.compute_stress_free_configuration()
# do the main simulation
simulator.simulate()
print("done")

# from ansys.heart.simulator.simulator import MechanicsSimulator


# to do


# remove cells that are not associated with parts.
