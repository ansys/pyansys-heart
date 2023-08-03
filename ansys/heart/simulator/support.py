"""Provide methods to support workflows."""
import os
import pathlib as Path

from ansys.heart.custom_logging import LOGGER
import ansys.heart.preprocessor.models_new as models


def preprocess_model(info: models.ModelInfo, model_type: str = None, clean_workdir=True):
    """Preprocess a model using model info as input."""
    path_to_model = os.path.join(info.workdir, "heart_model.pickle")
    info.path_to_model = path_to_model

    LOGGER.info("##############################")
    LOGGER.info("## Launching preprocessor: ###")
    LOGGER.info("##############################")
    LOGGER.info("## Remeshing: {:<13} ##".format(str(True)))
    LOGGER.info("##      Size: {:<13} ##".format(info.mesh_size))
    LOGGER.info("##############################")
    LOGGER.info("## working directory: %s" % info.workdir)
    LOGGER.info("## store model in: %s" % info.path_to_model)
    LOGGER.info("##############################")

    if not os.path.isdir(info.workdir):
        os.makedirs(info.workdir)

    if model_type == "BiVentricle":
        model = models.BiVentricle(info)
    elif model_type == "LeftVentricle":
        model = models.LeftVentricle(info)
    elif model_type == "FourChamber":
        model = models.FourChamber(info)
    elif model_type == "FullHeart":
        model = models.FullHeart(info)

    model._input.as_single_polydata.save(os.path.join(info.workdir, "input_polydata.vtp"))

    model.load_input()
    model.mesh_volume()
    model._update_parts()
    model.dump_model()

    model.print_info()
    if clean_workdir:
        model.info.clean_workdir([".stl", ".vtk", ".jou", ".log"])

    return model


def get_input_geom_and_part_defintions_from_public_database(
    mesh_path, model_type: str = "BiVentricle", database: str = "Strocchi2020"
):
    """Get the input geometry and part definitiosn from strocchi."""
    import copy

    from ansys.heart.preprocessor.mesh.connectivity import face_tetra_connectivity
    import numpy as np
    import pyvista as pv

    # read database labels
    if database == "Strocchi2020":
        from ansys.heart.preprocessor.database_labels import Strocchi2020 as database_labels
    else:
        print("No other databases included yet.")
        return

    # read full unstructured grid.
    mesh: pv.UnstructuredGrid = pv.read(mesh_path)
    if isinstance(mesh, pv.MultiBlock):
        mesh = mesh[0]

    if not isinstance(mesh, pv.UnstructuredGrid):
        raise TypeError("Expecting unstructured grid. Check inputs.")
    else:
        mesh: pv.UnstructuredGrid = mesh

    # remove caps and spaces in keys
    labels = {"-".join(k.lower().split()): v for k, v in database_labels.items()}
    labels_original = copy.deepcopy(labels)

    tag_to_label = {v: "-".join(k.lower().split(" ")) for k, v in labels.items()}

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
        surface_id = np.max(list(labels.values())) + 1
        labels[name] = surface_id

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
    tag_to_label = {v: "-".join(k.split(" ")) for k, v in labels.items()}

    # mesh as multiple parts
    part_definitions = {}
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
    for original_label, original_tag in labels_original.items():
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

    if model_type == "BiVentricle":
        del part_definitions1["Left atrium"]
        del part_definitions1["Right atrium"]
        del part_definitions1["Aorta"]
        del part_definitions1["Pulmonary artery"]

    return (geom_all_int, part_definitions1)


def _deprecated_run_preprocessor(
    model_type: models.HeartModel,
    database: str,
    path_original_mesh: Path,
    work_directory: Path,
    path_to_model: Path = None,
    mesh_size: float = 2.0,
    add_blood_pool: bool = False,
    clean_workdir: bool = True,
):
    """Run the preprocessor with the given input arguments.

    Parameters
    ----------
    model_type : models.HeartModel
        Type of model. Valid values include: LeftVentricle, BiVentricle, FourChamber, FullHeart
    database : str
        Name of the database. Either "Strocchi2020" or "Cristobal2021"
    path_original_mesh : Path
        Path to the input mesh file
    work_directory : Path
        Working directory
    path_to_model : Path, optional
        Path to the model, by default None, writes as "heart_model.pickle"
    mesh_size : float, optional
        Size used for remeshing the volume, by default 2.0
    clean_workdir : bool, optional
        Flag indicating whether to clean the working directory, by default True
    """
    DeprecationWarning("This method will be removed in the future.")
    if not path_to_model:
        path_to_model = os.path.join(work_directory, "heart_model.pickle")

    LOGGER.info("##############################")
    LOGGER.info("## Launching preprocessor: ###")
    LOGGER.info("##############################")
    LOGGER.info("## Model type: {:<12} ##".format(model_type.__name__))
    LOGGER.info("## Database: {:<14} ##".format(database))
    LOGGER.info("## Remeshing: {:<13} ##".format(str(True)))
    LOGGER.info("##      Size: {:<13} ##".format(mesh_size))
    LOGGER.info("##############################")
    LOGGER.info("## path mesh: %s" % path_original_mesh)
    LOGGER.info("## working directory: %s" % work_directory)
    LOGGER.info("## store model in: %s" % path_to_model)
    LOGGER.info("##############################")

    if not os.path.isdir(work_directory):
        os.makedirs(work_directory)

    # instantiate model information
    info = models.ModelInfo(
        database=database,
        _deprecated_path_to_case=path_original_mesh,
        work_directory=work_directory,
        path_to_model=path_to_model,
        add_blood_pool=add_blood_pool,
    )

    info.mesh_size = mesh_size

    info.clean_workdir(remove_all=True)
    info.create_workdir()
    info.dump_info()

    if model_type == models.LeftVentricle:
        model = models.LeftVentricle(info)

    elif model_type == models.BiVentricle:
        model = models.BiVentricle(info)

    elif model_type == models.FourChamber:
        model = models.FourChamber(info)

    elif model_type == models.FullHeart:
        model = models.FullHeart(info)

    model.extract_simulation_mesh()
    model.dump_model(path_to_model)
    model.print_info()
    if clean_workdir:
        model.info.clean_workdir([".stl", ".vtk", ".jou", ".log"])

    return model
