"""Some helper methods to process cases from Strocchi and Rodero databases."""

import logging

LOGGER = logging.getLogger("pyheart_global.simulator")


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
        from ansys.heart.preprocessor.models.v0_2.database_labels_to_id import (
            Strocchi2020 as database_labels,
        )
    elif database == "Rodero2021":
        from ansys.heart.preprocessor.models.v0_2.database_labels_to_id import (
            Rodero2021 as database_labels,
        )
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

    if database == "Rodero2021":
        mesh.rename_array("ID", "tags", preference="cell")

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

        # LOGGER.debug(pd)

        polydatas += [pd]

    # combine polydata's into one.
    all_p = polydatas[0]
    for p in polydatas[1:]:
        all_p += p

    # re-compute tag to label dict
    tag_to_label = {v: "-".join(k.split(" ")) for k, v in labels.items()}

    # mesh as multiple parts
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

            if len(sub_sub_geoms) == 3 and "left-ventricle" in label:
                names = ["septum", "endocardium", "epicardium"]
            elif len(sub_sub_geoms) == 2:
                names = ["endocardium", "epicardium"]
            elif len(sub_sub_geoms) > 3 and "left-ventricle" in label:
                LOGGER.debug("More surfaces than expected. Naming largest three")
                names = ["unknown-surface"] * (len(sub_sub_geoms) - 3) + [
                    "septum",
                    "endocardium",
                    "epicardium",
                ]
            elif len(sub_sub_geoms) > 2:
                names = ["unknown-surface"] * (len(sub_sub_geoms) - 2) + [
                    "endocardium",
                    "epicardium",
                ]

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

    tags = copy.deepcopy(geom_all.cell_data["tags"])
    geom_all.cell_data.set_scalars(name="surface-id", scalars=np.array(tags, dtype=int))

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

    # remove plane and inlet parts from the part definitions.
    part_definitions1 = {
        k: v
        for k, v in part_definitions.items()
        # if "myocardium" in k or "aorta" in k or "ventricle" in k or "pulmonary-artery" in k
        if not any(x in k for x in ["plane", "inlet"])
    }

    # rename:
    part_definitions1["Left ventricle"] = part_definitions1.pop("left-ventricle-myocardium")
    part_definitions1["Right ventricle"] = part_definitions1.pop("right-ventricle-myocardium")
    part_definitions1["Left atrium"] = part_definitions1.pop("left-atrium-myocardium")
    part_definitions1["Right atrium"] = part_definitions1.pop("right-atrium-myocardium")
    part_definitions1["Aorta"] = part_definitions1.pop("aorta-wall")
    part_definitions1["Pulmonary artery"] = part_definitions1.pop("pulmonary-artery-wall")

    # merge border/vein parts into atria
    part_merge_map = {
        "Left atrium": [
            "left-atrial-appendage-border",
            "left-superior-pulmonary-vein-border",
            "left-inferior-pulmonary-vein-border",
            "right-inferior-pulmonary-vein-border",
            "right-superior-pulmonary-vein-border",
        ],
        "Right atrium": ["superior-vena-cava-border", "inferior-vena-cava-border"],
    }
    for target_part, source_parts in part_merge_map.items():
        for source_part in source_parts:
            part_definitions1[target_part]["enclosed_by_boundaries"].update(
                part_definitions1[source_part]["enclosed_by_boundaries"]
            )
            del part_definitions1[source_part]

    # rename septum
    part_definitions1["Left ventricle"]["enclosed_by_boundaries"][
        "right-ventricle-septum"
    ] = part_definitions1["Left ventricle"]["enclosed_by_boundaries"].pop("left-ventricle-septum")

    # remove left atrial septal inlet boundary
    try:
        del part_definitions1["Left atrium"]["enclosed_by_boundaries"][
            "left-atrium-appendage-inlet"
        ]
    except KeyError:
        pass

    if model_type == "BiVentricle":
        del part_definitions1["Left atrium"]
        del part_definitions1["Right atrium"]
        del part_definitions1["Aorta"]
        del part_definitions1["Pulmonary artery"]

    if model_type == "FourChamber":
        del part_definitions1["Aorta"]
        del part_definitions1["Pulmonary artery"]

    return (geom_all_int, part_definitions1)
