# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Some helper methods to process cases from Strocchi and Rodero databases."""

import copy
import os
from typing import Literal

import numpy as np
import pyvista as pv

from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.utils.connectivity import face_tetra_connectivity


def _read_input_mesh(mesh_path: str, database: str) -> pv.UnstructuredGrid:
    """Read a mesh file from Rodero2021 or Strocchi2020.

    Parameters
    ----------
    mesh_path : str
        Path to the mesh file.
    database : str
        Database name.

    Returns
    -------
    pv.UnstructuredGrid
        Unstructured grid.

    Raises
    ------
    TypeError
        If the mesh fails to be imported as an UnstructuredGrid.
    """
    mesh: pv.UnstructuredGrid = pv.read(mesh_path)
    if isinstance(mesh, pv.MultiBlock):
        mesh = mesh[0]

    if not isinstance(mesh, pv.UnstructuredGrid):
        raise TypeError("Expecting unstructured grid. Check inputs.")

    else:
        mesh: pv.UnstructuredGrid = mesh

    if database == "Rodero2021":
        mesh.rename_array("ID", "tags", preference="cell")

    return mesh


def _get_original_labels(database: str, case_num: int = None) -> dict:
    """Import the original labels based on a database name.

    Parameters
    ----------
    database : str
        Database name.

    Returns
    -------
    dict
        Dictionary representing the label to ID map.
    """
    match database:
        case "Strocchi2020":
            if case_num in [12, 14]:
                database_labels = _Strocchi2020_case_12_14_labels
            else:
                database_labels = _Strocchi2020_labels

        case "Rodero2021":
            database_labels = _Rodero2021_labels

        case _:
            LOGGER.error(f"Database with name {database} not supported.")
            return None

    return database_labels


def _get_interface_surfaces(
    mesh: pv.UnstructuredGrid, labels: dict, tag_to_label: dict
) -> tuple[list[pv.PolyData], dict]:
    """Get the each of the interface surfaces as polydata.

    Parameters
    ----------
    mesh : pv.UnstructuredGrid
        Volume mesh.
    labels : dict
        Label dictionary to add the interface labels to.
    """
    tetras = np.reshape(mesh.cells, (mesh.n_cells, 5))[:, 1:]
    faces, c0, c1 = face_tetra_connectivity(tetras)
    # NOTE: if c0 != c1, face is an interface between two parts.
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

    interfaces = []
    # extract interfaces:
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

        interface = pv.PolyData(mesh.points, faces=faces_interface)
        interface.cell_data.set_scalars(
            name="surface-id", scalars=np.ones(interface.n_cells, dtype=int) * surface_id
        )

        interfaces += [interface]

    return interfaces, labels


def _find_endo_epicardial_regions(
    geom_all: pv.PolyData, tag_to_label: dict
) -> tuple[pv.PolyData, dict]:
    """Find the endo and epicardial regions from the surface geometry of the entire model.

    Parameters
    ----------
    geom_all : pv.PolyData
        Entire heart model.
    tag_to_label : dict
        Dictionary that maps the tags to the corresponding labels.
    """
    geom_all.cell_data["orig_ids"] = np.arange(0, geom_all.n_cells)

    tag_offset = max(tag_to_label.keys()) + 1

    new_tag_to_label = copy.deepcopy(tag_to_label)

    for tag, label in tag_to_label.items():
        # split myocardial surfaces
        if "myocardium" in label and "interface" not in label:
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

    return geom_all, new_tag_to_label


def _get_part_definitions(original_labels: dict, boundary_label_to_boundary_id: dict) -> dict:
    """Format the part definitions based on the original labels and the boundary labels.

    Parameters
    ----------
    original_labels : dict
        Dictionary with the original labels.
    boundary_label_to_boundary_id : dict
        Dictionary of the boundary label to boundary ID map.

    Returns
    -------
    dict
        Dictionary with the part definitions, which is the part ID and corresponding
        boundaries that enclose that part.
    """
    part_definitions = {}
    for original_label, original_tag in original_labels.items():
        # boundary_names = [original_label] + interface_keys
        if "myocardium" in original_label:
            part_label = original_label.replace("-myocardium", "")
        else:
            part_label = original_label

        enclosed_by_boundaries = {
            label: int(boundary_label_to_boundary_id[label])
            for label in boundary_label_to_boundary_id
            if part_label in label
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
    part_definitions1["Left ventricle"]["enclosed_by_boundaries"]["right-ventricle-septum"] = (
        part_definitions1["Left ventricle"]["enclosed_by_boundaries"].pop("left-ventricle-septum")
    )

    # remove left atrial septal inlet boundary
    try:
        del part_definitions1["Left atrium"]["enclosed_by_boundaries"][
            "left-atrium-appendage-inlet"
        ]
    except KeyError:
        pass

    return part_definitions1


def _sort_edge_loop(edges) -> np.ndarray:
    """Sorts the points in an edge loop."""
    remaining_edges = edges
    next_edge = edges[0]
    sorted_edges = [next_edge]
    remaining_edges.pop(0)
    while len(remaining_edges) > 0:
        # find connected edge of last edge
        node = sorted_edges[-1][1]
        mask = np.array(edges) == node
        if np.sum(mask[:, 1]) == 1:
            flip = True
        elif np.sum(mask[:, 0]) == 1:
            flip = False
        else:
            raise ValueError("Expecting just one match")
        idx = np.where(np.any(mask, axis=1))[0][0]
        if flip:
            sorted_edges.append(np.flip(remaining_edges[idx]).tolist())
        else:
            sorted_edges.append(remaining_edges[idx])
        remaining_edges.pop(idx)
    return np.array(sorted_edges)


def _smooth_boundary_edges(
    surface_mesh: pv.PolyData,
    id_to_label_map,
    sub_label_to_smooth: str = "endocardium",
    window_size: int = 5,
    project_edge_loop: bool = True,
) -> tuple[pv.PolyData, list]:
    """Smooth edges of surfaces that match the label string.

    Parameters
    ----------
    surface_mesh : pv.PolyData
        Input surface mesh.
    id_to_label_map : dict
        ID to label map.
    sub_label_to_smooth : str, default: ``'endocardium'``
        Select labels where this sub string is present.
    window_size : int, default: 5
        Window size of the smoothing method.
    project_edge_loop : bool, default: True
        Whether to project the edge loop to a repesentative plane before smoothing.

    Returns
    -------
    Tuple[pv.PolyData, dict]
        Preprocessor-compatible polydata object and dictionary with part definitions.
    """
    surfaces_to_smooth = [
        id for id, label in id_to_label_map.items() if sub_label_to_smooth in label
    ]

    surface_mesh.point_data["original-point-ids"] = np.arange(0, surface_mesh.n_points)

    all_edges = []
    for surf_id in surfaces_to_smooth:
        print("Processing " + id_to_label_map[surf_id])

        mask = surface_mesh.cell_data["surface-id"] == surf_id
        sub_surface = surface_mesh.extract_cells(mask).extract_surface()
        # get edges
        edges = sub_surface.extract_feature_edges(
            boundary_edges=True, non_manifold_edges=False, feature_edges=False, manifold_edges=False
        )
        conn = edges.connectivity()
        for region_id in np.unique(conn.cell_data["RegionId"]):
            mask1 = conn.cell_data["RegionId"] == region_id
            edges = conn.extract_cells(mask1).extract_surface()

            # only project if we have a manifold edge.
            # NOTE: this doesn't seem to do much.
            if edges.is_manifold:
                # ensure points are ordered correctly.

                # use a window average to "smooth" edge loop
                # assumes ordering is ok.
                edges_array = np.reshape(edges.lines, (edges.n_cells, 3))[:, 1:].tolist()
                try:
                    sorted_edges_array = _sort_edge_loop(edges_array)
                except Exception:
                    print(f"Failed to sort edges for {id_to_label_map[surf_id]} region {region_id}")
                    continue

                # project points
                if project_edge_loop:
                    # fit points through plane.
                    _, plane_center, plane_normal = pv.fit_plane_to_points(
                        edges.points, return_meta=True
                    )
                    # project points to the plane.
                    new_points = (
                        pv.PolyData(edges.points)
                        .project_points_to_plane(plane_center, plane_normal)
                        .points
                    )
                    edges.points = new_points

                sorted_points = edges.points[sorted_edges_array[:, 0], :]

                # smooth with window size
                num_points_to_add = int((window_size - 1) / 2)
                sorted_points
                sorted_points = np.concatenate(
                    (
                        sorted_points[-num_points_to_add:],
                        sorted_points,
                        sorted_points[0:num_points_to_add],
                    )
                )
                offset = num_points_to_add
                for ii, node in enumerate(sorted_points[:-num_points_to_add]):
                    sorted_points[ii + offset] = np.mean(
                        sorted_points[ii : ii + window_size, :], axis=0
                    )

                sorted_points = sorted_points[num_points_to_add:-num_points_to_add]

                edges.points[sorted_edges_array[:, 0]] = sorted_points

                surface_mesh.points[edges.point_data["original-point-ids"]] = copy.deepcopy(
                    edges.points
                )

            all_edges.append(edges)

    return surface_mesh, all_edges


def get_compatible_input(
    mesh_path: str,
    model_type: Literal["FullHeart", "FourChamber", "BiVentricle", "LeftVentricle"] = "FullHeart",
    database: str = "Rodero2021",
) -> tuple[pv.PolyData, dict]:
    """Extract a preprocessor-compatible input surface.

    Parameters
    ----------
    mesh_path : str
        Path to the input mesh (UnstructuredGrid or MultiBlock).
    model_type : str, default: ``'FullHeart'``
        Type of model to extract. Options are ``'FullHeart'``, ``'FourChamber'``,
        ``'BiVentricle'``, and ``'LeftVentricle'``.
    database : str, default: ``'Rodero2021'``
        Database name. Options are ``'Rodero2021'`` and ``'Strocchi2020'``.

    Returns
    -------
    Tuple[pv.PolyData, dict]
        Preprocessor-compatible polydata object and dictionary with part definitions.
    """
    case_num = os.path.basename(mesh_path)
    case_num = int(case_num.replace(".case", "").replace(".vtk", ""))
    # get the original label <> id map
    database_labels = _get_original_labels(database, case_num)

    # read the mesh file.
    mesh = _read_input_mesh(mesh_path, database)

    # normalize label strings
    labels_to_tags = {"-".join(k.lower().split()): v for k, v in database_labels.items()}
    tags_to_label = {v: k for k, v in labels_to_tags.items()}

    labels_original = copy.deepcopy(labels_to_tags)

    # get interfaces between the different parts as surfaces
    interfaces, labels_to_tags = _get_interface_surfaces(mesh, labels_to_tags, tags_to_label)

    # combine polydata's into one.
    all_interfaces_as_polydata = interfaces[0]
    for interface in interfaces[1:]:
        all_interfaces_as_polydata += interface

    # Update tag to label dict
    tags_to_label = {v: "-".join(k.split(" ")) for k, v in labels_to_tags.items()}

    # extract surface of mesh - this is used to find the endo and epicardial
    # regions
    mesh_surface = mesh.extract_geometry()

    # find the endo and epicardial regions
    mesh_surface, new_tag_to_label = _find_endo_epicardial_regions(mesh_surface, tags_to_label)

    # update the the label to tag dictionary
    label_to_tag = {v: k for k, v in new_tag_to_label.items()}

    # Store surface "topology" in "surface-id"
    tags = copy.deepcopy(mesh_surface.cell_data["tags"])
    mesh_surface.cell_data.set_scalars(name="surface-id", scalars=np.array(tags, dtype=int))

    # combine interfaces and all surfaces of the model into single polydata.
    geom_with_interfaces = mesh_surface + all_interfaces_as_polydata

    # get the part definitions from the labels that are defined.
    part_definitions = _get_part_definitions(labels_original, label_to_tag)

    # delete parts of dictionary depending on the requested model.
    if model_type == "LeftVentricle":
        # For the LeftVentricle model the epicardium consists of the
        # following surfaces, so we need to merge these:
        # - left-ventricle-epicardium
        # - right-ventricle-septum
        # - interface_left-ventricle-myocardium_right-ventricle-myocardium
        epi_ids = [
            part_definitions["Left ventricle"]["enclosed_by_boundaries"][
                "interface_left-ventricle-myocardium_right-ventricle-myocardium"
            ],
            part_definitions["Left ventricle"]["enclosed_by_boundaries"]["right-ventricle-septum"],
            part_definitions["Left ventricle"]["enclosed_by_boundaries"][
                "left-ventricle-epicardium"
            ],
        ]
        mask = np.isin(geom_with_interfaces.cell_data["surface-id"], epi_ids)
        geom_with_interfaces.cell_data["surface-id"][mask] = epi_ids[-1]
        del part_definitions["Left ventricle"]["enclosed_by_boundaries"]["right-ventricle-septum"]
        del part_definitions["Left ventricle"]["enclosed_by_boundaries"][
            "interface_left-ventricle-myocardium_right-ventricle-myocardium"
        ]

        del part_definitions["Right ventricle"]
        del part_definitions["Left atrium"]
        del part_definitions["Right atrium"]
        del part_definitions["Aorta"]
        del part_definitions["Pulmonary artery"]

    if model_type == "BiVentricle":
        del part_definitions["Left atrium"]
        del part_definitions["Right atrium"]
        del part_definitions["Aorta"]
        del part_definitions["Pulmonary artery"]

    elif model_type == "FourChamber":
        del part_definitions["Aorta"]
        del part_definitions["Pulmonary artery"]

    return geom_with_interfaces, part_definitions


# Dictionaries to map labels to ID for Strocchi et al 2020 and Rodero et al 2021 datasets."""
_Strocchi2020_labels = {
    "Left ventricle myocardium": 1,
    "Right ventricle myocardium": 2,
    "Left atrium myocardium": 3,
    "Right atrium myocardium": 4,
    "Aorta wall": 5,
    "Pulmonary artery wall": 6,
    "Left atrial appendage border": 7,
    "Left superior pulmonary vein border": 8,
    "Left inferior pulmonary vein border": 9,
    "Right inferior pulmonary vein border": 10,
    "Right superior pulmonary vein border": 11,
    "Superior vena cava border": 12,
    "Inferior vena cava border": 13,
    "Mitral valve plane": 14,
    "Tricuspid valve plane": 15,
    "Aortic valve plane": 16,
    "Pulmonary valve plane": 17,
    "Left atrium appendage inlet": 18,
    "Left superior pulmonary vein inlet": 19,
    "Left inferior pulmonary vein inlet": 20,
    "Right inferior pulmonary vein inlet": 21,
    "Right superior pulmonary vein inlet": 22,
    "Superior vena cava inlet": 23,
    "Inferior vena cava inlet": 24,
}
#! NOTE: Strocchi 2020 cases 12 and 14 have different labeling.
_Strocchi2020_case_12_14_labels = {
    "Left ventricle myocardium": 1,
    "Right ventricle myocardium": 2,
    "Left atrium myocardium": 3,
    "Right atrium myocardium": 4,
    "Aorta wall": 5,
    "Pulmonary artery wall": 6,
    "Left atrial appendage border": 7,
    "Left superior pulmonary vein border": 8,
    "Left inferior pulmonary vein border": 9,
    "Right inferior pulmonary vein border": 9,
    "Right superior pulmonary vein border": 10,
    "Superior vena cava border": 11,
    "Inferior vena cava border": 12,
    "Mitral valve plane": 13,
    "Tricuspid valve plane": 14,
    "Aortic valve plane": 15,
    "Pulmonary valve plane": 16,
    "Left atrium appendage inlet": 17,
    "Left superior pulmonary vein inlet": 18,
    "Left inferior pulmonary vein inlet": 18,
    "Right inferior pulmonary vein inlet": 20,
    "Right superior pulmonary vein inlet": 19,
    "Superior vena cava inlet": 21,
    "Inferior vena cava inlet": 22,
}
_Rodero2021_labels = {
    "Left ventricle myocardium": 1,
    "Right ventricle myocardium": 2,
    "Left atrium myocardium": 3,
    "Right atrium myocardium": 4,
    "Aorta wall": 5,
    "Pulmonary artery wall": 6,
    "Left atrial appendage border": 18,
    "Left superior pulmonary vein border": 21,
    "Left inferior pulmonary vein border": 20,
    "Right inferior pulmonary vein border": 19,
    "Right superior pulmonary vein border": 22,
    "Superior vena cava border": 23,
    "Inferior vena cava border": 24,
    "Mitral valve plane": 7,
    "Tricuspid valve plane": 8,
    "Aortic valve plane": 9,
    "Pulmonary valve plane": 10,
    "Left atrium appendage inlet": 11,
    "Left superior pulmonary vein inlet": 12,
    "Left inferior pulmonary vein inlet": 13,
    "Right inferior pulmonary vein inlet": 14,
    "Right superior pulmonary vein inlet": 15,
    "Superior vena cava inlet": 16,
    "Inferior vena cava inlet": 17,
}
