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

"""Module containing methods for mesh operations related to the VTK library."""

import copy

import numpy as np
import pyvista as pv
import vtk

from ansys.health.heart import LOG as LOGGER


def compute_surface_nodal_area_pyvista(surface: pv.PolyData) -> np.ndarray:
    """Compute an average nodal area by summing surface areas of connected elements.

    Parameters
    ----------
    vtk_surface : vtk.vtkPolyData
        VTK object describing the object.

    Returns
    -------
    np.array
        NumPy array with nodal areas of length number of points.

    Notes
    -----
    This method adds the partial areas of connected elements/cells to each node.

    """
    num_points = surface.n_points
    nodal_area = np.zeros(num_points)
    # compute area of all cells
    n_cells = surface.n_cells
    cell_area = np.zeros(n_cells)
    for icell in range(n_cells):
        cell_area[icell] = surface.GetCell(icell).ComputeArea()
        # cell_area[icell] = vtk_surface.GetCell(icell).ComputeArea()

    # tris = get_tri_info_from_polydata(surface)[1]
    tris = np.reshape(surface.faces, (surface.n_cells, 4))[:, 1:]

    ii = 0
    for points, area in zip(tris, cell_area):
        nodal_area[points] += area / 3
        ii += 1
    return nodal_area


def extrude_polydata(
    surface: pv.PolyData,
    extrude_by: float = 1,
    extrude_direction: np.ndarray = np.empty(0),
) -> pv.PolyData:
    """Extrude a given polydata surface in a given direction.

    Parameters
    ----------
    surface : pv.PolyData
        Surface to extrude
    extrude_by : float, default: 1
        Amount to extrude.
    extrude_direction : np.array, default: np.empty(0)
        Direction of extrusion, which should have three components. If no
        components are specified, it extrudes in the normal direction.

    Returns
    -------
    pv.PolyData
        Extruded polydata object.
    """
    extrude_normal = False
    if len(extrude_direction) == 0:
        extrude_normal = True

    # NOTE: pyvista does not have a extrusion in cell normal direction. Hence
    # resort to VTK method.
    surface = surface.compute_normals()

    extrude = vtk.vtkLinearExtrusionFilter()
    extrude.CappingOn()

    if extrude_normal:
        extrude.SetExtrusionTypeToNormalExtrusion()
    else:
        extrude.SetExtrusionTypeToVectorExtrusion()
        extrude.SetVector(extrude_direction[0], extrude_direction[1], extrude_direction[2])

    extrude.SetInputData(surface)
    extrude.SetScaleFactor(extrude_by)
    extrude.Update()

    return pv.PolyData(extrude.GetOutput())


def cell_ids_inside_enclosed_surface(
    source: pv.UnstructuredGrid | pv.PolyData, surface: pv.PolyData, tolerance: float = 1e-9
) -> np.ndarray:
    """Get IDs of cells with centroids that are inside a given surface.

    Parameters
    ----------
    source : pv.UnstructuredGrid
        Source object to check which cells are inside/outside
        the specified surface.
    surface : pv.PolyData
        Surface to check whether cells are inside/outside.
    tolerance : float, default: 1e-9
        Tolerance for the ``select_enclosed_points`` filter.

    Notes
    -----
    This method also accepts a source that represents the cell centroids.
    In this case, computing the cell centers is skipped.

    Returns
    -------
    np.ndarray
        Array with cell IDs that are inside the enclosed surface.
    """
    surface = surface.compute_normals()
    if isinstance(source, pv.PolyData) and source.n_verts == source.n_points:
        centroids = source
    else:
        centroids = source.cell_centers()
    centroids = centroids.select_enclosed_points(surface, tolerance=tolerance, check_surface=True)
    cell_ids_inside = np.where(centroids.point_data["SelectedPoints"] == 1)[0]
    return cell_ids_inside


def find_cells_close_to_nodes(
    mesh: pv.UnstructuredGrid, node_ids: list[int], radius: float = 2
) -> np.ndarray:
    """Find cell IDs close to nodes.

    Parameters
    ----------
    mesh : pv.UnstructuredGrid
        Target mesh.
    node_ids : list[int]
        Node IDs.
    radius : float, default: 2
        Influence radius.

    Returns
    -------
    np.ndarray
        Cell IDs.
    """
    # Get coordinates of the given node IDs
    points = mesh.points[node_ids]

    # Create a list to store cells within the sphere radius
    selected_cells = []

    # Iterate through each point and find cells within the sphere
    for point in points:
        # Create a sphere at the given point
        sphere = pv.Sphere(radius=radius, center=point)

        # Use Boolean intersection to find cells that intersect with the sphere
        selection = mesh.select_enclosed_points(sphere, tolerance=0.0)

        # Get the indices of the cells
        selected_points = selection.point_data["SelectedPoints"].nonzero()[0]
        cells_within_sphere = mesh.extract_points(selected_points).cell_data["vtkOriginalCellIds"]

        # Store unique cell indices
        selected_cells.extend(cells_within_sphere)

    # Return unique cell indices
    return np.unique(selected_cells)


def get_boundary_edges(surface: pv.PolyData) -> pv.MultiBlock:
    """Get the boundary edges from an input surface.

    Parameters
    ----------
    surface : pv.PolyData
        Surface to check for boundary edges.

    Returns
    -------
    pv.MultiBlock
        Multi-block data, where each block represents connected edges.
    """
    surface1 = copy.deepcopy(surface)
    edge_group = surface1.extract_feature_edges(
        boundary_edges=True, non_manifold_edges=False, feature_edges=False, manifold_edges=False
    )
    # NOTE: is line ordering ensured for closed loops?
    # use connectivity filter to find connected edges
    edge_group = edge_group.connectivity()
    # split by connectivity:
    edge_groups = edge_group.split_bodies("RegionId")

    return edge_groups


def get_boundary_edge_loops(
    surface: pv.PolyData, remove_open_edge_loops: bool = True, return_types: bool = False
) -> dict | tuple[dict, dict]:
    """Get the closed/open boundary edge loops of a surface mesh.

    Parameters
    ----------
    surface : pv.PolyData
        Surface mesh to check for boundary edges.
    remove_open_edge_loops : bool, default: True
        Whether to remove open edge loops from the returned dictionary.

    Returns
    -------
    dict
        Dictionary with the edges that make up the open/closed loop.
    """
    # NOTE: Perhaps more consistent to return a pyvista polydata.

    # add cell and point IDs to keep track of IDs.
    surface1 = copy.deepcopy(surface)
    surface1.cell_data["original-cell-ids"] = np.arange(0, surface1.n_cells)
    surface1.point_data["original-point-ids"] = np.arange(0, surface1.n_points)

    # get boundary edges separated by connectivity
    edges_block = get_boundary_edges(surface1)

    # lines formed with original point IDs
    edge_groups = {
        k: edges.point_data["original-point-ids"][edges.cells_dict[3]]
        for k, edges in enumerate(edges_block)
    }

    # check if it is a closed or open edge loop. Remove open ones.
    group_types = {}
    closed_edge_groups = {}
    for k, edge_group in edge_groups.items():
        counts = np.unique(edge_group, return_counts=True)[1]
        if np.all(counts == 2):
            group_types[k] = "closed"
            closed_edge_groups[k] = edge_group
        else:
            group_types[k] = "open"

    num_open_edge_groups = len(edge_groups) - len(closed_edge_groups.keys())

    if remove_open_edge_loops and num_open_edge_groups > 0:
        LOGGER.debug(f"Removed {num_open_edge_groups} edge groups")
        return closed_edge_groups
    else:
        if return_types:
            return edge_groups, group_types
        else:
            return edge_groups


def get_patches_delaunay(surface: pv.PolyData, closed_only: bool = True) -> list[pv.PolyData]:
    """Patch boundary edges with the Delaunay algorithm.

    Parameters
    ----------
    surface : pv.PolyData
        Surface with boundary edges to find patches for.
    closed_only : bool, default: True
        Whether to return patches for closed loops of boundary edges.

    Returns
    -------
    List[pv.PolyData]
        List of patches that close the open surface.
    """
    surface1 = copy.deepcopy(surface)
    surface1.cell_data["original-cell-ids"] = np.arange(0, surface1.n_cells)
    surface1.point_data["original-point-ids"] = np.arange(0, surface1.n_points)

    # edges_block = get_boundary_edges(surface1)
    if closed_only:
        edge_groups = get_boundary_edge_loops(surface1, remove_open_edge_loops=True)
    else:
        edge_groups = get_boundary_edge_loops(surface1, remove_open_edge_loops=False)

    # reconstruct polydata objects from global ids in edge_groups
    edges_block = pv.MultiBlock()
    for edge_group in edge_groups.values():
        edges = np.array(edge_group, dtype=int)
        lines = np.hstack([np.array([2] * edges.shape[0])[:, None], edges])
        edges_block.append(pv.PolyData(surface1.points, lines=lines).clean())

    # generate patches
    patches = []
    for edges in edges_block:
        cloud = pv.PolyData(edges.points)
        patch = cloud.delaunay_2d()
        patches.append(patch)

    return patches


def get_patches_with_centroid(
    surface: pv.PolyData, closed_only: bool = True
) -> list[pv.PolyData] | None:
    """Patch boundary edges with a custom algorithm using a central node.

    Parameters
    ----------
    surface : pv.PolyData
        Surface with boundary edges to find patches for.
    closed_only : bool, default: True
        Whether to return patches for closed loops of boundary edges.

    Notes
    -----
    Edges must be sorted properly for this method to return sensible patches.

    Returns
    -------
    List[pv.PolyData]
        List of patches that close the open surface.
    """
    surface1 = copy.deepcopy(surface)
    surface1.cell_data["original-cell-ids"] = np.arange(0, surface1.n_cells)
    surface1.point_data["original-point-ids"] = np.arange(0, surface1.n_points)

    # edges_block = get_boundary_edges(surface1)
    if closed_only:
        edge_groups = get_boundary_edge_loops(surface1, remove_open_edge_loops=True)
    else:
        raise NotImplementedError
        return
        edge_groups = get_boundary_edge_loops(surface1, remove_open_edge_loops=False)

    # reconstruct polydata objects from global ids in edge_groups
    patches = []
    for edge_group in edge_groups.values():
        centroid = np.mean(surface1.points[np.unique(edge_group), :], axis=0)

        pv.PolyData(centroid)

        centroid_id = surface1.points.shape[0]

        surface1.points = np.vstack([surface1.points, centroid])

        triangles = np.hstack(
            [edge_group, np.ones(edge_group.shape[0], dtype=int)[:, None] * centroid_id]
        )
        # form input to pv.PolyData
        triangles = np.hstack([np.ones(edge_group.shape[0], dtype=int)[:, None] * 3, triangles])

        patch = pv.PolyData(surface1.points, triangles)
        patches.append(patch)

    return patches


def are_connected(
    mesh1: pv.PolyData | pv.UnstructuredGrid, mesh2: pv.PolyData | pv.UnstructuredGrid
) -> bool:
    """Check whether two polydata or unstructured grid objects are connected.

    Parameters
    ----------
    mesh1 : Union[pv.PolyData, pv.UnstructuredGrid]
        First mesh.
    mesh2 : Union[pv.PolyData, pv.UnstructuredGrid]
        Second mesh.
    """
    try:
        mesh1.cell_data.remove("RegionId")
        mesh2.cell_data.remove("RegionId")
    except KeyError:
        pass

    mesh1 = mesh1.clean().connectivity()
    n_regions_1 = len(np.unique(mesh1.cell_data["RegionId"]))
    mesh2 = mesh2.clean().connectivity()
    n_regions_2 = len(np.unique(mesh1.cell_data["RegionId"]))

    merged = (mesh1 + mesh2).clean()
    try:
        merged.cell_data.remove("RegionId")
    except KeyError:
        LOGGER.debug("RegionId not present in the cell data")

    merged = merged.connectivity()

    n_regions_merged = len(np.unique(merged.cell_data["RegionId"]))

    if n_regions_merged < (n_regions_1 + n_regions_2):
        return True
    else:
        return False


def add_solid_name_to_stl(filename, solid_name, file_type: str = "ascii") -> None:
    """Add the name of the solid to the STL file.

    Notes
    -----
    This method supports only a single block.

    """
    if file_type == "ascii":
        start_str = "solid"
        end_str = "endsolid"
        f = open(filename, "r")
        list_of_lines = f.readlines()
        f.close()
        list_of_lines[0] = "{0} {1}\n".format(start_str, solid_name)
        list_of_lines[-1] = "{0} {1}\n".format(end_str, solid_name)

        f = open(filename, "w")
        f.writelines(list_of_lines)
        f.close()
    # replace part name in binary file
    elif file_type == "binary":
        with open(filename, "r+b") as fid:
            fid.seek(0)  # Go to the start of the file
            string_replace = "{:<40}".format(solid_name).encode()  # Format and encode the string
            fid.write(string_replace)
        fid.close()
    return


def find_corresponding_points(
    first_surface: pv.PolyData, second_surface: pv.PolyData, distance: float = 20
) -> np.ndarray:
    """Find corresponding points between two surfaces.

    Parameters
    ----------
    first_surface : pv.PolyData
        First surface.
    second_surface : pv.PolyData
        Second surface.
    distance : float
        Approximate largest distance between two surfaces.

    Returns
    -------
    np.ndarray
        2*N array
        The first row contains node IDs of the first surface.
        The second row contains corresponding node IDs on the second surface.
        If no corresponding node is found, ``None`` is returned.

    Notes
    -----
    This method uses ray tracing.
    The two surfaces are assumed to be close and nearly parallel.
    As a result, the correspondence is not one to one. Some points might
    have no corresponding match, while others might share the same
    corresponding point.
    """
    # NOTE: Using UVCs lead to a shift in the
    # longitudinal direction from epicardium to endocardium and is thus not an option.

    # Compute normal of the first surface
    first_surface.compute_normals(inplace=True)

    points_m = first_surface.points
    normals_m = first_surface.point_data["Normals"]

    # corresponding points
    corresp_points = []
    tree_s = second_surface.find_closest_point

    # Find intersections
    for i in range(len(points_m)):
        start_point = points_m[i]
        direction = normals_m[i]
        # Cast a ray along the normal direction
        intersection, _ = second_surface.ray_trace(
            start_point - direction * distance, start_point + direction * distance
        )

        if len(intersection) == 1:
            corresp_points.append(tree_s(intersection[0]))
        elif len(intersection) > 1:
            # Take the closet intersection point
            x = np.argmin(np.linalg.norm(intersection - start_point, axis=1))
            corresp_points.append(tree_s(intersection[x]))
        else:
            corresp_points.append(None)  # fill None for no corresponding point

    return np.vstack((range(0, first_surface.n_points), corresp_points))


def generate_thickness_lines(
    surface1: pv.PolyData, surface2: pv.PolyData, corresponding_points: np.ndarray = None
) -> pv.PolyData:
    """
    Generate lines from points on surface 1 to corresponding points on surface 2.

    Parameters
    ----------
    surface1 : pv.PolyData
        First surface.
    surface2 : pv.PolyData
        Second surface.
    corresponding_points : np.ndarray, default: None
        Corresponding points array.

    Returns
    -------
    pv.PolyData
        Object contains cell data named ``thickness``.
    """
    if corresponding_points is None:
        corresponding_points = find_corresponding_points(surface1, surface2)

    points = []
    lines = []
    thickness = []
    idx = 0

    for i in range(corresponding_points.shape[1]):
        if corresponding_points[1, i] is not None:
            # make sure point pair exist
            p1 = surface1.points[corresponding_points[0, i]]
            p2 = surface2.points[corresponding_points[1, i]]
            points.append(p1)
            points.append(p2)

            # connectivity
            start_idx = idx * 2
            lines.append([2, start_idx, start_idx + 1])
            # thickness
            thickness.append(np.linalg.norm(p1 - p2))
            idx += 1

    # Create a PolyData object containing all lines
    lines = pv.PolyData(np.array(points), lines=np.hstack(lines))
    lines.cell_data["thickness"] = thickness

    return lines.copy()
