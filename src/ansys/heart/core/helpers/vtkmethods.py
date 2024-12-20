# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
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

"""Module contains methods for mesh operations related to the vtk library."""

import copy
from typing import List, Optional, Tuple, Union

from deprecated import deprecated
import numpy as np
import pyvista as pv
import vtk
from vtk.numpy_interface import dataset_adapter as dsa  # type: ignore # noqa
from vtk.util import numpy_support as VN  # type: ignore # noqa
from vtk.util.numpy_support import numpy_to_vtk  # type: ignore # noqa

from ansys.heart.core import LOG as LOGGER


# TODO: replace partially with pyvista objects for convenience.
def get_tetra_info_from_unstructgrid(
    vtk_grid: vtk.vtkUnstructuredGrid, get_all_data: bool = True, deep_copy: bool = False
) -> Tuple[np.ndarray, np.ndarray, dict, dict]:
    """Get tetrahedron nodes, connectivity and cell/point data."""
    LOGGER.warning("This method will be deprecated: can use pyvista objects instead.")
    LOGGER.debug("Extracting tetrahedron cell and point data...")
    # read nodes into numpy array
    nodes = VN.vtk_to_numpy(vtk_grid.GetPoints().GetData())
    nodes = pv.UnstructuredGrid(vtk_grid).points

    # gets number of cells
    num_cells = vtk_grid.GetNumberOfCells()

    connect1 = VN.vtk_to_numpy(vtk_grid.GetCells().GetData())
    # test vtk exists only tetra type element
    assert len(connect1) / num_cells == 5.0

    # remove ID
    idx_remove = np.arange(0, len(connect1), 5)
    tetra = np.delete(connect1, idx_remove)
    # return with  (nelem,4)
    tetra = np.reshape(tetra, (num_cells, 4))

    # get cell and node data from vtk polydata
    cell_data = np.zeros((num_cells))
    point_data = np.zeros((num_cells))

    # get number of cell data arrays:
    num_cell_arrays = vtk_grid.GetCellData().GetNumberOfArrays()
    num_point_arrays = vtk_grid.GetPointData().GetNumberOfArrays()

    # grab cell and point data if any
    cell_data = {}
    point_data = {}

    if get_all_data:
        for ii in np.arange(0, num_cell_arrays, 1):
            array_name = vtk_grid.GetCellData().GetArrayName(ii)
            cell_data[array_name] = VN.vtk_to_numpy(vtk_grid.GetCellData().GetArray(array_name))

        for ii in np.arange(0, num_point_arrays, 1):
            array_name = vtk_grid.GetPointData().GetArrayName(ii)
            point_data[array_name] = VN.vtk_to_numpy(vtk_grid.GetPointData().GetArray(array_name))

    else:
        LOGGER.debug("Not implemented reading specific cell data arrays based on name yet")

    if not deep_copy:
        return nodes, tetra, cell_data, point_data
    else:
        return (
            copy.deepcopy(nodes),
            copy.deepcopy(tetra),
            copy.deepcopy(cell_data),
            copy.deepcopy(point_data),
        )


# TODO: replace with pyvista.
@deprecated(reason="This method will be deprecated: can use pyvista methods instead.")
def add_vtk_array(
    polydata: Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid],
    data: np.array,
    name: str,
    data_type: str = "cell",
    array_type: Union[int, float] = float,
):
    """Add vtk array to vtk polydata or unstructured grid object.

    Parameters
    ----------
    polydata : Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid]
        vtk object
    data : np.array
        input data. Can be either 1d array or 2d array
    name : str
        name of data field
    data_type : str, optional
        Type of data; either "cell" or "point", by default "cell"
    array_type : Union[int, float], optional
        Type of array to add, by default float
    """
    LOGGER.warning("This method will be deprecated: can use pyvista objects instead.")
    num_dim = len(data.shape)
    num_rows = data.shape[0]

    if data_type == "point":
        num_points = polydata.GetNumberOfPoints()
    elif data_type == "cell":
        num_points = polydata.GetNumberOfCells()
    else:
        raise ValueError("Expecting 'point' or 'cell'")

    assert num_rows == num_points

    # determine type of array:
    if array_type is int:
        vtk_array = vtk.vtkIntArray()
        data = np.array(data, int)
    elif array_type is float:
        vtk_array = vtk.vtkDoubleArray()

    vtk_array.SetName(name)

    # for multi dimensional array
    if num_dim > 1:
        num_cols = data.shape[1]
        vtk_array.SetNumberOfComponents(num_rows)
        vtk_array.SetNumberOfTuples(num_cols)
        for x in zip(range(num_points), data.T):
            vtk_array.SetTuple(*x)

    elif num_dim == 1:
        vtk_array.SetNumberOfValues(num_rows)
        for x in zip(range(num_points), data.T):
            vtk_array.SetValue(*x)

    if data_type == "cell":
        polydata.GetCellData().AddArray(vtk_array)
    elif data_type == "point":
        polydata.GetPointData().AddArray(vtk_array)

    return


# TODO: replace with pyvista.
@deprecated(reason="This method will be deprecated: can use pyvista methods instead.")
def create_vtk_polydata_from_points(points: np.ndarray) -> vtk.vtkPolyData:
    """Create VTK PolyData object from set of points.

    Parameters
    ----------
    points : np.array
        Point coordinates Nx3

    Returns
    -------
    vtk.vtkPolyData
        vtkPolyData object

    Notes
    -----
    To visualize in ParaView render the points as Gaussian Points
    """
    if len(points.shape) < 2 or points.shape[1] != 3:
        raise ValueError("Expecting an array of dimension Nx3")

    vtk_poly = vtk.vtkPolyData()
    poly = dsa.WrapDataObject(vtk_poly)
    poly.Points = points
    nodeids = np.arange(0, points.shape[0], 1)

    # add point array with node ids
    add_vtk_array(poly.VTKObject, nodeids, name="node-ids", data_type="point")

    return poly.VTKObject


def compute_surface_nodal_area_pyvista(surface: pv.PolyData) -> np.ndarray:
    """Compute an average nodal area by summing surface areas of connected elements.

    Parameters
    ----------
    vtk_surface : vtk.vtkPolyData
        Vtk object describing the object

    Returns
    -------
    np.array
        Numpy array with nodal areas of length number of points

    Notes
    -----
    Adds the partial areas of connected elements/cells to each node.

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


# TODO: replace with pyvista.
@deprecated(reason="This method will be deprecated: can use pyvista methods instead.")
def add_normals_to_polydata(
    vtk_polydata: vtk.vtkPolyData, return_normals: bool = False
) -> Union[vtk.vtkPolyData, Optional[Tuple[np.ndarray, np.ndarray]]]:
    """Add normals to the vtk.vtkPolyData object.

    Parameters
    ----------
    vtk_polydata : vtk.vtkPolyData
        Input surface.
    return_normals : bool, optional
        Return the cell and point normals as numpy arrays, by default False.

    Returns
    -------
    vtk_polydata : vtk.vtkPolyData
        Vtk surface with cell and point normals added.
    (cell_normals, point_normals) : (np.ndarray, np.ndarray), optional
        Cell normals and point normals, only provided if return_normals=True
    """
    # compute normals
    normal_filter = vtk.vtkPolyDataNormals()
    normal_filter.SetInputData(vtk_polydata)
    normal_filter.ComputeCellNormalsOn()
    normal_filter.ComputePointNormalsOn()
    # normal_filter.AutoOrientNormalsOn()
    # normal_filter.ConsistencyOff()
    # normal_filter.NonManifoldTraversalOff()
    # normal_filter.SetSplitting(0)

    normal_filter.SplittingOff()
    normal_filter.SetFeatureAngle(30)

    normal_filter.Update()

    if return_normals:
        normal_filter_dsa = dsa.WrapDataObject(normal_filter.GetOutput())
        return np.array(normal_filter_dsa.CellData["Normals"]), np.array(
            normal_filter_dsa.PointData["Normals"]
        )
    else:
        return normal_filter.GetOutput()


def extrude_polydata(
    vtk_surface: vtk.vtkPolyData,
    extrude_by: float = 1,
    extrude_direction: np.array = np.empty(0),
) -> vtk.vtkPolyData:
    """Extrude a given polydata surface in a given direction.

    Parameters
    ----------
    vtk_surface : vtk.vtkPolyData
        Surface to extrude
    extrude_by : float, optional
        Extrude by this much, by default 1
    extrude_direction : np.array, optional
        Direction of extrusion, should have three components if not specified
        extrudes in normal direction

    Returns
    -------
    vtk.vtkPolyData
        Extruded vtkPolyData object
    """
    extrude_normal = False
    if len(extrude_direction) == 0:
        extrude_normal = True

    vtk_surface = add_normals_to_polydata(vtk_surface)

    extrude = vtk.vtkLinearExtrusionFilter()
    extrude.CappingOn()

    if extrude_normal:
        extrude.SetExtrusionTypeToNormalExtrusion()
    else:
        extrude.SetExtrusionTypeToVectorExtrusion()
        extrude.SetVector(extrude_direction[0], extrude_direction[1], extrude_direction[2])

    extrude.SetInputData(vtk_surface)
    extrude.SetScaleFactor(extrude_by)
    extrude.Update()
    extruded_polydata = extrude.GetOutput()

    return pv.PolyData(extruded_polydata)


# TODO: replace with pyvista.
@deprecated(reason="This method will be deprecated: can use pyvista methods instead.")
def create_vtk_surface_triangles(
    points: np.ndarray, triangles: np.ndarray, clean=True
) -> vtk.vtkPolyData:
    """Create vtkPolyData object from array of points and array of triangles.

    Parameters
    ----------
    points : np.array
        Nx3 array of point coordinates
    triangles : np.array
        Mx3 array of triangle definitions
    clean : Boolean, True by default
        use vtkCleanPolyData Filter to remove unused nodes, etc.
        But may have unexpected behavior...

    Returns
    -------
    vtk.vtkPolyData
        VTK Object PolyData object describing the surface
    """
    LOGGER.warning("This method will be deprecated: can initialize pyvista objects directly.")
    num_points = points.shape[0]
    points_vtk = vtk.vtkPoints()
    points_vtk.SetNumberOfPoints(num_points)
    points_vtk.SetData(numpy_to_vtk(np.asarray(points, order="C", dtype=float), deep=1))

    triangles_vtk = vtk.vtkCellArray()
    for tri in triangles:
        triangle_vtk = vtk.vtkTriangle()
        triangle_vtk.GetPointIds().SetId(0, tri[0])
        triangle_vtk.GetPointIds().SetId(1, tri[1])
        triangle_vtk.GetPointIds().SetId(2, tri[2])
        triangles_vtk.InsertNextCell(triangle_vtk)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points_vtk)
    polydata.SetPolys(triangles_vtk)
    polydata.Modified()
    # polydata.Update()

    if clean:
        # clean polydata
        clean_filter = vtk.vtkCleanPolyData()
        clean_filter.SetInputData(polydata)
        clean_filter.PointMergingOn()
        clean_filter.Update()

        polydata_cleaned = clean_filter.GetOutput()
        return polydata_cleaned

    else:
        return polydata


# TODO: replace with pyvista.
@deprecated(reason="This method will be deprecated: can use pyvista methods instead.")
def cell_ids_inside_enclosed_surface(
    vtk_source: vtk.vtkUnstructuredGrid, vtk_surface: vtk.vtkPolyData
) -> np.ndarray:
    """Tag any cells that are inside an enclosed surface.

    Parameters
    ----------
    vtk_source : vtk.vtkUnstructuredGrid
        Source VTK object of which to check the whether the cells are inside/outside
        the specified surface
    vtk_surface : vtk.vtkPolyData
        Enclosed surface

    Returns
    -------
    vtk.vtkUnstructuredGrid
        VTK object with additional cell data indicating whether
        the cell is in/outside the provided surface
    """
    vtk_surface = add_normals_to_polydata(vtk_surface)
    points = pv.UnstructuredGrid(vtk_source).points
    tetra = pv.UnstructuredGrid(vtk_source).cells_dict[pv.CellType.TETRA]

    centroids = np.mean(points[tetra, :], axis=1)

    vtk_centroids = create_vtk_polydata_from_points(centroids)

    select = vtk.vtkSelectEnclosedPoints()
    select.SetInputData(vtk_centroids)
    select.SetSurfaceData(vtk_surface)
    select.CheckSurfaceOn()
    select.SetTolerance(1e-9)  # fraction of diagonal of bounding box!
    select.Update()

    output = select.GetOutput()

    output_dsa = dsa.WrapDataObject(output)

    cell_ids_inside = np.where(output_dsa.PointData["SelectedPoints"] == 1)[0]

    return cell_ids_inside


def find_cells_close_to_nodes(
    mesh: pv.UnstructuredGrid, node_ids: list[int], radius: float = 2
) -> np.ndarray:
    """Find cell IDs close to nodes.

    Parameters
    ----------
    mesh : pv.UnstructuredGrid
        target mesh
    node_ids : list[int]
        node IDs
    radius : float, optional
        influence radius, by default 2

    Returns
    -------
    np.ndarray
        cell IDs
    """
    # Get coordinates of the given node IDs
    points = mesh.points[node_ids]

    # Create a list to store cells within the sphere radius
    cells_within_sphere = []

    # Iterate through each point and find cells within the sphere
    for point in points:
        # Create a sphere at the given point
        sphere = pv.Sphere(radius=radius, center=point)

        # Use boolean intersection to find cells that intersect with the sphere
        selection = mesh.select_enclosed_points(sphere, tolerance=0.0)

        # Get the indices of the cells
        selected_points = selection.point_data["SelectedPoints"].nonzero()[0]
        selected_cells = mesh.extract_points(selected_points).cell_data["vtkOriginalCellIds"]

        # Store unique cell indices
        cells_within_sphere.extend(selected_cells)

    # Return unique cell indices
    return np.unique(cells_within_sphere)


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
) -> dict:
    """Get the closed/open boundary edge loops of a surface mesh.

    Parameters
    ----------
    surface : pv.PolyData
        Surface mesh to check for boundary edges
    remove_open_edge_loops : bool, optional
        Removes open edge loops from the return dictionary, by default True

    Returns
    -------
    dict
        dictionary with the edges that make up the open/closed loop
    """
    # NOTE: Perhaps more consistent to return a pyvista polydata.

    # add cell and point ids to keep track of ids.
    surface1 = copy.deepcopy(surface)
    surface1.cell_data["original-cell-ids"] = np.arange(0, surface1.n_cells)
    surface1.point_data["original-point-ids"] = np.arange(0, surface1.n_points)

    # get boundary edges separated by connectivity
    edges_block = get_boundary_edges(surface1)

    # lines formed with original point ids
    edge_groups = {
        k: edges.point_data["original-point-ids"][edges.cells_dict[3]]
        for k, edges in enumerate(edges_block)
    }

    # check if it is a closed or open edge-loop, remove open ones.
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
        LOGGER.warning(f"Removed {num_open_edge_groups} edge groups")
        return closed_edge_groups
    else:
        if return_types:
            return edge_groups, group_types
        else:
            return edge_groups


def get_patches_delaunay(surface: pv.PolyData, closed_only: bool = True) -> List[pv.PolyData]:
    """Patch boundary edges with a delaunay algorithm.

    Parameters
    ----------
    surface : pv.PolyData
        Surface with boundary edges for which to find patches.
    closed_only : bool
        Flag indicating whether to return patches for closed loops of boundary edges,
        by default True

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


def get_patches_with_centroid(surface: pv.PolyData, closed_only: bool = True) -> List[pv.PolyData]:
    """Patch boundary edges with a custom algorithm using a central node.

    Parameters
    ----------
    surface : pv.PolyData
        Surface with boundary edges for which to find patches.
    closed_only : bool
        Flag indicating whether to return patches for closed loops of boundary edges,
        by default True

    Notes
    -----
    Edges need to be sorted properly for this method to return sensible patches.

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
    mesh1: Union[pv.PolyData, pv.UnstructuredGrid], mesh2: Union[pv.PolyData, pv.UnstructuredGrid]
):
    """Check whether two PolyData or UnstructuredGrids are connected.

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
        LOGGER.warning("RegionId was not present in the cell data")

    merged = merged.connectivity()

    n_regions_merged = len(np.unique(merged.cell_data["RegionId"]))

    if n_regions_merged < (n_regions_1 + n_regions_2):
        return True
    else:
        return False


def add_solid_name_to_stl(filename, solid_name, file_type: str = "ascii") -> None:
    """Add name of solid to stl file.

    Notes
    -----
    Supports only single block.

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


if __name__ == "__main__":
    print()
