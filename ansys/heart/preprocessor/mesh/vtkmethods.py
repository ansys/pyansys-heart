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
from typing import Optional, Tuple, Union

from ansys.heart.core import LOG as LOGGER
from ansys.heart.preprocessor.mesh.fluenthdf5 import add_solid_name_to_stl
import numpy as np
import pyvista as pv
import vtk
from vtk.numpy_interface import dataset_adapter as dsa  # type: ignore # noqa
from vtk.util import numpy_support as VN  # type: ignore # noqa
from vtk.util.numpy_support import numpy_to_vtk  # type: ignore # noqa


def write_vtkdata_to_vtkfile(vtk_data: Union[vtk.vtkUnstructuredGrid, vtk.vtkPolyData], fname: str):
    """Write a vtk unstructured grid object to vtk file."""
    DeprecationWarning("This method is deprecated: use pyvista objects instead.")
    writer = vtk.vtkDataSetWriter()
    writer.SetFileName(fname)
    writer.SetInputData(vtk_data)
    writer.SetFileTypeToBinary()
    writer.Write()
    return


def vtk_surface_to_stl(
    vtk_data: Union[vtk.vtkUnstructuredGrid, vtk.vtkPolyData], filename: str, solid_name: str = None
) -> None:
    """Write stl from vtk surface mesh (polydata)."""
    DeprecationWarning("This method is deprecated: use pyvista objects instead.")

    vtk_data_to_write = vtk_data

    if type(vtk_data) != type(vtk.vtkPolyData()):
        # convert to polydata
        surface_filter = vtk.vtkDataSetSurfaceFilter()
        surface_filter.SetInputData(vtk_data)
        surface_filter.Update()
        vtk_data_to_write = surface_filter.GetOutput()

    writer = vtk.vtkSTLWriter()
    writer.SetFileName(filename)
    writer.SetInputData(vtk_data_to_write)
    writer.SetFileTypeToBinary()
    writer.Write()

    if solid_name:
        add_solid_name_to_stl(filename, solid_name, "binary")

    return


def get_tetra_info_from_unstructgrid(
    vtk_grid: vtk.vtkUnstructuredGrid, get_all_data: bool = True, deep_copy: bool = False
) -> Tuple[np.ndarray, np.ndarray, dict, dict]:
    """Get tetrahedron nodes, connectivity and cell/point data."""
    LOGGER.warning("This method will be deprecated: can use pyvista objects instead.")
    LOGGER.debug("Extracting tetrahedron cell and point data...")
    # read nodes into numpy array
    nodes = VN.vtk_to_numpy(vtk_grid.GetPoints().GetData())

    # gets number of cells
    num_cells = vtk_grid.GetNumberOfCells()
    num_points = vtk_grid.GetNumberOfPoints()

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


def get_tri_info_from_polydata(
    vtk_polydata: vtk.vtkPolyData, get_all_data: bool = True, deep_copy: bool = False
) -> Tuple[np.ndarray, np.ndarray, dict, dict]:
    """Get connectivity, celldata and point data info from polydata object.

    Notes
    -----
    Assumes triangular elements
    """
    LOGGER.warning("This method will be deprecated: can use pyvista objects instead.")
    # logger.debug("Extracting triangle cell and point data...")

    vtk_polydata_obj = dsa.WrapDataObject(vtk_polydata)

    nodes = np.array(vtk_polydata_obj.Points)
    polys = np.array(vtk_polydata_obj.Polygons)

    if not np.all(polys[np.arange(0, len(polys), 4)]):
        raise TypeError("Some polygons are not triangular")

    tris = polys.reshape(-1, 4)
    tris = np.delete(tris, 0, axis=1)

    # store cell/point data in dictionary
    cell_data = {}
    point_data = {}

    for key in vtk_polydata_obj.CellData.keys():
        cell_data[key] = np.array(vtk_polydata_obj.CellData[key])

    for key in vtk_polydata_obj.PointData.keys():
        point_data[key] = np.array(vtk_polydata_obj.PointData[key])

    # try to add global ids
    try:
        global_ids = VN.vtk_to_numpy(vtk_polydata.GetPointData().GetGlobalIds())
        point_data["global_ids"] = global_ids
    except:
        LOGGER.debug("Global Ids were not added to point data...")

    if not deep_copy:
        return nodes, tris, cell_data, point_data
    else:
        return (
            copy.deepcopy(nodes),
            copy.deepcopy(tris),
            copy.deepcopy(cell_data),
            copy.deepcopy(point_data),
        )


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
    LOGGER.warning("This method will be deprecated: can use pyvista objects instead.")
    if len(points.shape) < 2 or points.shape[1] != 3:
        raise ValueError("Expecting an array of dimension Nx3")

    vtk_poly = vtk.vtkPolyData()
    poly = dsa.WrapDataObject(vtk_poly)
    poly.Points = points
    nodeids = np.arange(0, points.shape[0], 1)

    # add point array with node ids
    add_vtk_array(poly.VTKObject, nodeids, name="node-ids", data_type="point")

    return poly.VTKObject


def remove_duplicate_nodes(
    nodes: np.ndarray, elements: np.ndarray, tolerance: float = 1e-7
) -> Tuple[np.ndarray, np.ndarray]:
    """Find and removes duplicate nodes and remaps element definitions.

    Parameters
    ----------
    nodes : np.ndarray
        Array with nodal coordinates
    elements : np.ndarray
        Array with element definition
    tolerance: float
        Tolerance - the same for each coordinate
    """
    nodes_rounded = np.array(np.round(nodes * 1 / tolerance), dtype=int)
    unique_nodes, indices, inverse_indices, counts = np.unique(
        nodes_rounded,
        axis=0,
        return_index=True,
        return_inverse=True,
        return_counts=True,
    )
    unique_nodes = nodes[indices, :]
    elements = inverse_indices[elements]

    return unique_nodes, elements


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
    LOGGER.warning("This method will be deprecated: can use pyvista .compute_normals() instead.")
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

    return extruded_polydata


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
    LOGGER.warning("This method will be deprecated: can use pyvista methods instead.")
    vtk_surface = add_normals_to_polydata(vtk_surface)
    points, tetra, _, _ = get_tetra_info_from_unstructgrid(vtk_source)

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


def get_connected_regions(
    nodes: np.ndarray, triangles: np.ndarray, return_vtk_object: bool = False
) -> np.ndarray:
    """Find the connected regions.

    Parameters
    ----------
    nodes : np.ndarray
        NumNodes x 3 array with point coordinates
    triangles : np.ndarray
        NumTriangles x 3 array with triangle definitions
    return_vtk_object : bool, optional
        Flag indicating whether to return the vtk (surface) object, by default False

    Returns
    -------
    np.array
        Array with region ids
    vtk.vtkPolyData
        VTK Object with region ids
    """
    LOGGER.warning("This method will be deprecated: can use pyvista objects instead.")
    vtk_surface = create_vtk_surface_triangles(nodes, triangles)

    # connectivity filter to extract all connected regions
    connectivity0 = vtk.vtkPolyDataConnectivityFilter()
    connectivity0.SetExtractionModeToAllRegions()
    connectivity0.SetColorRegions(1)
    connectivity0.SetInputData(vtk_surface)
    connectivity0.Update()
    vtk_surface_with_regions = connectivity0.GetOutput()

    vtk_surface_dsa = dsa.WrapDataObject(vtk_surface_with_regions)
    region_ids = vtk_surface_dsa.PointData["RegionId"]

    # add separate array to facilitate point2cell filter. Somehow region ids are protected from
    # this filter
    add_vtk_array(
        vtk_surface_with_regions, region_ids + 1, "regions", data_type="point", array_type=int
    )
    # map point data to cell data
    point2cell = vtk.vtkPointDataToCellData()
    point2cell.SetInputData(vtk_surface_with_regions)
    point2cell.Update()
    vtk_surface_with_regions = point2cell.GetOutput()

    # get cell data from triangular polydata
    cell_data = get_tri_info_from_polydata(vtk_surface_with_regions)[2]

    region_ids = cell_data["regions"]
    # cast to ints
    region_ids = np.array(region_ids, dtype=int)
    if return_vtk_object:
        return region_ids, vtk_surface_with_regions
    else:
        return region_ids


def vtk_cutter(vtk_polydata: vtk.vtkPolyData, cut_plane) -> vtk.vtkPolyData:
    """
    Cut a vtk polydata by a plane.

    Parameters
    ----------
    vtk_polydata: vtk polydata
    cut_plane: dictionary contains key: 'center' and 'normal'

    Returns
    -------
    vtkpolydata
    """
    LOGGER.warning("This method will be deprecated: can use pyvista objects instead.")
    # create a plane to cut
    plane = vtk.vtkPlane()
    plane.SetOrigin(cut_plane["center"][0], cut_plane["center"][1], cut_plane["center"][2])
    plane.SetNormal(cut_plane["normal"][0], cut_plane["normal"][1], cut_plane["normal"][2])

    # create cutter
    cutter = vtk.vtkCutter()
    # cutter.SetNumberOfContours(20) #how to control number of points?
    cutter.SetCutFunction(plane)
    cutter.SetInputData(vtk_polydata)
    cutter.Update()

    return cutter.GetOutput()


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


if __name__ == "__main__":
    print()
