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

from ansys.heart.core import LOG as LOGGER
from ansys.heart.preprocessor.mesh.fluenthdf5 import add_solid_name_to_stl
import numpy as np
import pyvista as pv
import vtk
from vtk.numpy_interface import dataset_adapter as dsa  # type: ignore # noqa
from vtk.util import numpy_support as VN  # type: ignore # noqa
from vtk.util.numpy_support import numpy_to_vtk  # type: ignore # noqa


def read_vtk_unstructuredgrid_file(path_to_vtk: str) -> vtk.vtkUnstructuredGrid:
    """Read vtk unstructured grid file."""
    DeprecationWarning("This method is deprecated: use pyvista objects instead.")
    reader = vtk.vtkUnstructuredGridReader()  # noqa
    reader.SetFileName(path_to_vtk)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.ReadAllFieldsOn()
    reader.Update()
    return reader.GetOutput()


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


def vtk_surface_filter(
    vtk_grid: vtk.vtkUnstructuredGrid, keep_global_ids: bool = False
) -> vtk.vtkPolyData:
    """Extract surface from a vtk object (polydata or unstructured grid)."""
    DeprecationWarning("This method is deprecated: use pyvista objects instead.")
    LOGGER.debug("Extracting surface from vtk unstructured grid...")
    # make sure global id will be kept
    if keep_global_ids:
        with_id = vtk.vtkGenerateGlobalIds()  # noqa
        with_id.SetInputData(vtk_grid)
        with_id.Update()
        # extract surface
        surface = vtk.vtkDataSetSurfaceFilter()  # noqa
        surface.SetInputData(with_id.GetOutput())
        surface.Update()
    else:
        surface = vtk.vtkDataSetSurfaceFilter()  # noqa
        surface.SetInputData(vtk_grid)
        surface.Update()

    return surface.GetOutput()


def get_info_from_vtk(
    vtk_grid: Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid]
) -> Tuple[np.ndarray, np.ndarray, dict, dict]:
    """Use numpy support to get points, cell connectivity, cell data, point data from vtk object.

    Parameters
    ----------
    vtk_grid : Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid]
        Vtk object from which to extract the information

    Returns
    -------
    List
        List of points, elements, cell data, and point data

    Raises
    ------
    ValueError
        _description_
    """
    LOGGER.warning("This method will be deprecated: can use pyvista objects instead.")
    vtk_obj = dsa.WrapDataObject(vtk_grid)

    num_cells = vtk_obj.GetNumberOfCells()
    points = vtk_obj.Points  # nodes
    cells_conn = vtk_obj.Cells  # cell connectivity matrix
    cell_types = vtk_obj.CellTypes

    cell_data = {}
    point_data = {}

    for key in vtk_obj.CellData.keys():
        cell_data[key] = vtk_obj.CellData[key]

    for key in vtk_obj.PointData.keys():
        cell_data[key] = vtk_obj.PointData[key]

    elements = {}
    if np.all(cell_types == 10):
        # 10 = tetrahedron
        elements = cells_conn.reshape(num_cells, 4 + 1)
        elements = elements[:, 1:]  # remove one column
    elif np.all(cell_types == 1):
        # 1 = vertex
        elements = cells_conn.reshape(num_cells, 1 + 1)
        elements = elements[:, 1:]  # remove one column
    else:
        raise ValueError("Not supporting anything other than tetrahedrons and vertices")
    # np.savetxt("source_points.csv", points, delimiter=",")
    return points, elements, cell_data, point_data


def vtk_map_discrete_cell_data(
    vtk_object_source: Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid],
    vtk_object_target: Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid],
    data_name: str,
) -> Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid]:
    """
    Map discrete values from a source to a target.

    Notes
    -----
    Uses linear interpolation with
    1 closest point. Note that this computes the centroids of the cells first, creates
    a new vtk poly data object as point cloud and uses that for interpolating the target
    cell data field.

    Parameters
    ----------
    vtk_object_source : Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid]
        Source vtk object
    vtk_object_target : Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid]
        Target vtk object
    """
    DeprecationWarning("This method is deprecated: use pyvista objects instead.")
    source_points, source_tets, source_cdata, source_pdata = get_info_from_vtk(vtk_object_source)
    target_points, target_tets, target_cdata, target_pdata = get_info_from_vtk(vtk_object_target)

    centroids_source = np.mean(source_points[source_tets, :], axis=1)
    centroids_target = np.mean(target_points[target_tets, :], axis=1)

    tags_source = source_cdata[data_name]

    vtk_poly_source = vtk.vtkPolyData()

    # use numpy interface for more compact syntax
    poly_source = dsa.WrapDataObject(vtk_poly_source)
    poly_source.Points = centroids_source

    # add tags array to poly data
    add_vtk_array(poly_source.VTKObject, tags_source, name=data_name, data_type="point")

    # set up target poly data
    vtk_poly_target = vtk.vtkPolyData()
    poly_target = dsa.WrapDataObject(vtk_poly_target)
    poly_target.Points = centroids_target

    # interpolate tags from source(input) to target(source)
    # interpolator linear kernel with 1 closest point
    linearKernel = vtk.vtkLinearKernel()
    linearKernel.SetKernelFootprintToNClosest()
    linearKernel.SetNumberOfPoints(1)

    interpolator = vtk.vtkPointInterpolator()
    # interpolates from input > source?
    interpolator.SetInputData(poly_target.VTKObject)
    interpolator.SetSourceData(poly_source.VTKObject)  # interpolates from input > source?
    # interpolator.SetKernel(voronoiKernel)
    interpolator.SetKernel(linearKernel)
    interpolator.Update()

    poly_target = dsa.WrapDataObject(interpolator.GetOutput())

    tags_target = poly_target.PointData["tags"]

    # adds the target tags as a cell array to the target
    add_vtk_array(vtk_object_target, tags_target, name="tags", data_type="cell")

    return vtk_object_target


def vtk_map_continuous_data(
    source: Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid],
    target: Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid],
    normalize_vectors: bool = True,
    array_names_to_include: list = [],
) -> Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid]:
    """Map cell and point data from source to target.

    Parameters
    ----------
    input : Union[vtk.PolyData, vtk.UnstructuredGrid]
        Input which to use as reference
    source : Union[vtk.PolyData, vtk.UnstructuredGrid]
        Target object onto which to interpolate data
    array_names_to_include : list
        List of array names to include for interpolation.
        If empty all cell and point arrays will be interpolated

    Notes
    -----
    Makes use of VoronoiKernel and mapping cell to point data
    consequently point data is mapped back to the cell.

    Modifies the underlying data of the target vtk object and overwrites if
    a data field with the same name is already present.
    """
    # NOTE: could use AddExcludeArray to exclude specific arrays from the interpolation
    LOGGER.warning("This method will be deprecated: can use pyvista objects instead.")

    if array_names_to_include == []:
        include_all = True
    else:
        include_all = False

    # get info on the data arrays
    source_obj = dsa.WrapDataObject(source)
    target_obj = dsa.WrapDataObject(target)

    cell_array_names_source = source_obj.CellData.keys()
    point_array_names_source = source_obj.PointData.keys()

    cell_array_names_target = target_obj.CellData.keys()
    point_array_names_target = target_obj.PointData.keys()

    # convert to point data
    cell2point = vtk.vtkCellDataToPointData()
    cell2point.SetInputData(source)
    cell2point.PassCellDataOff()
    cell2point.SetContributingCellOption(0)  # try 0: All; 1: Patch; 2: DataSetMax
    cell2point.Update()
    source_interpolator = cell2point.GetOutput()

    # interpolator: voronoi kernel:
    voronoiKernel = vtk.vtkVoronoiKernel()

    interpolator = vtk.vtkPointInterpolator()
    interpolator.SetInputData(
        target
    )  # input data is data on which to interpolate. NOTE Paraview naming seems wrong
    interpolator.SetSourceData(source_interpolator)
    interpolator.SetKernel(voronoiKernel)
    # interpolator.SetKernel( linearKernel )

    # add list of excluded array names:
    # does not distinguish between cell or point data
    if not include_all:
        excludes = list(
            set(cell_array_names_source + point_array_names_source) - set(array_names_to_include)
        )
        for exclude in excludes:
            interpolator.AddExcludedArray(exclude)
        LOGGER.debug("Excluding %d array names" % len(excludes))

    interpolator.Update()

    # convert all data to cell data: NOTE here an interpolation error occurs
    # again
    point2cell = vtk.vtkPointDataToCellData()
    point2cell.SetInputData(interpolator.GetOutput())
    point2cell.PassPointDataOn()
    point2cell.Update()
    target_updated = point2cell.GetOutput()

    # clean up vtk arrays: compare with datafields originally
    # present in the source/target
    ii = 0
    while True:
        name = target_updated.GetCellData().GetArrayName(ii)
        if name is None:
            break

        if name not in cell_array_names_source and name not in cell_array_names_target:
            LOGGER.debug("Removing cell data..." + name)
            vtk_remove_arrays(target_updated, array_name=name, data_type="cell_data")

        else:
            ii = ii + 1

    ii = 0
    while True:
        name = target_updated.GetPointData().GetArrayName(ii)
        if name is None:
            break

        if name not in point_array_names_source and name not in point_array_names_target:
            LOGGER.debug("Removing point data..." + name)
            vtk_remove_arrays(target_updated, array_name=name, data_type="point_data")
        else:
            ii = ii + 1

    # NOTE 2: normalizes vector data
    target_updated1 = dsa.WrapDataObject(target_updated)

    # normalize cell data

    # NOTE "fiber" data is normalized and added correctly by "add_vtk_array"
    # however "sheet" data somehow introduces a bug in the vtk file and
    # results in a file which cannot be read with Paraview.

    normalize_vectors = False
    LOGGER.warning("Normalization of vectors is buggy and turned off")

    if normalize_vectors:
        for key in target_updated1.CellData.keys():
            data = target_updated1.CellData[key]
            num_dim = len(data.shape)
            # if vector data normalize
            if num_dim > 1:
                normalize = True
            else:
                normalize = False

            if normalize:
                LOGGER.debug("Normalizing data: " + key)
                norm = np.linalg.norm(data, axis=1)
                data = data / norm[:, None]
                data = np.array(data)
                add_vtk_array(target_updated1.VTKObject, data, key, "cell", float)

    LOGGER.warning("Removed returning cell / point data")

    return target_updated1.VTKObject


def vtk_remove_arrays(
    vtk_grid: Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid],
    array_name: str = "",
    data_type: str = "cell_data",
    remove_all: bool = False,
    except_array_names: List[str] = [],
) -> Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid]:
    """Remove all or specific data arrays from vtk object."""
    LOGGER.warning("This method will be deprecated: can use pyvista objects instead.")
    if data_type not in ["cell_data", "point_data", "both"]:
        raise ValueError("Data type not valid")

    if data_type == "both" and len(except_array_names) > 0:
        raise ValueError("Please specify either point_data or cell_data when using exception list")

    remove_cell_data = False
    remove_point_data = False

    if data_type in ["cell_data", "both"]:
        remove_cell_data = True
    if data_type in ["point_data", "both"]:
        remove_point_data = True

    num_cell_data = vtk_grid.GetCellData().GetNumberOfArrays()
    num_point_data = vtk_grid.GetPointData().GetNumberOfArrays()

    expected_number_of_arrays = len(except_array_names)
    if remove_all:
        if remove_cell_data:
            idx = 0
            while num_cell_data > expected_number_of_arrays:
                name_array_remove = vtk_grid.GetCellData().GetArrayName(idx)
                if name_array_remove in except_array_names:
                    idx = idx + 1
                    continue
                vtk_grid.GetCellData().RemoveArray(idx)
                num_cell_data = vtk_grid.GetCellData().GetNumberOfArrays()

        if remove_point_data:
            idx = 0
            while num_point_data > expected_number_of_arrays:
                name_array_remove = vtk_grid.GetPointData().GetArrayName(idx)
                if name_array_remove in except_array_names:
                    idx = idx + 1
                    continue
                vtk_grid.GetPointData().RemoveArray(idx)
                num_point_data = vtk_grid.GetPointData().GetNumberOfArrays()

    else:
        if remove_cell_data:
            for ii in range(num_cell_data):
                if array_name == vtk_grid.GetCellData().GetArrayName(ii):
                    vtk_grid.GetCellData().RemoveArray(ii)

        elif remove_point_data:
            for ii in range(num_point_data):
                if array_name == vtk_grid.GetPointData().GetArrayName(ii):
                    vtk_grid.GetPointData().RemoveArray(ii)

    return vtk_grid


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
    tris = np.reshape(surface.faces, (surface.n_faces, 4))[:, 1:]

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

    return pv.PolyData(extruded_polydata)


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


def smooth_polydata(vtk_polydata: vtk.vtkPolyData) -> vtk.vtkPolyData:
    """Use Laplacian smoothing to smooth the vtk polydata object."""
    LOGGER.warning("This method will be deprecated: can use pyvista objects instead.")
    smooth_filter = vtk.vtkSmoothPolyDataFilter()
    smooth_filter.SetInputData(vtk_polydata)
    smooth_filter.SetNumberOfIterations(15)
    smooth_filter.SetRelaxationFactor(0.01)
    smooth_filter.SetConvergence(0.0)
    smooth_filter.FeatureEdgeSmoothingOn()
    # smooth_filter.FeatureEdgeSmoothingOff() # smooths feature edges
    smooth_filter.BoundarySmoothingOn()
    smooth_filter.Update()

    # Update normals on newly smoothed polydata
    normal_gen = vtk.vtkPolyDataNormals()
    normal_gen.SetInputData(smooth_filter.GetOutput())
    normal_gen.ComputePointNormalsOn()
    normal_gen.ComputeCellNormalsOn()
    normal_gen.Update()
    vtk_polydata_smooth = normal_gen.GetOutput()
    return vtk_polydata_smooth


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

    num_open_edge_groups = len(edge_groups) != len(closed_edge_groups.keys())

    if remove_open_edge_loops:
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
    edges_block = pv.MultiBlock()
    patches = []
    for edge_group in edge_groups.values():
        centroid = np.mean(surface1.points[np.unique(edge_group), :], axis=0)

        edge_points = pv.PolyData(surface1.points[np.unique(edge_group), :])
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


if __name__ == "__main__":
    print()
