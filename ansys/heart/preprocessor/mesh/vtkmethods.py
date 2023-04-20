"""Module contains methods for mesh operations related to the vtk library."""
import copy
import os
from pathlib import Path
from typing import List, Optional, Tuple, Union

from ansys.heart.custom_logging import LOGGER
from ansys.heart.preprocessor.mesh.misc import add_solid_name_to_stl
import meshio
import numpy as np
import pyvista
import vtk
from vtk.numpy_interface import dataset_adapter as dsa  # type: ignore # noqa
from vtk.util import numpy_support as VN  # type: ignore # noqa
from vtk.util.numpy_support import numpy_to_vtk  # type: ignore # noqa


def read_ensight_file(path_to_ensight: str) -> vtk.vtkUnstructuredGrid:
    """Read ensight file."""
    file_path = os.path.dirname(path_to_ensight)
    filename = os.path.basename(path_to_ensight)

    ens = vtk.vtkGenericEnSightReader()
    ens.SetByteOrderToLittleEndian()
    ens.SetFilePath(file_path)
    ens.SetCaseFileName(filename)
    ens.ReadAllVariablesOn()
    ens.Update()
    return ens.GetOutput().GetBlock(0)


def read_vtk_unstructuredgrid_file(path_to_vtk: str) -> vtk.vtkUnstructuredGrid:
    """Read vtk unstructured grid file."""
    reader = vtk.vtkUnstructuredGridReader()  # noqa
    reader.SetFileName(path_to_vtk)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.ReadAllFieldsOn()
    reader.Update()
    return reader.GetOutput()


def read_vtk_polydata_file(path_to_vtk: str) -> vtk.vtkPolyData:
    """Read vtk PolyData file."""
    reader = vtk.vtkUnstructuredGridReader()  # noqa
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(path_to_vtk)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.ReadAllFieldsOn()
    reader.Update()
    return reader.GetOutput()


def vtk_read_mesh_file(path_to_mesh: str) -> vtk.vtkUnstructuredGrid():
    """Read either ensight format or vtk format into vtk polydata."""
    mesh_extension = Path(path_to_mesh).suffix

    # reads .case (ensight):
    if mesh_extension == ".case":
        LOGGER.debug("Reading ensight file...{0}".format(path_to_mesh))
        vtk_data = read_ensight_file(path_to_mesh)
    # reads .vtk
    elif mesh_extension == ".vtk":
        LOGGER.debug("Reading vtk file...{0}".format(path_to_mesh))
        vtk_data = read_vtk_unstructuredgrid_file(path_to_mesh)

    return vtk_data


def write_vtkdata_to_vtkfile(vtk_data: Union[vtk.vtkUnstructuredGrid, vtk.vtkPolyData], fname: str):
    """Write a vtk unstructured grid object to vtk file."""
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


def get_vtk_points(vtk_object: Union[vtk.vtkUnstructuredGrid, vtk.vtkPolyData]) -> np.ndarray:
    """Return the points of a vtk unstructured grid or polydata object."""
    return VN.vtk_to_numpy(vtk_object.GetPoints().GetData())


def get_tetra_info_from_unstructgrid(
    vtk_grid: vtk.vtkUnstructuredGrid, get_all_data: bool = True, deep_copy: bool = False
) -> Tuple[np.ndarray, np.ndarray, dict, dict]:
    """Get tetrahedron nodes, connectivity and cell/point data."""
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


def threshold_vtk_data(
    vtk_obj: Union[vtk.vtkUnstructuredGrid, vtk.vtkPolyData],
    lower_limit: Union[float, int],
    upper_limit: Union[float, int],
    data_name: str,
    epsilon: float = 1e-3,
    data_type: str = "CellData",
) -> Tuple[Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid], np.ndarray]:
    """Use the vtk thresholding filter to extract a part of the model.

    Parameters
    ----------
    vtk_obj : Union[vtk.vtkUnstructuredGrid, vtk.vtkPolyData]
        Vtk object of the model
    lower_limit : Union[float, int]
        Lower limit
    upper_limit : Union[float, int]
        Upper limit
    data_name : str
        Name of the cell data field to processes
    epsilon : _type_, optional
        Allowed tolerance for filter, by default 1e-3
    data_type: str, optional
        Type of data to filter. Either "CellData" or "PointsData"
    """
    if data_type not in ["CellData", "PointData"]:
        raise ValueError("Please specify either 'CellData' or 'PointData'")

    with_id = vtk.vtkGenerateGlobalIds()  # noqa
    with_id.SetInputData(vtk_obj)
    with_id.Update()

    threshold = vtk.vtkThreshold()  # noqa
    threshold.SetInputData(with_id.GetOutput())
    if data_type == "CellData":
        threshold.SetInputArrayToProcess(
            0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, data_name  # noqa
        )
    elif data_type == "PointData":
        threshold.SetInputArrayToProcess(
            0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, data_name  # noqa
        )

    threshold.SetLowerThreshold(lower_limit - epsilon)
    threshold.SetUpperThreshold(upper_limit + epsilon)
    threshold.AllScalarsOn()
    threshold.Update()
    result = threshold.GetOutput()
    ids = VN.vtk_to_numpy(result.GetPointData().GetGlobalIds())
    return result, ids


def threshold_vtk_data_integers(
    vtk_ugrid: vtk.vtkUnstructuredGrid,
    ints_for_thresholding: list,
    data_name: str,
) -> Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid]:
    """Threshold the vtk object using integers."""
    # find integer groups
    tag_diff = np.diff(ints_for_thresholding, prepend=ints_for_thresholding[0])
    slice_locs = np.where(tag_diff != 1)[0]

    tag_groups = []
    for ii, start_idx in enumerate(slice_locs[:]):
        if ii < len(slice_locs) - 1:
            end_idx = slice_locs[ii + 1]
            tag_groups.append(ints_for_thresholding[start_idx:end_idx])
        else:
            tag_groups.append(ints_for_thresholding[start_idx:])

    # use thresholding to extract tag groups and combine into vtk object.
    append_filter = vtk.vtkAppendFilter()
    # appendFilter.MergePointsOn() # avoids duplicate points
    for int_group in tag_groups:
        vtk_tmp, _ = threshold_vtk_data(
            vtk_ugrid,
            lower_limit=int_group[0],
            upper_limit=int_group[-1],
            data_name=data_name,
        )
        append_filter.AddInputData(vtk_tmp)

    # update with all appended data
    append_filter.Update()

    return append_filter.GetOutput()


def vtk_surface_filter(
    vtk_grid: vtk.vtkUnstructuredGrid, keep_global_ids: bool = False
) -> vtk.vtkPolyData:
    """Extract surface from a vtk object (polydata or unstructured grid)."""
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


def get_surface_info(surface: vtk.vtkPolyData) -> dict:
    """Get data from a vtk polydata surface object (filter)."""
    surface_data = {
        "points": VN.vtk_to_numpy(surface.GetPoints().GetData()),
        "connect": (VN.vtk_to_numpy(surface.GetPolys().GetConnectivityArray())).reshape(-1, 3),
        "ids_to_volume": VN.vtk_to_numpy(surface.GetPointData().GetGlobalIds()),
    }
    return surface_data


def convert_vtk_into_tetra_only(path_to_vtkfile: str):
    """Extract tetrahedrons from the source vtk file."""
    mesh = meshio.read(path_to_vtkfile)
    points = mesh.points
    elems = mesh.get_cells_type("tetra")
    meshio.write_points_cells(path_to_vtkfile, points, [("tetra", elems)])

    LOGGER.debug("Number of nodes generated: {0}".format(len(points)))
    LOGGER.debug("Number of elements generated: {0}".format(len(elems)))
    return


def get_vtk_data_field(
    vtk_grid: Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid],
    field_name: str,
    data_type: str,
) -> np.ndarray:
    """Get data field from vtk polydata or unstructured grid object.

    Parameters
    ----------
    vtk_grid : Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid]
        Vtk object from which to extract the data
    field_name : str
        Name of data field
    data_type : str
        Type of data to extract - either cell_data or point_data

    Returns
    -------
    np.array
        Numpy array of data field

    Raises
    ------
    ValueError
        Specified wrong data type
    """
    if data_type not in ["cell_data", "point_data"]:
        raise ValueError("Expecting cell_data or point_data ")

    if data_type == "cell_data":
        data = VN.vtk_to_numpy(vtk_grid.GetCellData().GetArray(field_name))
    elif data_type == "point_data":
        data = VN.vtk_to_numpy(vtk_grid.GetPointData().GetArray(field_name))

    return data


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

    Note
    ----
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

    Note
    ----
    Makes use of VoronoiKernel and mapping cell to point data
    consequently point data is mapped back to the cell.

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
        Modifies the underlying data of the target vtk object and overwrites if
        a data field with the same name is already present.
    """
    # NOTE: could use AddExcludeArray to exclude specific arrays from the interpolation

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


def rename_vtk_array(
    vtkobject: Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid],
    old_array_name: str,
    new_array_name: str,
    data_type: str = "both",
) -> Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid]:
    """Rename cell or point array of vtk object.

    Parameters
    ----------
    vtkobject : Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid]
        vtk object
    old_array_name : str
        Old name of the data array
    new_array_name : str
        New name of the data array
    data_type : str, optional
        Data types to search. Allowed options include: "cell_data", "point_data" or "both",
        by default "both"
    """
    num_cell_data = vtkobject.GetCellData().GetNumberOfArrays()
    num_point_data = vtkobject.GetPointData().GetNumberOfArrays()

    replaced = 0

    # replace cell data names
    if data_type == "cell_data" or data_type == "both":
        for ii in range(num_cell_data):
            vtk_array = vtkobject.GetCellData().GetArray(ii)
            array_name = vtkobject.GetCellData().GetArray(ii).GetName()
            if array_name == old_array_name:
                LOGGER.debug(
                    "Replacing old cell data name '{0}' with new name: '{1}'".format(
                        array_name, new_array_name
                    )
                )
                vtkobject.GetCellData().GetArray(ii).SetName(new_array_name)
                replaced += 1

    # replace point data names
    if data_type == "point_data" or data_type == "both":
        for ii in range(num_point_data):
            array_name = vtkobject.GetPointData().GetArrayName(ii)
            if array_name == old_array_name:
                LOGGER.debug(
                    "Replacing old point data name '{0}' with new name: '{1}'".format(
                        array_name, new_array_name
                    )
                )
                vtkobject.GetPointData().GetArray(ii).SetName(new_array_name)
                replaced += 1
    if replaced == 0:
        LOGGER.debug("No array names replaced")

    return vtkobject


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


def find_duplicate_elements(elements_ref: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """Find duplicate elements.

    Parameters
    ----------
    elements_ref : np.ndarray
        Array of reference elements - these are not changed
    elements : np.ndarray
        Array of elements where to check if any is already defined in the reference element array

    Returns
    -------
    duplicate_indices : np.ndarray
        Indices of the elements that are already defined in the reference array
    """
    if elements.shape[1] != elements_ref.shape[1]:
        raise ValueError("Element dimension not the same - check if 2nd dimension is consistent")

    elements_sort = np.sort(elements)
    elements_ref_sort = np.sort(elements_ref)

    # TODO: can this be improved?
    mask = np.zeros((elements_sort.shape[0], elements_ref_sort.shape[0]), dtype=bool)
    for ii, element in enumerate(elements_sort):
        mask[ii, :] = np.all(element == elements_ref_sort, axis=1)

    duplicate_indices = np.where(np.any(mask, axis=1))[0]

    return duplicate_indices


def compute_volume_stl(stl_name: str) -> float:
    """Compute the volume of a vtk PolyData Object."""
    # this avoids having to write the stl first, but directly
    # computes the volume from the polydata
    reader = vtk.vtkSTLReader()  # noqa
    reader.SetFileName(stl_name)
    reader.Update()
    # if not isinstance(vtk_object, vtk.vtkPolyData ):
    #     raise TypeError("Expecting poly data object ")

    mass_property = vtk.vtkMassProperties()  # noqa
    mass_property.SetInputData(reader.GetOutput())
    mass_property.Update()
    volume = mass_property.GetVolume()  # mm3

    return float(volume)


def vtk_unstructured_grid_to_numpy(vtk_object: vtk.vtkUnstructuredGrid):
    """[BROKEN] Extract cells and points from an arbitrary unstructured grid.

    Parameters
    ----------
    vtk_object : vtk.vtkUnstructuredGrid
        Vtk object of the unstructured grid

    Raises
    ------
    ValueError
        _description_

    Notes
    -----
    This is not supported yet
    """
    vtk_data = dsa.WrapDataObject(vtk_object)

    cells = vtk_data.Cells
    num_cells = len(vtk_data.CellTypes)
    cell_types = vtk_data.CellTypes

    supported_cell_types = {
        vtk.VTK_VERTEX: {"name": "vertex", "num_nodes": 1},
        vtk.VTK_TRIANGLE: {"name": "triangle", "num_nodes": 3},
        vtk.VTK_TETRA: {"name": "tetra", "num_nodes": 4},
    }

    cell_num_nodes = copy.deepcopy(np.array(cell_types))
    cell_num_nodes[cell_types == 10] = 4
    num_nodes_per_cell = 4
    # determine cell locations
    indices = np.arange(0, cells.shape[0], 1)
    indices_to_remove = np.arange(0, cells.shape[0], 5)

    np.delete(indices, indices_to_remove)

    np.unique(cell_types, return_index=True, return_counts=True)

    # find the start indices of the different (continuous) cell blocks
    num_detected_cell_types = len(np.unique(cell_types))
    if num_detected_cell_types > 1:
        start_indices_cellblock = np.where(np.diff(cell_types) != 0)[0] + 1
        start_indices_cellblock = np.append(0, start_indices_cellblock)
        end_indices_block = np.append(start_indices_cellblock[1:], len(cells))

    elif num_detected_cell_types == 1:
        start_indices_cellblock = np.array([0], dtype=int)
        end_indices_cellblock = np.array([len(cell_types)], dtype=int)

    num_blocks = len(start_indices_cellblock)

    all_cells = np.ones((num_cells, 5), dtype=int) * np.nan  # array with all cells

    offset = 0
    offset1 = 0

    cell_types_all = np.ones(num_cells) * np.nan

    for block in range(0, num_blocks):
        # get cell type
        cell_type = cell_types[start_indices_cellblock[block]]
        num_nodes_per_cell = supported_cell_types[cell_type]["num_nodes"]

        LOGGER.debug("\tReading {0} block".format(supported_cell_types[cell_type]["name"]))

        start_idx = start_indices_cellblock[block] + offset
        end_idx = end_indices_cellblock[block] * (num_nodes_per_cell + 1) + offset

        cell_types_all[start_idx:end_idx] = cell_type

        num_cells_in_block = (end_idx - start_idx) / (num_nodes_per_cell + 1)

        if not num_cells_in_block.is_integer():
            raise ValueError("Number of cells is not an integer value ")
        else:
            num_cells_in_block = int(num_cells_in_block)

        cells_np = np.reshape(
            cells[start_idx:end_idx],
            (num_cells_in_block, num_nodes_per_cell + 1),
        )
        cells_np = cells_np[:, 1:]

        all_cells[offset1 : num_cells_in_block + offset1, 0:num_nodes_per_cell] = cells_np

        offset = offset + (end_idx - start_idx) + 1
        offset1 = offset1 + num_cells_in_block + 1

    return


def vtk_compute_cell_area(vtk_surface: vtk.vtkPolyData) -> np.ndarray:
    """Compute area of each surface element in a polydata object.

    Parameters
    ----------
    vtk_surface : vtk.vtkPolyData
        Vtk surface

    Returns
    -------
    np.array
        Array with cell area's
    """
    n_cells = vtk_surface.GetNumberOfCells()
    cell_area = np.zeros(n_cells)
    for icell in range(n_cells):
        cell_area[icell] = vtk_surface.GetCell(icell).ComputeArea()

    return cell_area


def compute_surface_nodal_area(vtk_surface: vtk.vtkPolyData) -> np.array:
    """Compute an average nodal area by summing surface areas of connected elements.

    Parameters
    ----------
    vtk_surface : vtk.vtkPolyData
        Vtk object describing the object

    Returns
    -------
    np.array
        Numpy array with nodal areas of length number of points

    Note
    ----
    Adds the partial areas of connected elements/cells to each node.

    """
    num_points = vtk_surface.GetNumberOfPoints()
    nodal_area = np.zeros(num_points)
    cell_area = vtk_compute_cell_area(vtk_surface)
    ii = 0

    tris = get_tri_info_from_polydata(vtk_surface)[1]

    for points, area in zip(tris, cell_area):
        nodal_area[points] += area / 3
        ii += 1
    return nodal_area


def compute_surface_nodal_area_pyvista(surface: pyvista.PolyData) -> np.array:
    """Compute an average nodal area by summing surface areas of connected elements.

    Parameters
    ----------
    vtk_surface : vtk.vtkPolyData
        Vtk object describing the object

    Returns
    -------
    np.array
        Numpy array with nodal areas of length number of points

    Note
    ----
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


def find_points_inside_polydata(vtk_surface: vtk.vtkPolyData, points: np.ndarray):
    """Return indices of points that are inside the polydata object."""
    # set points
    tolerance = 1e-4
    points = vtk.vtkPolyData()
    points.SetPoints(points)

    # mark points with filter
    enclosed_points_filter = vtk.vtkSelectEnclosedPoints()
    enclosed_points_filter.SetSurfaceData(vtk_surface)
    enclosed_points_filter.SetInputData(points)
    enclosed_points_filter.SetTolerance(tolerance)
    enclosed_points_filter.Update()
    return


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


def mark_elements_inside_surfaces(
    volume_mesh: vtk.vtkUnstructuredGrid, surfaces: List[vtk.vtkPolyData]
) -> np.ndarray:
    """Mark cells based on whether they are inside the provided list of surfaces."""
    import tqdm as tqdm

    # grab centroids of each cell
    nodes, tetra, _, _ = get_tetra_info_from_unstructgrid(volume_mesh)

    centroids = np.mean(nodes[tetra, :], axis=1)

    cell_tags = np.zeros(tetra.shape[0], dtype=int) - 1

    for ii, surface in enumerate(surfaces):
        cell_ids_inside = cell_ids_inside_enclosed_surface(volume_mesh, surface)
        cell_tags[cell_ids_inside] = ii + 1

    # NOTE: very slow!
    # cell_tags of value -1 were outside all the surfaces - find the first connected tetrahedron
    LOGGER.debug("%d cells not enclosed by any of the given surfaces" % np.sum(cell_tags == -1))
    LOGGER.debug("Assigning data of closest cells")
    for cell_id in tqdm.tqdm(np.where(cell_tags == -1)[0], ascii=True):
        # use data from closest tetrahedron
        centroid = centroids[cell_id]
        sorted_cell_ids = np.argsort(np.linalg.norm(centroid - centroids, axis=1))
        closed_cell_id = sorted_cell_ids[np.argwhere(cell_tags[sorted_cell_ids] > -1).flatten()[0]]

        cell_tags[cell_id] = cell_tags[closed_cell_id]

    return cell_tags


def convert_to_polydata(vtk_ugrid: vtk.vtkUnstructuredGrid) -> vtk.vtkPolyData:
    """Use geometry filter to convert unstructured grid to polydata object.

    Parameters
    ----------
    vtk_ugrid : vtk.vtkUnstructuredGrid
        Unstructured grid object

    Returns
    -------
    vtk.vtkPolyData
        Polydata object
    """
    geom = vtk.vtkGeometryFilter()
    geom.SetInputData(vtk_ugrid)
    geom.Update()
    return geom.GetOutput()


def append_vtk_polydata_files(files: list, path_to_merged_vtk: str, substrings: List[str] = []):
    """Append a list of polydata vtk files into a single vtk file.

    Parameters
    ----------
    files : list
        List of vtk files of PolyData type
    path_to_merged_vtk : str
        Path to output vtk
    substrings : List[str], Optional
        Tags the cells using this list of substrings. Default []
    """
    from ansys.heart.preprocessor.mesh.vtkmethods import add_vtk_array
    import vtk

    # append vtk surfaces
    reader = vtk.vtkPolyDataReader()
    append = vtk.vtkAppendPolyData()
    for file in files:
        if not os.path.isfile(file):
            print("File not found...")
            continue
        reader.SetFileName(file)
        reader.Update()
        polydata = vtk.vtkPolyData()
        polydata.ShallowCopy(reader.GetOutput())

        # add cell data based on any substrings that are found
        if substrings != []:
            cell_tag = 0
            for ii, substring in enumerate(substrings):
                if substring in Path(file).name:
                    cell_tag = ii + 1

            cell_tags = np.ones(polydata.GetNumberOfCells()) * cell_tag
            add_vtk_array(
                polydata=polydata, data=cell_tags, name="tags", data_type="cell", array_type=int
            )

        append.AddInputData(polydata)

    append.Update()

    write_vtkdata_to_vtkfile(append.GetOutput(), path_to_merged_vtk)
    return


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


if __name__ == "__main__":
    print()
