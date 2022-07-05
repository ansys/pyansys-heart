"""This example shows how to extract all necessary information
from the original source data to create a simulation-ready
left-ventricle mesh, bi-ventricle mesh and four-chamber mesh
"""

import os
import numpy as np

from pathlib import Path

from ansys.heart.preprocessor.heart_model import HeartModel
from ansys.heart.preprocessor.model_information import ModelInformation
from ansys.heart.custom_logging import logger
from ansys.heart.preprocessor.vtk_module import (
    read_vtk_unstructuredgrid_file,
    write_vtkdata_to_vtkfile,
    vtk_read_mesh_file,
    vtk_map_continuous_data,
)
from vtk.numpy_interface import dataset_adapter as dsa

BASE_WORK_DIR = os.path.join(Path(__file__).parent.absolute(), "..", "workdir")


def run_preprocessor(
    model_type: str,
    database_name: str,
    path_original_mesh: str,
    work_directory: str,
    mesh_size: float = 2.0,
    remesh: bool = True,
):
    """Runs the preprocessor with the given model information

    Parameters
    ----------
    model_type : str
        Type of model
    database_name : str
        Database name. Either Strocchi2020, Cristobal2021, Strocchi2020_Modified,
        or Cristobal2021_Modified
    path_original_mesh : str
        Path to input mesh (vtk)
    work_directory : str
        Path to work directory
    """
    if not os.path.isdir(work_directory):
        os.makedirs(work_directory)

    # create model
    model_info = ModelInformation(
        model_type=model_type,
        database_name=database_name,
        path_original_mesh=path_original_mesh,
        working_directory=work_directory,
    )
    model_info.mesh_size = mesh_size

    model_info_path = os.path.join(work_directory, "model_info.json")

    model = HeartModel(model_info)
    model.extract_simulation_mesh_from_simplified_geometry()

    model.dump_model(model_info_path, clean_working_directory=False)
    return


def append_vtk_files(files: list, path_to_merged_vtk: str):
    """Appends a list of vtk files into a single vtk file

    Parameters
    ----------
    files : list
        List of vtk files of PolyData type
    path_to_merged_vtk : str
        Path to output vtk
    """

    import vtk
    from ansys.heart.preprocessor.vtk_module import add_vtk_array

    # append vtk surfaces
    reader = vtk.vtkPolyDataReader()
    append = vtk.vtkAppendFilter()
    for file in files:
        if not os.path.isfile(file):
            print("File not found...")
            continue
        reader.SetFileName(file)
        reader.Update()
        polydata = vtk.vtkPolyData()
        polydata.ShallowCopy(reader.GetOutput())

        # add cell data
        # NOTE. not general
        if "LV" in Path(file).name:
            cell_tag = 1
        elif "RV" in Path(file).name:
            cell_tag = 2
        else:
            cell_tag = 0

        cell_tags = np.ones(polydata.GetNumberOfCells()) * cell_tag
        add_vtk_array(
            polydata=polydata, data=cell_tags, name="tags", data_type="cell", array_type=int
        )
        append.AddInputData(polydata)

    append.Update()

    write_vtkdata_to_vtkfile(append.GetOutput(), path_to_merged_vtk)


if __name__ == "__main__":
    path_to_case = os.path.join(
        Path(__file__).parents[3], "downloads", "Strocchi2020_simplified", "p05.vtk"
    )

    if not os.path.isfile(path_to_case):
        append_files = True
    else:
        append_files = False

    append_files = True
    # appends vtk files
    if append_files:
        # append files
        files = []
        files.append(
            os.path.join(
                Path(__file__).parents[3],
                "downloads",
                "Strocchi2020_simplified",
                "p05_LV_volume.vtk",
            )
        )

        files.append(
            os.path.join(
                Path(__file__).parents[3],
                "downloads",
                "Strocchi2020_simplified",
                "p05_RV_volume.vtk",
            )
        )
        append_vtk_files(files, path_to_case)

    # prepare input for preprocessor
    work_directory = os.path.join(Path(path_to_case).parent, "workdir")
    model_type = "BiVentricle"
    database_name = "Strocchi2020_simplified"
    mesh_size = 1.5

    run_preprocessor(
        model_type, database_name, path_to_case, work_directory, mesh_size, remesh=True
    )

    ## add uvc point data to simulation mesh
    path_to_source = os.path.join(
        Path(__file__).parents[3], "downloads", "Strocchi2020", "05", "05.case"
    )
    path_to_target = os.path.join(
        Path(__file__).parents[3],
        "downloads",
        "Strocchi2020_simplified",
        "workdir",
        "simulation_mesh.vtk",
    )
    source_vtk = vtk_read_mesh_file(path_to_source)
    target_vtk = vtk_read_mesh_file(path_to_target)

    source_dsa = dsa.WrapDataObject(source_vtk)
    source_point_data_names = source_dsa.PointData.keys()
    uvc_array_names = [k for k in source_point_data_names if "uvc_" in k]
    target_vtk1 = vtk_map_continuous_data(
        source_vtk, target_vtk, False, array_names_to_include=uvc_array_names
    )
    write_vtkdata_to_vtkfile(target_vtk1, path_to_target)

    logger.info("** DONE **")
