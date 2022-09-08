"""This example shows how to extract all necessary information
from the original source data to create a simulation-ready
left-ventricle mesh, bi-ventricle mesh and four-chamber mesh
"""

import os
from pathlib import Path

from ansys.heart.custom_logging import LOGGER
from ansys.heart.preprocessor._deprecated_heart_model import HeartModel
from ansys.heart.preprocessor._deprecated_model_information import ModelInformation
from ansys.heart.preprocessor.vtk_module import (
    read_vtk_unstructuredgrid_file,
    vtk_map_continuous_data,
    vtk_read_mesh_file,
    write_vtkdata_to_vtkfile,
)
import numpy as np
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


if __name__ == "__main__":
    import glob as glob

    base_folder = os.path.join(Path(__file__).parents[3], "downloads", "Strocchi2020_Demo1.2")

    cases = glob.glob(os.path.join(base_folder, "*.vtk"))

    for path_to_case in cases[8:]:

        work_directory = os.path.join(Path(path_to_case).parent, Path(path_to_case).name[0:3])
        case_num = Path(path_to_case).name[1:3]

        if not os.path.isdir(work_directory):
            os.makedirs(work_directory)

        model_type = "BiVentricle"
        database_name = "Strocchi2020_simplified"
        mesh_size = 1.5

        run_preprocessor(
            model_type, database_name, path_to_case, work_directory, mesh_size, remesh=True
        )

        ## add uvc point data to simulation mesh
        path_to_source = os.path.join(
            Path(__file__).parents[3], "downloads", "Strocchi2020", case_num, case_num + ".case"
        )
        path_to_target = os.path.join(work_directory, "simulation_mesh.vtk")
        source_vtk = vtk_read_mesh_file(path_to_source)
        target_vtk = vtk_read_mesh_file(path_to_target)

        source_dsa = dsa.WrapDataObject(source_vtk)
        source_point_data_names = source_dsa.PointData.keys()
        uvc_array_names = [k for k in source_point_data_names if "uvc_" in k]
        target_vtk1 = vtk_map_continuous_data(
            source_vtk, target_vtk, False, array_names_to_include=uvc_array_names
        )
        source_cell_data_names = ["fiber", "sheet"]
        target_vtk1 = vtk_map_continuous_data(
            source_vtk, target_vtk1, False, array_names_to_include=source_cell_data_names
        )
        write_vtkdata_to_vtkfile(target_vtk1, path_to_target)

    LOGGER.info("** DONE **")
