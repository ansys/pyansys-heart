"""This module provides some useful functions to support workflows
"""
import os
from ansys.heart.preprocessor._deprecated_heart_model import HeartModel
from ansys.heart.preprocessor._deprecated_model_information import ModelInformation
from ansys.heart.custom_logging import LOGGER


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
    LOGGER.info("##############################")
    LOGGER.info("## Launching preprocessor: ###")
    LOGGER.info("##############################")
    LOGGER.info("## Model type: {:<12} ##".format(model_type))
    LOGGER.info("## Database: {:<14} ##".format(database_name))
    LOGGER.info("## Remeshing: {:<13} ##".format(str(remesh)))
    if remesh:
        LOGGER.info("##      Size: {:<13} ##".format(mesh_size))
    LOGGER.info("##############################")
    LOGGER.info("## path mesh: %s" % path_original_mesh)
    LOGGER.info("## working directory: %s" % work_directory)
    LOGGER.info("##############################")

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
    model.extract_simulation_mesh(remesh=remesh)
    model.dump_model(model_info_path, clean_working_directory=False)

    return


def run_preprocessor_improved(
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
    LOGGER.info("##############################")
    LOGGER.info("## Launching preprocessor: ###")
    LOGGER.info("##############################")
    LOGGER.info("## Model type: {:<12} ##".format(model_type))
    LOGGER.info("## Database: {:<14} ##".format(database_name))
    LOGGER.info("## Remeshing: {:<13} ##".format(str(remesh)))
    if remesh:
        LOGGER.info("##      Size: {:<13} ##".format(mesh_size))
    LOGGER.info("##############################")
    LOGGER.info("## path mesh: %s" % path_original_mesh)
    LOGGER.info("## working directory: %s" % work_directory)
    LOGGER.info("##############################")

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
    model.extract_simulation_mesh_improved(remesh=remesh)
    model.dump_model(model_info_path, clean_working_directory=False)

    return
