"""Collect some pre-sets for mesh extraction"""
import pathlib
import os

from ansys.heart.preprocessor._deprecated_model_information import ModelInformation
from ansys.heart.preprocessor._deprecated_heart_model import HeartModel
from conftest import get_assets_folder, clean_directory, create_directory


def workflow_extract_mesh(model_type: str, working_directory: pathlib.Path) -> pathlib.Path:
    """Uses the preprocessor the extracts a mesh of specific model type

    Parameters
    ----------
    model_type : str
        Type of model. Options include: 'LeftVentricle', 'BiVentricle', 'FourChamber'
    working_directory : pathlib.Path
        Working directory where to write model output

    Returns
    -------
    pathlib.Path
        Path to directory where all files are written
    """
    mesh_size = 2.0

    output_directory = os.path.join(working_directory, model_type)
    create_directory(output_directory)
    clean_directory(output_directory)

    case_path = os.path.join(get_assets_folder(), "cases", "strocchi2020", "01", "01.case")

    # create model info
    model_info = ModelInformation(
        model_type=model_type,
        database_name="Strocchi2020",
        path_original_mesh=case_path,
        working_directory=output_directory,
    )

    model_info.mesh_size = mesh_size

    model = HeartModel(model_info)
    model.extract_simulation_mesh()
    model.dump_model(
        os.path.join(output_directory, "model_info.json"), clean_working_directory=False
    )

    return output_directory
