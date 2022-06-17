"""This examples creates LS-DYNA files from existing
model_info.json files and its dependencies.
Run "example_extract_simulation_mesh.py" first to create
the required files.
"""

import os
from pathlib import Path

from ansys.heart.preprocessor.heart_model import HeartModel
from ansys.heart.preprocessor.model_information import ModelInformation
from ansys.heart.custom_logging import logger
from ansys.heart.writer.dynawriter import MechanicsDynaWriter

from examples.heart.preprocessor.example_extract_simulation_mesh import clean_directory

ABS_BASE_PATH = Path(__file__).parent.absolute()


def create_ls_dyna_files(path_to_model_info: str, export_directory: str = None):
    """Creates the LS-DYNA files for the model specified
    """
    if not export_directory:
        export_directory = os.path.join(Path(path_model_info).parent.absolute(), "lsdyna_files")

    clean_directory(export_directory)

    # instantiate model ifnormation and heart model
    model_info = ModelInformation()
    heart_model = HeartModel(model_info)

    # load processed heart model from info file
    heart_model.load_model(path_model_info)

    # create dyna writer object and pass loaded model
    dyna_writer = MechanicsDynaWriter(heart_model)
    dyna_writer.update()
    dyna_writer.export(export_directory)

    logger.info("** DONE ** ")

    return


if __name__ == "__main__":

    run_all = False
    if run_all:
        models_to_run = [
            "LeftVentricle",
            "BiVentricle",
            "FourChamber"
            ]
    else:
        models_to_run = []

    if "LeftVentricle" in models_to_run:
        # left ventricle model:
        path_model_info = os.path.join(
            ABS_BASE_PATH, "..", "workdir", "left_ventricle_model", "model_info.json"
        )

        create_ls_dyna_files(path_model_info)

    if "BiVentricle" in models_to_run:
        # bi-ventricle model
        path_model_info = os.path.join(
            ABS_BASE_PATH, "..", "workdir", "bi_ventricle_model", "model_info.json"
        )
        create_ls_dyna_files(path_model_info)

    if "FourChamber" in models_to_run:
        # four chamber model
        path_model_info = os.path.join(
            ABS_BASE_PATH, "..", "workdir", "four_chamber_model", "model_info.json"
        )
        create_ls_dyna_files(path_model_info)

    # four chamber model with original mesh
    path_model_info = os.path.join(
        ABS_BASE_PATH, "..", "workdir", "four_chamber_model_original", "model_info.json"
    )
    create_ls_dyna_files(path_model_info)