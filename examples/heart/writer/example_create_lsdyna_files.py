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
from ansys.heart.writer.dynawriter import (
    FiberGenerationDynaWriter,
    MechanicsDynaWriter,
    ZeroPressureMechanicsDynaWriter,
)
from ansys.heart.general import clean_directory

ABS_BASE_PATH = Path(__file__).parent.absolute()


def create_ls_dyna_files(path_to_model_info: str, writer_type: str, export_directory: str = None):
    """Creates the LS-DYNA files for the model specified"""

    if writer_type not in ["Mechanics", "ZeroPressure", "FiberGeneration"]:
        logger.error("Writer type %s not valid" % writer_type)
        return

    if not export_directory:
        export_directory = os.path.join(
            Path(path_to_model_info).parent.absolute(), "lsdyna_files", writer_type.lower()
        )

    clean_directory(export_directory)

    # instantiate model information and heart model
    model_info = ModelInformation()
    heart_model = HeartModel(model_info)

    # load processed heart model from info file
    heart_model.load_model(path_model_info)

    # create dyna writer based on type of writer requested
    if writer_type == "Mechanics":
        dyna_writer = MechanicsDynaWriter(
            heart_model, system_model_name="ConstantPreloadWindkesselAfterload"
        )

    elif writer_type == "ZeroPressure":
        dyna_writer = ZeroPressureMechanicsDynaWriter(heart_model)
    elif writer_type == "FiberGeneration":
        dyna_writer = FiberGenerationDynaWriter(heart_model)

    dyna_writer.update()
    dyna_writer.export(export_directory)

    logger.info("** DONE ** ")

    return


if __name__ == "__main__":

    models_to_run = ["LeftVentricle", "BiVentricle", "FourChamber", "FourChamberOriginal"]
    models_to_run = ["BiVentricle"]
    database = "Strocchi2020"

    for model in models_to_run:
        path_model_info = os.path.join(
            ABS_BASE_PATH, "..", "workdir", database, model, "model_info.json"
        )
        # path_model_info = "D:\\development\\pyheart-lib\\pyheart-lib\\downloads\\Strocchi2020\\05\\workdir\\model_info.json"
        create_ls_dyna_files(path_model_info, writer_type="Mechanics")
        create_ls_dyna_files(path_model_info, writer_type="ZeroPressure")
        create_ls_dyna_files(path_model_info, writer_type="FiberGeneration")
