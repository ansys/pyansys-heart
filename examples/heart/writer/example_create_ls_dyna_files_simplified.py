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

    # instantiate model ifnormation and heart model
    model_info = ModelInformation()
    heart_model = HeartModel(model_info)

    # load processed heart model from info file
    heart_model.load_model(path_model_info)

    # create dyna writer based on type of writer requested
    if writer_type == "Mechanics":
        dyna_writer = MechanicsDynaWriter(heart_model)
    elif writer_type == "ZeroPressure":
        dyna_writer = ZeroPressureMechanicsDynaWriter(heart_model)
    elif writer_type == "FiberGeneration":
        dyna_writer = FiberGenerationDynaWriter(heart_model)

    dyna_writer.update()
    dyna_writer.export(export_directory)

    logger.info("** DONE ** ")

    return


if __name__ == "__main__":

    path_model_info = "D:\\development\\pyheart-lib\\pyheart-lib\\downloads\\Strocchi2020_simplified\\workdir\\model_info.json"

    create_ls_dyna_files(path_model_info, writer_type="Mechanics")
    create_ls_dyna_files(path_model_info, writer_type="ZeroPressure")
    create_ls_dyna_files(path_model_info, writer_type="FiberGeneration")
