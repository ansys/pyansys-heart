import os
from pathlib import Path

from ansys.heart.preprocessor._deprecated_heart_model import HeartModel
from ansys.heart.preprocessor._deprecated_model_information import ModelInformation
from ansys.heart.custom_logging import LOGGER
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
        LOGGER.error("Writer type %s not valid" % writer_type)
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
        dyna_writer = MechanicsDynaWriter(
            heart_model, system_model_name="ConstantPreloadWindkesselAfterload"
        )
    elif writer_type == "ZeroPressure":
        dyna_writer = ZeroPressureMechanicsDynaWriter(heart_model)
    elif writer_type == "FiberGeneration":
        dyna_writer = FiberGenerationDynaWriter(heart_model)

    dyna_writer.update()
    dyna_writer.export(export_directory)

    LOGGER.info("** DONE ** ")

    return


if __name__ == "__main__":

    base_folder = os.path.join(Path(__file__).parents[3], "downloads", "Strocchi2020_Demo1.2")
    case_folders = [
        # "p04",
        "p05",
        # "p06",
        # "p08",
        # "p09",
        # "p10",
        # "p12",
        # "p13",
        # "p14",
        # "p16",
        # "p19",
        # "p21",
        # "p22",
        # "p23",
        # "p24",
    ]
    for case_folder in case_folders:

        path_model_info = os.path.join(base_folder, case_folder, "model_info.json")

        create_ls_dyna_files(path_model_info, writer_type="Mechanics")
        create_ls_dyna_files(path_model_info, writer_type="ZeroPressure")
        create_ls_dyna_files(path_model_info, writer_type="FiberGeneration")
