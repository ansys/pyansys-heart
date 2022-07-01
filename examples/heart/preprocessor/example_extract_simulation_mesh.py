"""This example shows how to extract all necessary information
from the original source data to create a simulation-ready
left-ventricle mesh, bi-ventricle mesh and four-chamber mesh
"""


import os
from pathlib import Path

from ansys.heart.preprocessor.heart_model import HeartModel
from ansys.heart.preprocessor.model_information import ModelInformation
from ansys.heart.custom_logging import logger
from ansys.heart.general import clean_directory

# some useful global variables:
ASSET_PATH = os.path.join(
    Path(__file__).parent.absolute(), "..", "..", "..", "tests", "heart", "assets"
)
CASE_PATH_STROCCHI = os.path.join(
    os.path.abspath(ASSET_PATH), "cases", "strocchi2020", "01", "01.case"
)
CASE_PATH_CRISTOBAL = os.path.join(
    os.path.abspath(ASSET_PATH), "cases", "cristobal2021", "01", "01.vtk"
)

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
    model.extract_simulation_mesh(remesh=remesh)

    model.dump_model(model_info_path, clean_working_directory=False)
    return


if __name__ == "__main__":

    models_to_run = ["LeftVentricle", "BiVentricle", "FourChamber", "FourChamberOriginal"]
    # databases_to_run = ["Strocchi2020", "Cristobal2021", "Strocchi2020_Simplified"]
    # models_to_run = ["LeftVentricle"]
    models_to_run = ["FourChamberOriginal"]
    databases_to_run = ["Strocchi2020"]

    for database in databases_to_run:
        if database == "Strocchi2020":
            case_path = CASE_PATH_STROCCHI
        elif database == "Cristobal2021":
            case_path = CASE_PATH_CRISTOBAL

        for model_type in models_to_run:
            logger.info("***************************")
            work_directory = os.path.join(BASE_WORK_DIR, database, model_type)

            if model_type in ["LeftVentricle", "BiVentricle", "FourChamber"]:
                do_remesh = True
            elif model_type in ["FourChamberOriginal"]:
                do_remesh = False
                model_type = "FourChamber"

            run_preprocessor(
                model_type=model_type,
                database_name=database,
                path_original_mesh=case_path,
                work_directory=work_directory,
                remesh=do_remesh,
            )

            logger.info("***************************")

    logger.info("** DONE **")
