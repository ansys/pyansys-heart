"""This example shows how to extract all necessary information
from the original source data to create a simulation-ready
left-ventricle mesh, bi-ventricle mesh and four-chamber mesh
"""

import os
from pathlib import Path

from ansys.heart.custom_logging import LOGGER
from ansys.heart.workflow.support import run_preprocessor_improved

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


# DOWNLOAD_PATH = os.path.join(Path(__file__).parents[3], "downloads")
# CASE_PATH_STROCCHI = os.path.join(DOWNLOAD_PATH, "Strocchi2020", "05", "05.case")

if __name__ == "__main__":

    models_to_run = ["LeftVentricle", "BiVentricle", "FourChamber", "FourChamberOriginal"]
    # databases_to_run = ["Strocchi2020", "Cristobal2021", "Strocchi2020_Simplified"]
    models_to_run = ["BiVentricle"]
    models_to_run = ["FullHeart"]
    # models_to_run = ["FourChamberOriginal"]
    databases_to_run = ["Strocchi2020"]

    for database in databases_to_run:
        if database == "Strocchi2020":
            case_path = CASE_PATH_STROCCHI
        elif database == "Cristobal2021":
            case_path = CASE_PATH_CRISTOBAL

        for model_type in models_to_run:
            LOGGER.info("***************************")
            # work_directory = os.path.join(Path(case_path).parent, "workdir")
            work_directory = os.path.join(BASE_WORK_DIR, model_type)

            if model_type in ["LeftVentricle", "BiVentricle", "FourChamber", "FullHeart"]:
                do_remesh = True
            elif model_type in ["FourChamberOriginal"]:
                do_remesh = False
                model_type = "FourChamber"

            run_preprocessor_improved(
                model_type=model_type,
                database_name=database,
                path_original_mesh=case_path,
                work_directory=work_directory,
                remesh=do_remesh,
            )

            LOGGER.info("***************************")

    LOGGER.info("** DONE **")
