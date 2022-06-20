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
CASE_PATH = os.path.join(os.path.abspath(ASSET_PATH), "cases", "cristobal2021", "01", "01.vtk")

BASE_WORK_DIR = os.path.join(Path(__file__).parent.absolute(), "..", "workdir_cristobal")

REMOVE_INTERMEDIATE_FILES = False  # flag indicating whether to remove intermediate files


def extract_leftventricle_mesh():
    """Extracts a left-ventricle mesh"""
    work_dir = os.path.join(BASE_WORK_DIR, "left_ventricle_model")
    clean_directory(work_dir)
    model_info_path = os.path.join(work_dir, "model_info.json")

    # create model
    lv_ventricle_info = ModelInformation(
        model_type="LeftVentricle",
        database_name="Cristobal2021",
        path_original_mesh=CASE_PATH,
        working_directory=work_dir,
    )

    lv_ventricle_info.mesh_size = 2.0

    left_ventricle_model = HeartModel(lv_ventricle_info)
    left_ventricle_model.extract_simulation_mesh()

    left_ventricle_model.get_model_characteristics()

    left_ventricle_model.dump_model(
        model_info_path, clean_working_directory=REMOVE_INTERMEDIATE_FILES
    )  # Toggle clean to keep/remove all intermediate files )
    return


def extract_biventricle_mesh():
    """Extracts a bi-ventricle mesh"""
    work_dir = os.path.join(BASE_WORK_DIR, "bi_ventricle_model")
    clean_directory(work_dir)
    model_info_path = os.path.join(work_dir, "model_info.json")

    # create model
    bi_ventricle_info = ModelInformation(
        model_type="BiVentricle",
        database_name="Cristobal2021",
        path_original_mesh=CASE_PATH,
        working_directory=work_dir,
    )

    bi_ventricle_info.mesh_size = 2.0
    biventricle_model = HeartModel(bi_ventricle_info)
    biventricle_model.extract_simulation_mesh()
    biventricle_model.get_model_characteristics()

    biventricle_model.dump_model(
        model_info_path, clean_working_directory=REMOVE_INTERMEDIATE_FILES
    )  # Toggle clean to keep/remove all intermediate files )
    return


def extract_fourchamber_mesh():
    """Extracts a four chamber mesh"""
    work_dir = os.path.join(BASE_WORK_DIR, "four_chamber_model")
    clean_directory(work_dir)
    model_info_path = os.path.join(work_dir, "model_info.json")

    # create model
    fourchamber_info = ModelInformation(
        model_type="FourChamber",
        database_name="Cristobal2021",
        path_original_mesh=CASE_PATH,
        working_directory=work_dir,
    )
    fourchamber_info.mesh_size = 2.0

    four_chamber_model = HeartModel(fourchamber_info)
    four_chamber_model.extract_simulation_mesh()
    four_chamber_model.get_model_characteristics()

    four_chamber_model.dump_model(
        model_info_path, clean_working_directory=REMOVE_INTERMEDIATE_FILES
    )  # Toggle clean to keep/remove all intermediate files
    return


def extract_fourchamber_mesh_noremesh():
    """Extracts a four chamber mesh"""
    work_dir = os.path.join(BASE_WORK_DIR, "four_chamber_model_original")
    clean_directory(work_dir)
    model_info_path = os.path.join(work_dir, "model_info.json")

    # create model
    fourchamber_info = ModelInformation(
        model_type="FourChamber",
        database_name="Cristobal2021",
        path_original_mesh=CASE_PATH,
        working_directory=work_dir,
    )

    four_chamber_model = HeartModel(fourchamber_info)
    four_chamber_model.extract_simulation_mesh( remesh = False )
    four_chamber_model.get_model_characteristics()

    four_chamber_model.dump_model(
        model_info_path, clean_working_directory=REMOVE_INTERMEDIATE_FILES
    )
    return

if __name__ == "__main__":

    if not os.path.isdir(BASE_WORK_DIR):
        os.mkdir(BASE_WORK_DIR)

    run_all = True
    if run_all:
        models_to_run = [
            # "LeftVentricle",
            "BiVentricle",
            # "FourChamber",
            # "FourChamberOriginal"
            ]
    else:
        models_to_run = []

    # extract left ventricle mesh
    if "LeftVentricle" in models_to_run:
        logger.info("***************************")
        extract_leftventricle_mesh()

    if "BiVentricle" in models_to_run:
        # extract biventricle mesh
        logger.info("***************************")
        extract_biventricle_mesh()

    if "FourChamber" in models_to_run:
        # extract four chamber mesh
        logger.info("***************************")
        extract_fourchamber_mesh()

    if "FourChamberOriginal" in models_to_run:
        logger.info("***************************")
        extract_fourchamber_mesh_noremesh()

    logger.info("***************************")
    logger.info("** DONE **")
