"""Extracts a bi-ventricle mesh"""

import os
import pathlib

import ansys.heart.preprocessor as preproc
from ansys.heart.preprocessor._deprecated_heart_model import HeartModel
from ansys.heart.preprocessor._deprecated_model_information import ModelInformation

# get path to case
case_path = os.path.join(
    pathlib.Path(preproc.__file__).parents[3],
    "tests",
    "heart",
    "assets",
    "cases",
    "strocchi2020",
    "01",
    "01.case",
)

# specify working directory
work_directory = os.path.join(pathlib.Path(__file__).parents[0], "workdir", "bi_ventricle_model")
model_info_path = os.path.join(work_directory, "model_info.json")

# create work directory
if not os.path.isdir(work_directory):
    os.makedirs(work_directory)

# create model
model_info = ModelInformation(
    model_type="BiVentricle",
    database_name="Strocchi2020",
    path_original_mesh=case_path,
    working_directory=work_directory,
)

model_info.mesh_size = 2.0

model = HeartModel(model_info)
model.extract_simulation_mesh()
model.get_model_characteristics()

model.dump_model(model_info_path, clean_working_directory=True)
