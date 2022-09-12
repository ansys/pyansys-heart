"""Extracts a bi-ventricle mesh"""

import os
import pathlib

import ansys.heart.preprocessor as preproc
import ansys.heart.preprocessor.models as models

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


# create model
info = models.ModelInfo(
    database="Strocchi2020",
    path_to_case=case_path,
    work_directory=work_directory,
    mesh_size=2.0,
)

info.create_workdir()
info.clean_workdir(remove_all=True)
info.dump_info()

model = models.BiVentricle(info)
model.extract_simulation_mesh()
model.print_info()
model.dump_model(os.path.join(work_directory, "bi_ventricle.pickle"))
model.info.clean_workdir([".stl", ".vtk", ".msh.h5", ".jou", ".log", ".trn"])
