"""Extract a bi-ventricle mesh from the public database of 24 pathological hearts
by Strocchi et al (2020)."""

import os
import pathlib

from ansys.heart.misc.downloader import download_case, unpack_case
import ansys.heart.preprocessor.models as models

# download case from remote repository
case_num = 1
case_path: pathlib.Path = download_case(database="Strocchi2020", case_number=case_num)
unpack_case(case_path)
case_path = os.path.join(pathlib.Path(case_path).parent, f"{case_num:02d}", f"{case_num:02d}.case")

# specify working directory
work_directory = os.path.join(pathlib.Path(__file__).parents[1], "workdir", "bi_ventricle_model")

# create model
info = models.ModelInfo(
    database="Strocchi2020",
    path_to_case=case_path,
    work_directory=work_directory,
    mesh_size=2.0,
)

# create working directory
info.create_workdir()
info.clean_workdir(remove_all=True)
info.dump_info()

# extract simulation mesh and dump model
model = models.BiVentricle(info)
model.extract_simulation_mesh()
model.print_info()
model.dump_model(os.path.join(work_directory, "bi_ventricle.pickle"))
model.info.clean_workdir([".stl", ".vtk", ".msh.h5", ".jou", ".log", ".trn"])
