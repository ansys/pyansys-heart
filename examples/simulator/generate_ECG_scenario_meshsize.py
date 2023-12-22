import os
from pathlib import Path

from ansys.heart.misc.downloader import download_case, unpack_case
import ansys.heart.preprocessor.models as models

os.environ["USE_OLD_HEART_MODELS"] = "1"

case_file = str(
    Path(Path(__file__).resolve().parents[2], "downloads", "Strocchi2020", "01", "01.case")
)
download_folder = str(Path(Path(__file__).resolve().parents[2], "downloads"))

workdir = str(
    Path(
        Path(__file__).resolve().parents[2],
        "downloads",
        "Strocchi2020",
        "01",
        "biventricle_scenario_2",
    )
)

path_to_model = os.path.join(workdir, "heart_model.pickle")
# path_to_simulation_mesh=os.path.join(workdir, "simulation_mesh.vtk")

# if not os.path.isfile(path_to_model):
#     raise FileExistsError(f"{path_to_model} not found")

if not os.path.isfile(case_file):
    path_to_downloaded_file = download_case(
        "Strocchi2020", 1, download_folder=download_folder, overwrite=False
    )
    unpack_case(path_to_downloaded_file)

# Initialize the model
info = models.ModelInfo(
    database="Strocchi2020",
    path_to_case=case_file,
    work_directory=workdir,
    path_to_model=path_to_model,
    # path_to_simulation_mesh=path_to_simulation_mesh,
    add_blood_pool=False,
    mesh_size=2,
)

# create the working directory
info.create_workdir()

# clean the working directory
info.clean_workdir(extensions_to_remove=[".stl", ".vtk", ".msh.h5"])
# dump information to stdout

info.dump_info()

model = models.BiVentricle(info)

# extract the simulation mesh
model.extract_simulation_mesh()

# dump the model to disk for future use
model.dump_model(path_to_model)
# print the resulting information
model.print_info()
