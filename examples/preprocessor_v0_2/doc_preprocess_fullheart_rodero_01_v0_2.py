"""

Create a full-heart model
---------------------------------
This example shows you how to download a case from the Rodero et al (2020) database
and process that case into a simulation-ready four chamber heart model.
"""
import os

os.environ["ANSYS_HEART_MODEL_VERSION"] = "v0.2"


###############################################################################
# Example setup
# -------------
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory and generated model

# sphinx_gallery_start_ignore
# sphinx_gallery_thumbnail_path = '_static/images/four_chamber_mesh.png'
# sphinx_gallery_end_ignore

import json
from pathlib import Path

from ansys.heart.preprocessor.models.v0_2.database_preprocessor import get_compatible_input
import ansys.heart.preprocessor.models.v0_2.models as models

# specify necessary paths.
# Note that we need to cast the paths to strings to facilitate serialization.
# case_file = str(
#     Path(Path(__file__).resolve().parents[2], "downloads", "Rodero2021", "01", "01.vtk")
# )
# download_folder = str(Path(Path(__file__).resolve().parents[2], "downloads"))
file_path = r"pyansys-heart\downloads\Rodero2021\01\01.vtk"

workdir = str(
    Path(Path(__file__).resolve().parents[2], "downloads", "Rodero2021", "01", "FullHeart_v0_2")
)

if not os.path.isdir(workdir):
    os.makedirs(workdir)

path_to_model = str(Path(workdir, "heart_model.pickle"))

path_to_input = os.path.join(workdir, "input_model.vtp")
path_to_part_definitions = os.path.join(workdir, "part_definitions.json")

if os.path.isfile(path_to_input):
    # load the existing input model.
    import pyvista as pv

    input_geom = pv.PolyData(path_to_input)
    with open(path_to_part_definitions, "r") as f:
        part_definitions = json.load(f)

else:
    # otherwise get the input geometry and part definitions
    input_geom, part_definitions = get_compatible_input(
        file_path, model_type="FullHeart", database="Rodero2021"
    )
    # save for future use.
    input_geom.save(path_to_input)
    with open(path_to_part_definitions, "w") as f:
        json.dump(part_definitions, f, indent=True)


info = models.ModelInfo(
    input=input_geom,
    scalar="surface-id",
    part_definitions=part_definitions,
    work_directory=workdir,
    mesh_size=1.5,
)

# model = preprocess_model(
#     info=info, model_type="FullHeart", clean_workdir=False, use_wrapper=True
# )

path_to_model = os.path.join(info.workdir, "heart_model.pickle")
info.path_to_model = path_to_model

# initialize full heart model
model = models.FullHeart(info)

# clean working directory
model.info.clean_workdir([".stl", ".msh.h5", ".pickle"])

# load input model
model.load_input()

# save original input as polydata
model._input.as_single_polydata.save(os.path.join(info.workdir, "input_polydata.vtp"))

# mesh the volume
model.mesh_volume(use_wrapper=True)
# update all the parts
model._update_parts()
# dump the model to disk
model.dump_model()

model.mesh.save(os.path.join(model.info.workdir, "simulation-mesh.vtu"))

# print some info about the processed model.
model.print_info()

# clean working directory
model.info.clean_workdir([".stl", ".vtk", ".jou", ".log"])

# save cavities
for cavity in model.cavities:
    cavity.surface.save(os.path.join(model.info.workdir, "cavity_" + cavity.name + ".stl"))
