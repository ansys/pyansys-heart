"""

Create a full-heart model
---------------------------------
This example shows you how to download a case from the Rodero et al (2020) database
and process that case into a simulation-ready four chamber heart model.
"""

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

import os
from pathlib import Path
import json

from ansys.heart.misc.downloader import download_case, unpack_case
import ansys.heart.preprocessor.models_new as models
from ansys.heart.simulator.support import (
    get_input_geom_and_part_defintions_from_public_database,
    preprocess_model,
)

# sphinx_gallery_start_ignore
os.environ["USE_OLD_HEART_MODELS"] = "0"
# sphinx_gallery_end_ignore

# specify necessary paths.
# Note that we need to cast the paths to strings to facilitate serialization.
case_file = str(
    Path(Path(__file__).resolve().parents[2], "downloads", "Rodero2021", "01", "01.case")
)
download_folder = str(Path(Path(__file__).resolve().parents[2], "downloads"))
workdir = str(
    Path(Path(__file__).resolve().parents[2], "downloads", "Rodero2021", "01", "FullHeartNew")
)
path_to_model = str(Path(workdir, "heart_model.pickle"))

###############################################################################
# Download the case
# ~~~~~~~~~~~~~~~~~
# Download and unpack the case from the public repository of full hearts if it is
# not yet available. This model will be used as input for the preprocessor.
if not os.path.isfile(case_file):
    path_to_downloaded_file = download_case(
        "Rodero2021", 1, download_folder=download_folder, overwrite=False
    )
    file_path = unpack_case(path_to_downloaded_file)


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
    input_geom, part_definitions = get_input_geom_and_part_defintions_from_public_database(
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

# print some info about the processed model.
model.print_info()

# clean working directory
model.info.clean_workdir([".stl", ".vtk", ".jou", ".log"])

# save cavities
for cavity in model.cavities:
    cavity.surface.save(os.path.join(model.info.workdir, "cavity_" + cavity.name + ".stl"))