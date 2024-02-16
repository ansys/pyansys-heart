import os
from pathlib import Path

from ansys.heart.misc.downloader import download_case, unpack_case
import ansys.heart.preprocessor.models.v0_1.models as models

# sphinx_gallery_start_ignore
os.environ["USE_OLD_HEART_MODELS"] = "1"
# sphinx_gallery_end_ignore

# specify necessary paths.
# Note that we need to cast the paths to strings to facilitate serialization.
case_file = str(
    Path(Path(__file__).resolve().parents[2], "downloads", "Strocchi2020", "01", "01.case")
)
download_folder = str(Path(Path(__file__).resolve().parents[2], "downloads"))
workdir = str(
    Path(Path(__file__).resolve().parents[2], "downloads", "Strocchi2020", "01", "Biv_typical_ecgs_15")
)
path_to_model = str(Path(workdir, "heart_model.pickle"))

###############################################################################
# Download the case
# ~~~~~~~~~~~~~~~~~
# Download and unpack the case from the public repository of full hearts if it is
# not yet available. This model will be used as input for the preprocessor.
if not os.path.isfile(case_file):
    path_to_downloaded_file = download_case(
        "Strocchi2020", 1, download_folder=download_folder, overwrite=False
    )
    unpack_case(path_to_downloaded_file)

###############################################################################
# Set required information
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Set the right database to which this case belongs, and set other relevant
# information such as the desired mesh size.
info = models.ModelInfo(
    database="Strocchi2020",
    path_to_case=case_file,
    work_directory=workdir,
    path_to_model=path_to_model,
    add_blood_pool=False,
    mesh_size=1.5,
)

# create the working directory
info.create_workdir()
# clean the working directory
info.clean_workdir(extensions_to_remove=[".stl", ".vtk", ".msh.h5"])
# dump information to stdout
info.dump_info()


###############################################################################
# Initialize the heart model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize the desired heart model, and invoke the main method to
# extract the simulation mesh and dump the model to disk. This will extract
# the relevant parts from the original model and remesh the entire surface and
# volume. Moreover, relevant anatomical features are extracted.

# instantiate a four chamber model
model = models.FourChamber(info)

# extract the simulation mesh
model.extract_simulation_mesh()

# dump the model to disk for future use
model.dump_model(path_to_model)
# print the resulting information
model.print_info()