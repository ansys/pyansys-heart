"""

EPMECA-simulator example
--------------------
This example shows you how to consume a full heart model and
set it up for electromechanical simulations.
"""

###############################################################################
# Example setup
# -------------
# before computing the fiber orientation, purkinje network we need to load
# the required modules, load a heart model and set up the simulator.
#
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory, model, and ls-dyna executable.

# sphinx_gallery_start_ignore
# Note that we need to put the thumbnail here to avoid weird rendering in the html page.
# sphinx_gallery_thumbnail_path = '_static/images/thumbnails/fh_epmeca.png'
# sphinx_gallery_end_ignore


import os
from pathlib import Path

from ansys.heart.misc.downloader import download_case, unpack_case
import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.simulator import DynaSettings, EPMechanicsSimulator
from pint import Quantity

# sphinx_gallery_start_ignore
os.environ["USE_OLD_HEART_MODELS"] = "1"
# sphinx_gallery_end_ignore

# specify necessary paths.
# Note that we need to cast the paths to strings to facilitate serialization.
case_file = str(
    Path(Path(__file__).resolve().parents[2], "downloads", "Cristobal2021", "01", "01.vtk")
)
download_folder = str(Path(Path(__file__).resolve().parents[2], "downloads"))
workdir = str(
    Path(Path(__file__).resolve().parents[2], "downloads", "Cristobal2021", "01", "FullHeart2.0")
)
path_to_model = str(Path(workdir, "heart_model.pickle"))

###############################################################################
# Download the case
# ~~~~~~~~~~~~~~~~~
# Download and unpack the case from the public repository of full hearts if it is
# not yet available. This model will be used as input for the preprocessor.
if not os.path.isfile(case_file):
    path_to_downloaded_file = download_case(
        "Cristobal2021", 1, download_folder=download_folder, overwrite=False
    )
    unpack_case(path_to_downloaded_file)

###############################################################################
# Set required information
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Set the right database to which this case belongs, and set other relevant
# information such as the desired mesh size.
info = models.ModelInfo(
    database="Cristobal2021",
    path_to_case=case_file,
    work_directory=workdir,
    path_to_model=path_to_model,
    add_blood_pool=False,
    mesh_size=2.0,
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
model = models.FullHeart(info)

# extract the simulation mesh
model.extract_simulation_mesh()

# dump the model to disk for future use
model.dump_model(path_to_model)
# print the resulting information
model.print_info()


###############################################################################
# specify LS-DYNA path
lsdyna_path = r"."  # >= DEV-105630

if not os.path.isfile(lsdyna_path):
    raise FileExistsError(f"{lsdyna_path} not found.")

# load full heart model.
model: models.FullHeart = models.HeartModel.load_model(path_to_model)

# set base working directory
model.info.workdir = str(workdir)


###############################################################################
# Instantiate the simulator object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# instantiate the simulator and settings appropriately.

# instantaiate dyna settings of choice

dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path, dynatype="intelmpi", num_cpus=8, platform="wsl"
)

# instantiate simulator object
simulator = EPMechanicsSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "epmeca"),
)

simulator.settings.load_defaults()

# compute fibers
simulator.compute_fibers()
simulator.compute_left_atrial_fiber()
simulator.compute_right_atrial_fiber(appendage=[39, 29, 98])

# switch atrial to active
simulator.model.left_atrium.has_fiber = True
simulator.model.left_atrium.is_active = True

simulator.model.right_atrium.has_fiber = True
simulator.model.right_atrium.is_active = True

# stress-free-configuration
simulator.compute_stress_free_configuration()

# Conduction system
simulator.compute_purkinje()
simulator.compute_conduction_system()


###############################################################################
# Start main simulation
# ~~~~~~~~~~~~~~~~~~~~~
simulator.settings.mechanics.analysis.end_time = Quantity(800, "ms")
simulator.settings.mechanics.analysis.dt_d3plot = Quantity(10, "ms")

simulator.model.dump_model(os.path.join(workdir, "heart_fib_beam.pickle"))

###############################################################################
# .. warning::
#    Atrial is underlying a constant (End-of-Diastolic) pressure
#    No circulation system is couple with atrial parts.

# start main simulation
simulator.simulate()

###############################################################################
# Result in LS-PrePost

###############################################################################
# .. video:: ../../_static/images/doc_Christobal01_epmeca_fh.mp4
#   :width: 600
#   :loop:
#   :class: center
