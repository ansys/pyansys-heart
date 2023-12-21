"""

Full heart EP-mechanics
-----------------------
This example shows you how to consume a full heart model and
set it up for a coupled electromechanical simulation.
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

import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.simulator import DynaSettings, EPMechanicsSimulator
from pint import Quantity

###############################################################################
# Example setup
# -------------
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory and generated model

# sphinx_gallery_start_ignore
# sphinx_gallery_thumbnail_path = '/_static/images/full_heart_mesh.png'
# sphinx_gallery_end_ignore


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

##############################################################################
# Set required information
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Set the right database to which this case belongs, and set other relevant
# information such as the desired mesh size.
#
# .. note::
#    The example was tested with a constant mesh size of 2.0, the stress-free
#    step may not converge when refining the mesh size.

info = models.ModelInfo(
    database="Cristobal2021",
    path_to_case=case_file,
    work_directory=workdir,
    path_to_model=path_to_model,
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
# the relevant parts and features from the original model and remesh
# the entire surface and volume. Relevant anatomical features that are extracted include
# the endocardium, epicardium, apex, and valve annuli

# instantiate a four chamber model
model = models.FullHeart(info)

# extract the simulation mesh
model.extract_simulation_mesh()

# dump the model to disk for future use
model.dump_model(path_to_model)

# print information regarding the generated model
model.print_info()

###############################################################################
# Instantiate the simulator object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# instantiate the simulator and settings appropriately.

# instantaiate dyna settings of choice
lsdyna_path = r"mppdyna_d_sse2_linux86_64_intelmmpi_105630"
dyna_settings = DynaSettings(lsdyna_path=lsdyna_path, dynatype="intelmpi", num_cpus=8)

# instantiate simulator object
simulator = EPMechanicsSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "ep-mechanics"),
)

# load default simulation settings
simulator.settings.load_defaults()

# compute fiber orientation in the ventricles and atria
simulator.compute_fibers()
simulator.compute_left_atrial_fiber()
simulator.compute_right_atrial_fiber(appendage=[39, 29, 98])

# switch atria to active
simulator.model.left_atrium.has_fiber = True
simulator.model.left_atrium.is_active = True

simulator.model.right_atrium.has_fiber = True
simulator.model.right_atrium.is_active = True

# Estimate the stress-free-configuration
simulator.compute_stress_free_configuration()

# Compute the conduction system
simulator.compute_purkinje()
simulator.compute_conduction_system()


###############################################################################
# Start main simulation
# ~~~~~~~~~~~~~~~~~~~~~
simulator.settings.mechanics.analysis.end_time = Quantity(800, "ms")
simulator.settings.mechanics.analysis.dt_d3plot = Quantity(10, "ms")

simulator.model.dump_model(os.path.join(workdir, "heart_fib_beam.pickle"))

###############################################################################
# .. note::
#    A constant pressure is prescribed to the atria.
#    No circulation system is coupled with the atria.

# start main simulation
simulator.simulate()

###############################################################################
# Result in LS-PrePost

###############################################################################
# .. video:: ../../_static/images/doc_Christobal01_epmeca_fh.mp4
#   :width: 600
#   :loop:
#   :class: center
