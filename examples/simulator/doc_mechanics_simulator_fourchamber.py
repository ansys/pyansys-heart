"""

Mechanics-simulator example
--------------------
This example shows you how to consume a four-cavity heart model and
set it up for the main mechanical simulation. This examples demonstrates how
you can load a pre-computed heart model, compute the fiber direction, compute the
stress free configuration, and finally simulate the mechanical model.
"""

###############################################################################
# Example setup
# -------------
# before computing the fiber orientation, and stress free configuration we
# need to load the required modules, load a heart model and configure the
# mechanical simulator.
#
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory, model, and ls-dyna executable.

# sphinx_gallery_start_ignore
# Note that we need to put the thumbnail here to avoid weird rendering in the html page.
# sphinx_gallery_thumbnail_path = '_static/images/thumbnails/frame_0043.png'
# sphinx_gallery_end_ignore

import os
from pathlib import Path

import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.simulator import MechanicsSimulator

# set working directory and path to model.
workdir = Path(
    Path(__file__).resolve().parents[2], "downloads", "Strocchi2020", "01", "FourChamber"
)
path_to_model = os.path.join(workdir, "heart_model.pickle")

if not os.path.isfile(path_to_model):
    raise FileExistsError(f"{path_to_model} not found")

# specify LS-DYNA path
lsdyna_path = Path(Path(__file__).resolve().parents[4], "dyna-versions", "ls-dyna_smp_d.exe")

if not os.path.isfile(lsdyna_path):
    raise FileExistsError(f"{lsdyna_path} not found.")

# load four chamber heart model.
model: models.FourChamber = models.HeartModel.load_model(path_to_model)

if not isinstance(model, models.FourChamber):
    raise TypeError("Expecting a FourChamber heart model.")

# set base working directory
model.info.workdir = str(workdir)

###############################################################################
# Instantiate the simulator object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# instantiate simulator. Change options where necessary. Note that you may need
# to configure your environment variables if you choose `mpp`.

simulator = MechanicsSimulator(
    model=model,
    lsdynapath=lsdyna_path,
    dynatype="smp",
    num_cpus=4,
    simulation_directory=os.path.join(workdir, "simulation-mechanics"),
    platform="windows",
)

###############################################################################
# Load simulation settings
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Here we load the default settings.

simulator.settings.load_defaults()

###############################################################################
# Compute the fiber orientation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute fiber orientation and plot the computed fibers on the entire model.

simulator.compute_fibers()
simulator.model.plot_fibers(n_seed_points=2000)

###############################################################################
# .. image:: /_static/images/fibers.png
#   :width: 400pt
#   :align: center

###############################################################################
# Compute the stress free configuration
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute the stress free configuration. That is, when imaged under diastole
# we need to approximate the initial stress at `t=0`. The stress free configuration
# is computed through Rausch' method.

simulator.compute_stress_free_configuration()
# plot the updated model.
simulator.model.plot_mesh(show_edges=True)

###############################################################################
# .. image:: /_static/images/stress_free.png
#   :width: 400pt
#   :align: center

###############################################################################
# Start main simulation
# ~~~~~~~~~~~~~~~~~~~~~
# Start the main mechanical simulation. This uses the previously computed fiber orientation
# and stress free configuration and runs the final LS-DYNA heart model.

###############################################################################
# .. warning::
#    In 4-chamber model, atria is passive and only to provide realistic boundary conditions.

simulator.simulate()

###############################################################################
# We can plot deformation and active stress (history variable 22) in LS-PrePost

###############################################################################
# .. video:: ../../_static/images/meca_4cv.mp4
#   :width: 600
#   :loop:
#   :class: center
