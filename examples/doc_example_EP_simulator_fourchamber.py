"""

EP-simulator example
--------------------
This example shows you how to consume a four-cavity heart model and 
set it up for the main electropysiology simulation. This examples demonstrates how 
you can load a pre-computed heart model, compute the fiber direction, compute the 
purkinje network and conduction system and finally simulate the electrophysiology.
"""

# sphinx_gallery_thumbnail_path = 'images/purkinje.png'

###############################################################################
# Example setup
# -------------
# before computing the fiber orientation, purkinje network we need to load
# the required modules, load a heart model and set up the simulator.
#
# Peform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory, model, and ls-dyna executable.

import os
from pathlib import Path
import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.simulator import EPSimulator

# set working directory and path to model.
workdir = Path(Path(__file__).parents[3], "downloads", "Strocchi2020", "01", "FourChamber")
path_to_model = os.path.join(workdir, "heart_model.pickle")

if not os.path.isfile(path_to_model):
    raise FileExistsError(f"{path_to_model} not found")

# specify LS-DYNA path
lsdyna_path = Path(Path(__file__).parents[5], "dyna-versions", "ls-dyna_smp_d.exe")

if not os.path.isfile(lsdyna_path):
    raise FileExistsError(f"{lsdyna_path} not found.")

# load four chamber heart model.
model: models.FourChamber = models.HeartModel.load_model(path_to_model)

if not isinstance(model, models.FourChamber):
    raise TypeError("Expecting a FourChamber heart model.")

# set base working directoy
model.info.workdir = str(workdir)

###############################################################################
# Instantiate the simulator object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# instantiate simulator. Change options where necessary.

simulator = EPSimulator(
    model=model,
    lsdynapath=lsdyna_path,
    dynatype="smp",
    num_cpus=1,
    simulation_directory=os.path.join(workdir, "simulation-EP"),
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
# .. image:: /images/fibers.png
#   :width: 400pt
#   :align: center

###############################################################################
# Compute conduction system
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute conduction system and purkinje network and visualize.
# The action potential will propogate faster through this system
# compared to the rest of the model.

simulator.compute_purkinje()
simulator.compute_conduction_system()
simulator.model.plot_purkinje()

###############################################################################
# .. image:: /images/purkinje.png
#   :width: 400pt
#   :align: center

###############################################################################
# Start main simulation
# ~~~~~~~~~~~~~~~~~~~~~
# Start the main EP simulation. This uses the previously computed fiber orientation
# and purkinje network to set up and run the LS-DYNA model.

simulator.simulate()
