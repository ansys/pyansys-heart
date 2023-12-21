"""

Full-heart EP-simulator example
-------------------------------
This example shows you how to consume a full-heart model and
set it up for the main electropysiology simulation. This examples demonstrates how
you can load a pre-computed heart model, compute the fiber direction, compute the
purkinje network and conduction system and finally simulate the electrophysiology.
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

import os
from pathlib import Path

from ansys.heart.preprocessor.mesh.objects import Point
import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.simulator import DynaSettings, EPSimulator

# set working directory and path to model.
workdir = Path(Path(__file__).resolve().parents[2], "downloads", "Cristobal2021", "01", "FullHeart")

path_to_model = os.path.join(workdir, "heart_model.pickle")

if not os.path.isfile(path_to_model):
    raise FileExistsError(f"{path_to_model} not found")

# load four chamber heart model.
model: models.FullHeart = models.HeartModel.load_model(path_to_model)

# save model.
model.mesh.save(os.path.join(model.info.workdir, "simulation_model.vtu"))

# set base working directory
model.info.workdir = str(workdir)

###############################################################################
# Instantiate the simulator object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# instantiate the simulator and settings appropriately.

# instantaiate dyna settings of choice
lsdyna_path = r"mppdyna_d_sse2_linux86_64_intelmmpi_105630"
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path,
    dynatype="intelmpi",
    num_cpus=4,
)

# instantiate simulator. Change options where necessary.
simulator = EPSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "simulation-EP"),
)

###############################################################################
# Load simulation settings
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Here we load the default settings.

# Define electrode positions and add them to model
electrodes = [
    Point(name="V1", xyz=[76.53798632905277, 167.67667039945263, 384.3139099410445]),
    Point(name="V2", xyz=[64.97540262482013, 134.94983038904573, 330.4783062379255]),
    Point(name="V3", xyz=[81.20629301587647, 107.06245851801455, 320.58645260857344]),
    Point(name="V4", xyz=[85.04956217691463, 59.54502732121309, 299.2838953724169]),
    Point(name="V5", xyz=[42.31377680589025, 27.997010728192166, 275.7620409440143]),
    Point(name="V6", xyz=[-10.105919604515957, -7.176987485426985, 270.46379012676135]),
    Point(name="RA", xyz=[-29.55095501940962, 317.12543912177983, 468.91891094294414]),
    Point(name="LA", xyz=[-100.27895839242505, 135.64520460914244, 222.56688206809142]),
    Point(name="RL", xyz=[203.38825799615842, 56.19020893502452, 538.5052677637375]),
    Point(name="LL", xyz=[157.56391664248335, -81.66615972595032, 354.17867264210076]),
]
model.electrodes = electrodes

simulator.settings.load_defaults()

###############################################################################
# Compute the fiber orientation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute fiber orientation and plot the computed fibers on the entire model.

###############################################################################
# .. warning::
#    Atrial fiber orientation is approximated by apex-base direction in this model

simulator.compute_fibers()
simulator.model.plot_fibers(n_seed_points=2000)

###############################################################################
# .. image:: /_static/images/fibers.png
#   :width: 400pt
#   :align: center

###############################################################################
# Compute conduction system
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute conduction system and purkinje network and visualize.
# The action potential will propagate faster through this system
# compared to the rest of the model.

simulator.compute_purkinje()

# by calling this method, stimulation will be at the atrioventricular node
# if you skip it, the two apex regions of the ventricles will be stimulated
simulator.compute_conduction_system()

simulator.model.plot_purkinje()

###############################################################################
# .. image:: /_static/images/purkinje.png
#   :width: 400pt
#   :align: center

###############################################################################
# Start main simulation
# ~~~~~~~~~~~~~~~~~~~~~
# Start the main EP simulation. This uses the previously computed fiber orientation
# and purkinje network to set up and run the LS-DYNA model.

simulator.simulate()
