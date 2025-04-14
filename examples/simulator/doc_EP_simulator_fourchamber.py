# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""

Run a four-chamber heart EP simulation
--------------------------------------
This example shows how to consume a four-chamber heart model and
set it up for the main EP (electropysiology) simulation. It loads a pre-computed
heart model and computes the fiber direction, Purkinje network, and conduction system.
It then simulates the electrophysiology.
"""

###############################################################################
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory, heart model, and LS-DYNA executable file.

# sphinx_gallery_start_ignore
# Note that we need to put the thumbnail here to avoid weird rendering on the HTML page.
# sphinx_gallery_thumbnail_path = '_static/images/purkinje.png'
# sphinx_gallery_end_ignore

import os
from pathlib import Path

import ansys.health.heart.models as models
from ansys.health.heart.objects import Point
from ansys.health.heart.settings.settings import DynaSettings
from ansys.health.heart.simulator import EPSimulator

# Accept the DPF license agreement.
# https://dpf.docs.pyansys.com/version/stable/getting_started/licensing.html#ref-licensing
os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"

# Set the working directory and path to the model. This example assumes that there is a
# preprocessed model named ``heart_model.vtu`` in the working directory.
workdir = Path.home() / "pyansys-heart" / "downloads" / "Strocchi2020" / "01" / "FourChamber"
path_to_model = str(workdir / "heart_model.vtu")

# Specify the LS-DYNA path. (The last tested working version is ``intelmpi-linux-DEV-106117``.)
lsdyna_path = r"ls-dyna_msmpi.exe"

# Load the four-chamber heart model.
model: models.FourChamber = models.FourChamber(working_directory=workdir)
model.load_model_from_mesh(path_to_model, path_to_model.replace(".vtu", ".partinfo.json"))

# Define electrode positions and add them to model (Positions are for patient 01 only.)
# Positions were defined using a template torso geometry.
electrodes = [
    Point(name="V1", xyz=[-29.893285751342773, 27.112899780273438, 373.30865478515625]),
    Point(name="V2", xyz=[33.68170928955078, 30.09606170654297, 380.5427551269531]),
    Point(name="V3", xyz=[56.33562469482422, 29.499839782714844, 355.533935546875]),
    Point(name="V4", xyz=[100.25729370117188, 43.61333465576172, 331.07635498046875]),
    Point(name="V5", xyz=[140.29800415039062, 81.36004638671875, 349.69970703125]),
    Point(name="V6", xyz=[167.9899139404297, 135.89862060546875, 366.18634033203125]),
    Point(name="RA", xyz=[-176.06332397460938, 57.632076263427734, 509.14202880859375]),
    Point(name="LA", xyz=[133.84518432617188, 101.44053649902344, 534.9176635742188]),
    Point(name="RL", xyz=[203.38825799615842, 56.19020893502452, 538.5052677637375]),
    Point(name="LL", xyz=[128.9441375732422, 92.85327911376953, 173.07363891601562]),
]
model.electrodes = electrodes


if not isinstance(model, models.FourChamber):
    raise TypeError("Expecting a FourChamber heart model.")

# set base working directory
model.workdir = str(workdir)

###############################################################################
# Instantiate the simulator
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Instantiate the simulator and define settings.

# Instantiate DYNA settings.
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path, dynatype="smp", num_cpus=4, platform="windows"
)

# Instantiate simulator and modify options as needed.
simulator = EPSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "simulation-EP"),
)

###############################################################################
# Load simulation settings
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Load the default simulation settings.

simulator.settings.load_defaults()

###############################################################################
# Compute UVCs
# ~~~~~~~~~~~~
# Compute UNCs (Universal Ventricular Coordinates). The transmural coordinate
# is used to define the endo, mid, and epi layers.

###############################################################################

simulator.compute_uhc()

###############################################################################
# Compute the fiber orientation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute fiber orientation and plot the computed fibers on the entire model.

###############################################################################
# .. warning::
#    The atrial fiber orientation is approximated by the apex-base direction.
#    Development is ongoing.

# Compute ventricular fibers.
simulator.compute_fibers()

# Compute atrial fibers.
simulator.model.right_atrium.active = True
simulator.model.left_atrium.active = True
simulator.model.right_atrium.fiber = True
simulator.model.left_atrium.fiber = True

# Strocchi/Rodero data has a marked left atrium appendage point.
simulator.compute_left_atrial_fiber()
# Select the right atrium appendage point.
simulator.compute_right_atrial_fiber(appendage=[-33, 82, 417])

simulator.model.plot_fibers(n_seed_points=2000)

###############################################################################
# .. image:: /_static/images/fibers.png
#   :width: 400pt
#   :align: center

###############################################################################
# Compute the conduction system
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute the conduction system and Purkinje network, and then visualize the results.
# The action potential propagates faster through this system compared to the rest of the model.

simulator.compute_purkinje()

# By calling this method, stimulation occurs at the Atrioventricular node.
# If you do not call this method, stimulation occurs at the apex nodes of the two ventricles.
simulator.compute_conduction_system()

simulator.model.plot_purkinje()

###############################################################################
# .. image:: /_static/images/purkinje.png
#   :width: 400pt
#   :align: center

###############################################################################
# Start the main simulation
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Start the main EP simulation. This uses the previously computed fiber orientation
# and Purkinje network to set up and run the LS-DYNA model with different solver
# options.

simulator.simulate()
# The two following solves only work with LS-DYNA DEV-110013 or later
simulator.settings.electrophysiology.analysis.solvertype = "Eikonal"
simulator.simulate(folder_name="main-ep-Eikonal")
simulator.settings.electrophysiology.analysis.solvertype = "ReactionEikonal"
simulator.simulate(folder_name="main-ep-ReactionEikonal")


###############################################################################
# View a plot of the transmembrane potential in LS-PrePost.

###############################################################################
# .. only:: html
#
#     .. video:: ../../_static/images/ep_4cv.mp4
#       :width: 600
#       :loop:
#       :class: center
