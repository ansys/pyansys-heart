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

Run a full-heart EP simulation
------------------------------
This example shows how to consume a full-heart model and set it up for the
main EP (electrophysiology) simulation. It loads a pre-computed heart model
and computes the fiber orientation, Purkinje network, and conduction system. It
then simulates the electrophysiology.
"""

###############################################################################
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory, heart model, and LS-DYNA executable file.

import os
from pathlib import Path

from ansys.health.heart.examples import get_preprocessed_fullheart
import ansys.health.heart.models as models
from ansys.health.heart.simulator import DynaSettings, EPSimulator

# Accept the DPF license agreement.
# https://dpf.docs.pyansys.com/version/stable/getting_started/licensing.html#ref-licensing
# by setting the environment variable ``ANSYS_DPF_ACCEPT_LA`` to ``Y``.
# for instance by: os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"

# Set the working directory and path to the model. This example assumes that there is a

workdir = Path.home() / "pyansys-heart" / "downloads" / "Rodero2021" / "01" / "FullHeart"
path_to_model, path_to_partinfo, _ = get_preprocessed_fullheart(resolution="2.0mm")

###############################################################################
# Load the full-heart model
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the full-heart model.
model: models.FullHeart = models.FullHeart(working_directory=workdir)
model.load_model_from_mesh(path_to_model, path_to_partinfo)

# Save the model.
model.mesh.save(os.path.join(model.workdir, "simulation_model.vtu"))

###############################################################################
# Instantiate the simulator
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Instantiate the simulator and define settings.

# Specify the LS-DYNA path. (The last tested working version is ``intelmpi-linux-DEV-106117``.)
lsdyna_path = r"ls-dyna_msmpi.exe"

# Instantiate DYNA settings.
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path, dynatype="intelmpi", num_cpus=4, platform="wsl"
)

# Instantiate the simulator, modifying options as necessary.
simulator = EPSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "simulation-EP"),
)

###############################################################################
# Load simulation settings
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Load the default settings.
simulator.settings.load_defaults()

###############################################################################
# Compute fiber orientation
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute fiber orientation and plot the fibers on the entire model.

# Compute ventricular fibers.
simulator.compute_fibers(method="D-RBM")

# Compute atrial fibers.
simulator.model.right_atrium.active = True
simulator.model.left_atrium.active = True
simulator.model.right_atrium.fiber = True
simulator.model.left_atrium.fiber = True
simulator.compute_left_atrial_fiber()
simulator.compute_right_atrial_fiber(appendage=[39, 29, 98])
simulator.model.plot_fibers(n_seed_points=2000)

###############################################################################
# Compute the conduction system
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute the conduction system and Purkinje network, and then visualize the results.
# The action potential propagates faster through this system compared to the rest of the model.

simulator.compute_purkinje()

# By calling this method, stimulation is at the atrioventricular node.
# If you do not call this method, the two apex regions of the ventricles are stimulated.
simulator.compute_conduction_system()

simulator.model.plot_purkinje()

###############################################################################
# Start the main simulation
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Start the main EP simulation. This uses the previously computed fiber orientation
# and Purkinje network to set up and run the LS-DYNA model.


# Compute the Eikonal solution. This only computes the activation time.
simulator.settings.electrophysiology.analysis.solvertype = "ReactionEikonal"
simulator.simulate(folder_name="main-ep-ReactionEikonal")
