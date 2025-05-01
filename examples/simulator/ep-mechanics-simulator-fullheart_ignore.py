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

Run a full-heart EP mechanics simulation
----------------------------------------
This example shows how to consume a full-heart model and
set it up for a coupled electromechanical simulation.
"""

###############################################################################
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory, heart model, and LS-DYNA executable file.

import os
from pathlib import Path

from pint import Quantity

from ansys.health.heart.examples import get_preprocessed_fullheart
import ansys.health.heart.models as models
from ansys.health.heart.settings.material.ep_material import EPMaterial
from ansys.health.heart.settings.material.material import ISO, Mat295
from ansys.health.heart.simulator import DynaSettings, EPMechanicsSimulator

# Accept the DPF license agreement.
# https://dpf.docs.pyansys.com/version/stable/getting_started/licensing.html#ref-licensing
# by setting the environment variable ``ANSYS_DPF_ACCEPT_LA`` to ``Y``.
# for instance by: os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"

# Set the working directory and path to the model. This example assumes that there is a
# preprocessed model named ``heart_model.vtu`` in the working directory.
path_to_model, path_to_partinfo, _ = get_preprocessed_fullheart(resolution="2.0mm")
workdir = Path.home() / "pyansys-heart" / "downloads" / "Rodero2021" / "01" / "FullHeart"

###############################################################################
# Load the full-heart model
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the full-heart model.
model: models.FullHeart = models.FullHeart(working_directory=workdir)
model.load_model_from_mesh(path_to_model, path_to_partinfo)

###############################################################################
# Instantiate the simulator
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Instantiate DYNA settings.
lsdyna_path = r"your_dyna_exe"  # tested with DEV-111820
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path, dynatype="intelmpi", platform="wsl", num_cpus=4
)

# Instantiate the simulator.
simulator = EPMechanicsSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "ep-mechanics"),
)

# Load default simulation settings.
simulator.settings.load_defaults()

# Compute fiber orientation in the ventricles and atria.
simulator.compute_fibers()
simulator.compute_left_atrial_fiber()
simulator.compute_right_atrial_fiber(appendage=[39, 29, 98])

# Switch the atria to active.
simulator.model.left_atrium.fiber = True
simulator.model.left_atrium.active = True

simulator.model.right_atrium.fiber = True
simulator.model.right_atrium.active = True

## Optionally, you can create more anatomical details.

# Extract elements around atrial caps and assign as a passive material.
ring = simulator.model.create_atrial_stiff_ring(radius=5)
# Material is stiff, and the value is arbitrarily chosen.
stiff_iso = Mat295(rho=0.001, iso=ISO(itype=-1, beta=2, kappa=10, mu1=0.1, alpha1=2))
ring.meca_material = stiff_iso
# Assign default EP material as for atrial.
ring.ep_material = EPMaterial.Active()
# plot the mesh
simulator.model.plot_mesh()

# Compute UHCs (Universal Heart Coordinates).
simulator.compute_uhc()

# Extract elements around atrialvenricular valves and assign as a passive material.
simulator.model.create_stiff_ventricle_base(stiff_material=stiff_iso)

# Estimate the stress-free configuration.
simulator.compute_stress_free_configuration()

# Compute the conduction system.
simulator.compute_purkinje()
simulator.compute_conduction_system()
simulator.model.plot_purkinje()


###############################################################################
# Start the main simulation
# ~~~~~~~~~~~~~~~~~~~~~~~~~
simulator.settings.mechanics.analysis.end_time = Quantity(800, "ms")
simulator.settings.mechanics.analysis.dt_d3plot = Quantity(10, "ms")

simulator.model.save_model(os.path.join(workdir, "heart_fib_beam.vtu"))

###############################################################################
# .. note::
#    A constant pressure is prescribed to the atria.
#    No circulation system is coupled with the atria.

# Start main simulation.
simulator.simulate()

###############################################################################
# View results in LS-PrePost
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

###############################################################################
# .. only:: html
#
#     .. video:: ../../_static/images/doc_Christobal01_epmeca_fh.mp4
#       :width: 600
#       :loop:
#       :class: center
