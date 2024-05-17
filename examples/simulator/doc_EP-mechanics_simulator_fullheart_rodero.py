# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
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

# set this environment variable to ensure you are using v0.2 of the model
os.environ["ANSYS_HEART_MODEL_VERSION"] = "v0.2"

import ansys.heart.preprocessor.models.v0_2.models as models
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

# accept dpf license aggrement
# https://dpf.docs.pyansys.com/version/stable/getting_started/licensing.html#ref-licensing
os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"

# specify necessary paths.
# Note that we need to cast the paths to strings to facilitate serialization.
case_file = os.path.join("pyansys-heart", "downloads", "Rodero2021", "01", "01.vtk")
workdir = os.path.join(os.path.dirname(case_file), "FullHeart")
path_to_model = os.path.join(workdir, "heart_model.pickle")

###############################################################################
# Load the full heart model
# ~~~~~~~~~~~~~~~~~~~~~~~~~

# instantiate a four chamber model
model: models.FullHeart = models.HeartModel.load_model(path_to_model)

###############################################################################
# Instantiate the simulator object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# instantiate the simulator and settings appropriately.

# instantaiate dyna settings of choice
lsdyna_path = r"mppdyna_d_sse2_linux86_64_intelmmpi_105630"
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path, dynatype="intelmpi", platform="wsl", num_cpus=6
)

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
