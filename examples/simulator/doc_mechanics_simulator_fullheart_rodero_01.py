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

Full-heart mechanics
--------------------
This example shows you how to consume a preprocessed full heart model and
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

import os
from pathlib import Path

import ansys.heart.core.models as models
from ansys.heart.core.simulator import DynaSettings, MechanicsSimulator

# accept dpf license agreement
# https://dpf.docs.pyansys.com/version/stable/getting_started/licensing.html#ref-licensing
os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"

# set working directory and path to model. Note that we assume here that that there is a
# preprocessed model called "heart_model.vtu" available in the working directory.
workdir = Path.home() / "pyansys-heart" / "downloads" / "Rodero2021" / "01" / "FullHeart"
path_to_model = str(workdir / "heart_model.vtu")

# load the full heart.
model: models.FullHeart = models.FullHeart(working_directory=workdir)
model.load_model_from_mesh(path_to_model, path_to_model.replace(".vtu", ".partinfo.json"))

###############################################################################
# Instantiate the simulator object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# instantiate your DynaSettings and Simulator objects.
# Change options where necessary. Note that you may need to configure your environment
# variables if you choose to use a `mpi` version of LS-DYNA.

# instantiate dyna settings object
lsdyna_path = "lsdyna_intelmpi"
# instantiate dyna settings object
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path,
    dynatype="intelmpi",
    num_cpus=8,
)

# instantiate simulator object
simulator = MechanicsSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "simulation-mechanics"),
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

# # Plot the resulting fiber orientation
simulator.model.plot_fibers(n_seed_points=2000)

###############################################################################
# .. image:: /_static/images/full_heart_rodero_01_fibers.png
#   :width: 400pt
#   :align: center

###############################################################################
# Compute the stress free configuration
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute the stress free configuration. That is, when imaged under diastole
# we need to approximate the initial stress at `t=0`. The stress free configuration
# is computed through Rausch' method.

simulator.compute_stress_free_configuration()

###############################################################################
# Start main simulation
# ~~~~~~~~~~~~~~~~~~~~~
# Start the main mechanical simulation. This uses the previously computed fiber orientation
# and stress free configuration and runs the final LS-DYNA heart model.

simulator.simulate()


# sphinx_gallery_start_ignore
# Generate static images for docs.
#
from pathlib import Path

import pyvista as pv

# model.plot_fibers()
docs_images_folder = Path(Path(__file__).resolve().parents[2], "doc", "source", "_static", "images")

# Full mesh
filename = Path(docs_images_folder, "full_heart_rodero_01_fibers.png")

# fibers
plotter = pv.Plotter(off_screen=True)
plotter = pv.Plotter()
mesh = model.mesh
mesh = mesh.ctp()
streamlines = mesh.streamlines(vectors="fiber", source_radius=75, n_points=1000)
tubes = streamlines.tube()
plotter.add_mesh(mesh, opacity=0.5, color="white")
plotter.add_mesh(tubes, color="white")
plotter.camera.roll = -60
plotter.screenshot(filename)

# sphinx_gallery_end_ignore
