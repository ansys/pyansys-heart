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
Generate atrial fibers
----------------------
This example shows how to generate atrial fibers using the Laplace-Dirichlet Rule-Based
(LDRB) method.
"""

###############################################################################
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory, model, and LS-DYNA executable file. This example uses DEV-104373-g6d20c20aee.

# sphinx_gallery_start_ignore
# Note that we need to put the thumbnail here to avoid weird rendering on the HTML page.
# sphinx_gallery_thumbnail_path = '_static/images/thumbnails/atrial_fiber.png'
# sphinx_gallery_end_ignore

import os
from pathlib import Path

import numpy as np
import pyvista as pv

import ansys.health.heart.models as models
from ansys.health.heart.simulator import BaseSimulator, DynaSettings

# Specify the path to the working directory and heart model. The following path assumes
# that a preprocessed model is already available
workdir = Path.home() / "pyansys-heart" / "downloads" / "Strocchi2020" / "01" / "FourChamber"
path_to_model = str(workdir / "heart_model.vtu")

# Specify LS-DYNA path
lsdyna_path = r"ls-dyna_smp"

# Load heart model
model: models.FourChamber = models.FourChamber(working_directory=workdir)
model.load_model_from_mesh(path_to_model, path_to_model.replace(".vtu", ".partinfo.json"))


###############################################################################
# Instantiate the simulator
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Instantiate the simulator and modify options as needed.

###############################################################################
# .. note::
#    The ``DynaSettings`` object supports several LS-DYNA versions and platforms,
#    including ``smp``, ``intempi``, ``msmpi``, ``windows``, ``linux``, and ``wsl``.
#    Choose the one that works for your setup.

# instantiate LS-DYNA settings of choice
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path, dynatype="intelmpi", num_cpus=4, platform="windows"
)

simulator = BaseSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "simulation"),
)

simulator.settings.load_defaults()

# remove fiber/sheet information if it already exists
model.mesh.cell_data["fiber"] = np.zeros((model.mesh.n_cells, 3))
model.mesh.cell_data["sheet"] = np.zeros((model.mesh.n_cells, 3))

###############################################################################
# Compute atrial fibers
# ~~~~~~~~~~~~~~~~~~~~~

# compute left atrium fiber
la = simulator.compute_left_atrial_fiber()

# to compute right atrium fiber, give the appendage apex point
appendage_apex = [-33, 82, 417]
ra = simulator.compute_right_atrial_fiber(appendage_apex)

###############################################################################
# .. note::
#    You might need to define an appropriate point for the right atrial appendage.
#    The list specifies the x, y, and z coordinates close to the appendage.

###############################################################################
# Plot bundle selection results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
la.set_active_scalars("bundle")
la.plot()

###############################################################################
# .. image:: /_static/images/la_bundle.png
#   :width: 400pt
#   :align: center

ra.set_active_scalars("bundle")
ra.plot()

###############################################################################
# .. image:: /_static/images/ra_bundle.png
#   :width: 400pt
#   :align: center

###############################################################################
# Plot fibers
# ~~~~~~~~~~~
plotter = pv.Plotter()
mesh = la.ctp()
streamlines = mesh.streamlines(vectors="e_l", source_radius=50, n_points=50000)
tubes = streamlines.tube()
plotter.add_mesh(mesh, opacity=0.5, color="white")
plotter.add_mesh(tubes, color="red")
plotter.show()

###############################################################################
# .. image:: /_static/images/la_fiber.png
#   :width: 400pt
#   :align: center

plotter = pv.Plotter()
mesh = ra.ctp()
streamlines = mesh.streamlines(vectors="e_l", source_radius=50, n_points=50000)
tubes = streamlines.tube()
plotter.add_mesh(mesh, opacity=0.5, color="white")
plotter.add_mesh(tubes, color="red")
plotter.show()

###############################################################################
# .. image:: /_static/images/ra_fiber.png
#   :width: 400pt
#   :align: center

###############################################################################

# atrial fibers are automatically assigned to the heart model after computation

plotter = pv.Plotter()
mesh = model.mesh.ctp()
streamlines = mesh.streamlines(vectors="fiber", source_radius=100, n_points=50000)
tubes = streamlines.tube()
plotter.add_mesh(mesh, opacity=0.5, color="white")
plotter.add_mesh(tubes, color="red")
plotter.show()

###############################################################################
# .. image:: /_static/images/atrial_fiber_assign.png
#   :width: 400pt
#   :align: center
