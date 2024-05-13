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

UHC example
--------------------
This example shows how to compute universal heart coordinates (UHC) for
the ventricles.
"""

###############################################################################
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory, model, and ls-dyna executable (uses DEV-104373-g6d20c20aee).

# sphinx_gallery_start_ignore
# Note that we need to put the thumbnail here to avoid weird rendering in the html page.
# sphinx_gallery_thumbnail_path = '_static/images/thumbnails/uvc.png'
# sphinx_gallery_end_ignore

import copy
import os

import ansys.heart.preprocessor.models.v0_2.models as models
from ansys.heart.simulator.simulator import BaseSimulator, DynaSettings
import pyvista as pv

# set this environment variable to ensure you are using v0.2 of the model
os.environ["ANSYS_HEART_MODEL_VERSION"] = "v0.2"

# specify the path to the working directory and heart model
workdir = os.path.join("pyansys-heart", "downloads", "Strocchi2020", "01", "FourChamber")
path_to_model = os.path.join(workdir, "heart_model.pickle")

# specify LS-DYNA path
lsdyna_path = r"ls-dyna_smp"

# load heart model.
model: models.FourChamber = models.HeartModel.load_model(path_to_model)

# set base working directory
model.info.workdir = str(workdir)

###############################################################################
# Instantiate the simulator object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# instantiate simulator. Change options where necessary.

# instantiate dyna settings of choice
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path,
    dynatype="smp",
    num_cpus=1,
)

simulator = BaseSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "simulation"),
)

###############################################################################
# Compute UHC
# ~~~~~~~~~~~
# Compute UHC using Laplace Dirichlet method.

simulator.compute_uhc()

###############################################################################
# .. note::
#    There are several definitions for UHC (see https://github.com/KIT-IBT/Cobiveco).
#    Here, a simple approach is taken and the
#    Dirichlet conditions are shown below. At rotational direction, the start (pi), end (-pi)
#    and middle (0) points are defined from four-cavity long axis cut view.

###############################################################################
# .. image:: /_static/images/uvc_bc.png
#   :width: 600pt
#   :align: center

###############################################################################
# Visualization of UVCs
# ~~~~~~~~~~~~~~~~~~~~~

data_ventricles = pv.read(os.path.join(workdir, "simulation", "uvc", "uvc.vtk"))

plotter = pv.Plotter(shape=(1, 3))

plotter.subplot(0, 0)
plotter.add_mesh(data_ventricles, scalars="apico-basal")

plotter.subplot(0, 1)
plotter.add_mesh(copy.copy(data_ventricles), scalars="transmural")

plotter.subplot(0, 2)
plotter.add_mesh(copy.copy(data_ventricles), scalars="rotational")
plotter.show()

###############################################################################
# .. image:: /_static/images/uvc_result.png
#   :width: 600pt
#   :align: center

###############################################################################
# Assign data to full model
# ~~~~~~~~~~~~~~~~~~~~~~~~~

# UVC is assigned back to the full model automatically
# Atrial points are padded with NaN's
model.mesh.set_active_scalars("apico-basal")
model.mesh.plot()

###############################################################################
# .. image:: /_static/images/uvc_assign.png
#   :width: 600pt
#   :align: center
