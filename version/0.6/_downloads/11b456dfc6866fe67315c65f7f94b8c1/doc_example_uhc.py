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

import pyvista as pv

import ansys.heart.core.models as models
from ansys.heart.simulator.simulator import BaseSimulator, DynaSettings

# specify the path to the working directory and heart model
workdir = os.path.abspath(os.path.join("downloads", "Strocchi2020", "01", "FourChamber"))

# sphinx_gallery_start_ignore
# Overwrite with env variables: for testing purposes only. May be removed by user.
from pathlib import Path

try:
    case_file = str(Path(os.environ["PATH_TO_CASE_FILE"]))
    workdir = os.path.join(os.path.dirname(case_file), "FourChamber")
    path_to_dyna = str(Path(os.environ["PATH_TO_DYNA"]))

except KeyError:
    pass
# sphinx_gallery_end_ignore

path_to_model = os.path.join(workdir, "heart_model.vtu")

# specify LS-DYNA path
lsdyna_path = r"ls-dyna_smp"

# load heart model.
model: models.FourChamber = models.FourChamber(models.ModelInfo(work_directory=workdir))
model.load_model_from_mesh(path_to_model, path_to_model.replace(".vtu", ".partinfo.json"))
model._extract_apex()
model.compute_left_ventricle_anatomy_axis()
model.compute_left_ventricle_aha17()

# set base working directory
model.info.workdir = str(workdir)

###############################################################################
# Instantiate the simulator object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# instantiate simulator. Change options where necessary.

# instantiate dyna settings of choice
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path,
    dynatype="intelmpi",
    num_cpus=1,
)

# sphinx_gallery_start_ignore
# Overwrite with env variables: for testing purposes only. May be removed by user.
try:
    dyna_settings.lsdyna_path = path_to_dyna
    # assume we are in WSL if .exe not in path.
    if ".exe" not in path_to_dyna:
        dyna_settings.platform = "wsl"
except Exception:
    pass
# sphinx_gallery_end_ignore

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
# UVC is assigned back to the full model automatically
# Atrial points are padded with NaN's

plotter = pv.Plotter(shape=(1, 3))

plotter.subplot(0, 0)
plotter.add_mesh(simulator.model.mesh, scalars="apico-basal")

plotter.subplot(0, 1)
plotter.add_mesh(copy.copy(simulator.model.mesh), scalars="transmural")

plotter.subplot(0, 2)
plotter.add_mesh(copy.copy(simulator.model.mesh), scalars="rotational")
plotter.show()

###############################################################################
# .. image:: /_static/images/uvc_result.png
#   :width: 600pt
#   :align: center
