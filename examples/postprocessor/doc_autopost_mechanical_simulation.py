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

Postprocess the mechanical simulation folder
-----------------------------------------
This example shows how to use the postprocess script after a mechanical simulation.
"""

###############################################################################
# Perform required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~

# sphinx_gallery_start_ignore
# Note that we need to put the thumbnail here to avoid weird rendering in the html page.
# sphinx_gallery_thumbnail_path = '_static/images/thumbnails/pv.png'
# sphinx_gallery_end_ignore
import os
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

import ansys.health.heart.models as models
from ansys.health.heart.post.auto_process import mech_post
from ansys.health.heart.post.strain_calculator import AhaStrainCalculator
from ansys.health.heart.post.system_model_post import SystemModelPost

###############################################################################
# Set relevant paths
# ~~~~~~~~~~~~~~~~~~

path_to_model = r"D:\pyansys-heart\test_case\test_lv\model_with_fiber.pickle"

if not os.path.isfile(path_to_model):
    raise FileExistsError(f"{path_to_model} not found")

# load heart model.
model: models.LeftVentricle = models.HeartModel.load_model(path_to_model)

# set simulation path
meca_folder = pathlib.Path(r"D:\pyansys-heart\test_case\test_lv\main-mechanics")

###############################################################################
# Create PV loop
# ~~~~~~~~~~~~~~
# A PV (Pressure-Volume) loop figure is an important metric for heart function.
system = SystemModelPost(meca_folder)
fig = system.plot_pv_loop()
plt.show()

###############################################################################
# .. image:: /_static/images/pv.png
#   :width: 300pt
#   :align: center

# You can generate a series of PNG files by setting start and end times (in seconds)
for it, tt in enumerate(np.linspace(0.001, 3, 60)):
    # assume heart beat once per 1s
    fig = system.plot_pv_loop(t_start=0, t_end=tt)
    fig.savefig("pv_{0:d}.png".format(it))
    plt.close()
###############################################################################
# You can create an animation.

# `ffmpeg -f image2 -i pv_%d.png pv_loop.mp4`

###############################################################################
# .. only:: html
#
#     .. video:: ../../_static/images/pvloop.mp4
#       :width: 400
#       :loop:
#       :class: center


###############################################################################
# Compute myocardium wall strain
# ~~~~~~~~~~~~~~~~~~~~~~
# Compute left ventricle strain in longitudinal, radial, circumferential directions.

aha_evaluator = AhaStrainCalculator(model, d3plot_file=meca_folder / "d3plot")
# get LRC strain at a given time and export a file named LRC_10.vtk
strain17_at10 = aha_evaluator.compute_aha_strain_at(frame=10, out_dir=".")

# show generated vtk
aha = pv.read(r"LRC_10.vtk")
aha.set_active_scalars("AHA")
aha.plot()

###############################################################################
# .. image:: /_static/images/aha17.png
#   :width: 400pt
#   :align: center

# bulleye plot for strain
fig, ax = plt.subplots(figsize=(24, 16), nrows=1, ncols=3, subplot_kw=dict(projection="polar"))
fig.canvas.manager.set_window_title("Left Ventricle Bulls Eyes (AHA)")
for i in range(3):
    aha_evaluator.bullseye_plot(ax[i], strain17_at10[:, i])
ax[0].set_title("longitudinal")
ax[1].set_title("radial")
ax[2].set_title("circumferential")
plt.show()

###############################################################################
# .. image:: /_static/images/aha17_strain.png
#   :width: 400pt
#   :align: center

# get strain for all simulation frames, which takes a while
strain_table = aha_evaluator.compute_aha_strain(out_dir=".", write_vtk=False)

# plot
l_strain_base = np.mean(strain_table[:, 1:19:3], axis=1)
l_strain_mid = np.mean(strain_table[:, 19:37:3], axis=1)
l_strain_apical = np.mean(strain_table[:, 37::3], axis=1)

plt.plot(strain_table[:, 0], l_strain_base, label="Longitudinal strain @Basal")
plt.plot(strain_table[:, 0], l_strain_mid, label="Longitudinal strain @MidCavity")
plt.plot(strain_table[:, 0], l_strain_apical, label="Longitudinal strain @Apical")
plt.legend()
plt.show()

###############################################################################
# .. image:: /_static/images/l_strain_curve.png
#   :width: 400pt
#   :align: center

###############################################################################
# Run with default process scripts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# All of the preceding steps are encapsulated in one script:

mech_post(meca_folder, model)

###############################################################################
# You can open Paraview, load the state file
# :download:`post_main2.pvsm <../../_static/others/post_main2.pvsm>`,
# and specify the directory.

###############################################################################
# .. only:: html
#
#     .. video:: ../../_static/images/main_meca.mp4
#       :width: 600
#       :loop:
#       :class: center
