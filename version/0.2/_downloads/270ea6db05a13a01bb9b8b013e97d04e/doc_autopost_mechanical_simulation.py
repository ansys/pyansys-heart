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

Post process mechanical simulation folder
-----------------------------------------
This example shows you how to use post process script after mechanical simulation.
"""

###############################################################################
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules

# sphinx_gallery_start_ignore
# Note that we need to put the thumbnail here to avoid weird rendering in the html page.
# sphinx_gallery_thumbnail_path = '_static/images/thumbnails/pv.png'
# sphinx_gallery_end_ignore
import os
import pathlib

from ansys.heart.postprocessor.SystemModelPost import SystemModelPost
from ansys.heart.postprocessor.aha17_strain import AhaStrainCalculator
from ansys.heart.postprocessor.auto_process import mech_post
from ansys.heart.postprocessor.exporter import LVContourExporter
import ansys.heart.preprocessor.models as models
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

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
# Pressure-volume loop figure is an important metric for heart function
system = SystemModelPost(meca_folder)
fig = system.plot_pv_loop()
plt.show()

###############################################################################
# .. image:: /_static/images/pv.png
#   :width: 300pt
#   :align: center

# You can generate a series of png by setting start and end time (in second)
for it, tt in enumerate(np.linspace(0.001, 3, 60)):
    # assume heart beat once per 1s
    fig = system.plot_pv_loop(t_start=0, t_end=tt)
    fig.savefig("pv_{0:d}.png".format(it))
    plt.close()
###############################################################################
# An animation  can be created by

# `ffmpeg -f image2 -i pv_%d.png pv_loop.mp4`

###############################################################################
# .. video:: ../../_static/images/pvloop.mp4
#   :width: 400
#   :loop:
#   :class: center

###############################################################################
# Export left ventricle contour
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

exporter = LVContourExporter(os.path.join(meca_folder, "d3plot"), model)
# In case principle axis is not yet computed
model.compute_left_ventricle_anatomy_axis()

# cut from long axis 4 cavity view
cut_long = exporter.export_contour_to_vtk("l4cv", model.l4cv_axis)
# cut from short axis
cut_short = exporter.export_contour_to_vtk("short", model.short_axis)

# plot the first frame using pyvista
plotter = pv.Plotter()
plotter.add_mesh(exporter.lv_surfaces[0], opacity=0.6)
plotter.add_mesh(cut_long[0], line_width=3, color="red")
plotter.add_mesh(cut_short[0], line_width=3, color="green")
plotter.show()

###############################################################################
# .. image:: /_static/images/cut.png
#   :width: 400pt
#   :align: center

###############################################################################
# Myocardium wall strain
# ~~~~~~~~~~~~~~~~~~~~~~
# Compute left ventricle strain in longitudinal, radial, circumferential directions

# in case they are not pre-computed
model.compute_left_ventricle_anatomy_axis()
model.compute_left_ventricle_aha17()

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

# get strain for all simulation frames (this will take a while)
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
# All above steps are encapsulated in one script:

mech_post(meca_folder, model)

###############################################################################
# You can open Paraview and load the state file
# :download:`post_main2.pvsm <../../_static/others/post_main2.pvsm>`,
# and specify the folder.

###############################################################################
# .. video:: ../../_static/images/main_meca.mp4
#   :width: 600
#   :loop:
#   :class: center
