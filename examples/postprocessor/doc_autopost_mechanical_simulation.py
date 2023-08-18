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

from ansys.heart.postprocessor.SystemModelPost import SystemModelPost
from ansys.heart.postprocessor.auto_process import mech_post
from ansys.heart.postprocessor.exporter import LVContourExporter
import ansys.heart.preprocessor.models as models
import numpy as np
import pyvista as pv

###############################################################################
# Set relevant paths
# ~~~~~~~~~~~~~~~~~~

path_to_model = r"D:\pyheart-lib\test_case\test_lv\model_with_fiber.pickle"

if not os.path.isfile(path_to_model):
    raise FileExistsError(f"{path_to_model} not found")

# load heart model.
model: models.LeftVentricle = models.HeartModel.load_model(path_to_model)

# set simulation path
meca_folder = r"D:\pyheart-lib\test_case\test_lv\main-mechanics"

###############################################################################
# Create PV loop
# ~~~~~~~~~~~~~~
# Pressure-volume loop figure is an important metric for heart function
system = SystemModelPost(meca_folder)
fig = system.plot_pv_loop()
fig.savefig("pv.png")

###############################################################################
# .. image:: /_static/images/pv.png
#   :width: 300pt
#   :align: center

# You can generate a series of png by setting start and end time (in second)
for it, tt in enumerate(np.linspace(0.001, 3, 60)):
    # assume heart beat once per 1s
    fig = system.plot_pv_loop(t_start=0, t_end=tt)
    fig.savefig("pv_{0:d}.png".format(it))

###############################################################################
# An animation  can be created by

# `ffmpeg -f image2 -i pv_%d.png pv_loop.mp4`

###############################################################################
# .. video:: ../../_static/images/pvloop.mp4
#   :width: 400
#   :loop:

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
# myocardium wall strain
# ~~~~~~~~~~~~~~~~~~~~~~
# Compute left ventricle strain in longitudinal, radial, circumferential directions
# todo

###############################################################################
# Run with default process scripts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
