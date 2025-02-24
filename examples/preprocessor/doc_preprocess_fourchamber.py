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

Create a four chamber heart model
---------------------------------
This example shows you how to process a case file from the Strocchi2020 database
and process that into a simulation-ready full heart model.
"""

###############################################################################
# Example setup
# -------------
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory and generated model

# sphinx_gallery_start_ignore
# sphinx_gallery_thumbnail_path = '_static/images/four_chamber_mesh.png'
# sphinx_gallery_end_ignore

import json
import os

from ansys.heart.core.helpers.general import clean_directory
import ansys.heart.core.models as models
from ansys.heart.preprocessor.database_preprocessor import get_compatible_input

# Use Fluent 24.1 for meshing.
import ansys.heart.preprocessor.mesher as mesher

mesher._fluent_version = "24.1"

# specify necessary paths.
case_file = os.path.join("downloads", "Strocchi2020", "01", "01.case")

# sphinx_gallery_start_ignore
# Overwrite with env variables: for testing purposes only. May be removed by user.
from pathlib import Path

try:
    case_file = str(Path(os.environ["PATH_TO_CASE_FILE"]))
except KeyError:
    pass
# sphinx_gallery_end_ignore

workdir = os.path.join(os.path.dirname(case_file), "FourChamber")

if not os.path.isdir(workdir):
    os.makedirs(workdir)

path_to_model = os.path.join(workdir, "heart_model.vtu")
path_to_input = os.path.join(workdir, "input_model.vtp")
path_to_part_definitions = os.path.join(workdir, "part_definitions.json")

###############################################################################
# .. note::
#    You may need to (manually) download the .case or .vtk files from the Strocchi2020
#    and Rodero2021 databases first. See:
#
#    - https://zenodo.org/records/3890034
#    - https://zenodo.org/records/4590294
#
#    Alternatively you can make use of the download
#    module instead. See the download module.

###############################################################################
# Convert the .vtk file into compatible input
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
input_geom, part_definitions = get_compatible_input(
    case_file, model_type="FourChamber", database="Strocchi2020"
)

# Note that the input model and part definitions can be used for later use.
# save input geometry and part definitions:
input_geom.save(path_to_input)
with open(path_to_part_definitions, "w") as f:
    json.dump(part_definitions, f, indent=True)

###############################################################################
# Set required information
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Set the right database to which this case belongs, and set other relevant
# information such as the desired mesh size.

# create or clean working directory
if not os.path.isdir(workdir):
    os.makedirs(workdir)

clean_directory(workdir, [".stl", ".msh.h5", ".pickle"])

###############################################################################
# Initialize the heart model with info
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize the desired heart model with info.

# initialize four chamber heart model
model = models.FourChamber(working_directory=workdir)

# load input model generated in an earlier step.
model.load_input(input_geom, part_definitions, "surface-id")

# mesh the volume of all structural parts.
model.mesh_volume(use_wrapper=True, global_mesh_size=1.5)

# update the model and extract the required (anatomical) features
model._update_parts()

# dump the model to disk
model.save_model(path_to_model)

# Optionally save the simulation mesh as a vtk object for "offline" inspection
model.mesh.save(os.path.join(model.workdir, "simulation-mesh.vtu"))
model.save_model(os.path.join(model.workdir, "heart_model.vtu"))

# print some info about the processed model.
print(model)

# clean the working directory
clean_directory(workdir, extensions_to_remove=[".stl", ".vtk", ".msh.h5"])

# print part names
print(model.part_names)

###############################################################################
# Visualize results
# ~~~~~~~~~~~~~~~~~
# You can visualize and inspect the components of the model by accessing
# various properties/attributes and invoke methods.
print(f"Volume of LV cavity: {model.left_ventricle.cavity.volume} mm^3")
print(f"Volume of LV cavity: {model.left_atrium.cavity.volume} mm^3")

# plot the remeshed model
model.plot_mesh(show_edges=False)

###############################################################################
# .. image:: /_static/images/four_chamber_mesh.png
#   :width: 400pt
#   :align: center

# plot the endocardial surface of the left ventricle.
model.left_ventricle.endocardium.plot(show_edges=True, color="r")

###############################################################################
# .. image:: /_static/images/four_chamber_lv_endocardium.png
#   :width: 400pt
#   :align: center

# loop over all cavities and plot these in a single window.
import pyvista as pv

cavities = pv.PolyData()
for c in model.cavities:
    cavities += c.surface
cavities.plot(show_edges=True)

###############################################################################
# .. image:: /_static/images/four_chamber_cavities.png
#   :width: 400pt
#   :align: center

# sphinx_gallery_start_ignore
# Generate static images for docs.
#
from pathlib import Path

docs_images_folder = Path(Path(__file__).resolve().parents[2], "doc", "source", "_static", "images")

# Full mesh
filename = Path(docs_images_folder, "four_chamber_mesh.png")
plotter = pv.Plotter(off_screen=True)
model.mesh.set_active_scalars("_volume-id")
plotter.add_mesh(model.mesh, show_edges=False)
# plotter.show()
plotter.camera.roll = -60
plotter.screenshot(filename)

# Clipped full mesh

# left ventricle endocardium
filename = Path(docs_images_folder, "four_chamber_lv_endocardium.png")
plotter = pv.Plotter(off_screen=True)
model.mesh.set_active_scalars(None)
plotter.add_mesh(model.left_ventricle.endocardium, color="r", show_edges=True)
plotter.camera.roll = -60
plotter.screenshot(filename)

# Cavities
filename = Path(docs_images_folder, "four_chamber_cavities.png")
plotter = pv.Plotter(off_screen=True)
for c in model.cavities:
    plotter.add_mesh(c.surface, show_edges=True)
# plotter.show()
plotter.camera.roll = -60
plotter.screenshot(filename)
# sphinx_gallery_end_ignore
