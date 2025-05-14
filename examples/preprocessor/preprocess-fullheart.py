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

Create a full-heart model
-------------------------
This example shows how to process a case from Rodero et al. (2021) into
a simulation-ready heart model.
"""

###############################################################################
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory and generated model.

import json
import os
from pathlib import Path

import ansys.health.heart.models as models
from ansys.health.heart.pre.database_utils import get_compatible_input

# Use Fluent 2024 R1 for meshing
import ansys.health.heart.pre.mesher as mesher
from ansys.health.heart.utils.download import download_case_from_zenodo, unpack_case

mesher._fluent_version = "24.1"

# specify a download directory
download_folder = Path.home() / "pyansys-heart" / "downloads"

# Download a compatible case from the Zenodo database.
tar_file = download_case_from_zenodo("Rodero2021", 1, download_folder, overwrite=False)
# Unpack the case to get the input CASE or VTK file.
case_file = unpack_case(tar_file)

# Specify the working directory. This code uses the directory of the CASE file.
workdir = os.path.join(os.path.dirname(case_file), "FullHeart")

if not os.path.isdir(workdir):
    os.makedirs(workdir)

# Specify paths to the model, input, and part definitions.
path_to_model = os.path.join(workdir, "heart_model.vtu")
path_to_input = os.path.join(workdir, "input_model.vtp")
path_to_part_definitions = os.path.join(workdir, "part_definitions.json")

###############################################################################
# .. note::
#    You can also manually download the CASE or VTK files from the Strocchi 2020
#    and Rodero 2021 databases. For more information, see:
#
#    - `A Publicly Available Virtual Cohort of Four-chamber Heart Meshes for
#      Cardiac Electro-mechanics Simulations <https://zenodo.org/records/3890034>`_
#    - `Virtual cohort of adult healthy four-chamber heart meshes from CT images <https://zenodo.org/records/4590294>`_
#
#    Alternatively, you can simply click one of the buttons at the bottom of this page
#    to download a CASE file for the Rodero 2021 database in an IPYNB, PY, or ZIP format.

###############################################################################
# Convert the VTK file to a compatible input format
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
input_geom, part_definitions = get_compatible_input(
    case_file, model_type="FullHeart", database="Rodero2021"
)

# Note that the input model and part definitions can be saved for later use.
# Save input geometry and part definitions.
input_geom.save(path_to_input)
with open(path_to_part_definitions, "w") as f:
    json.dump(part_definitions, f, indent=True)

###############################################################################
# Create a heart model
# ~~~~~~~~~~~~~~~~~~~~
# Create a full-heart model.
model = models.FullHeart(working_directory=workdir)

# Load input model generated in an earlier step.
model.load_input(input_geom, part_definitions, "surface-id")

# Mesh the volume of all structural parts.
model.mesh_volume(use_wrapper=True, global_mesh_size=2.0, _global_wrap_size=2.0)

# Update the model and extract the required anatomical features.
model.update()

# Optionally save the simulation mesh as a VTK object for "offline" inspection.
model.mesh.save(os.path.join(model.workdir, "simulation-mesh.vtu"))
model.save_model(os.path.join(model.workdir, "heart_model.vtu"))

# Print some information about the processed model.
print(model)

# Print part names.
print(model.part_names)

###############################################################################
# Visualize results
# ~~~~~~~~~~~~~~~~~
# Visualize and inspect the components of the model by accessing
# various properties or attributes and invoking methods.
print(f"Volume of LV cavity: {model.left_ventricle.cavity.volume} mm^3")
print(f"Volume of LV cavity: {model.left_atrium.cavity.volume} mm^3")

# Plot the remeshed model.
model.plot_mesh(show_edges=False)

# Plot the endocardial surface of the left ventricle.
model.left_ventricle.endocardium.plot(show_edges=True, color="r")

# Loop over all cavities and plot them in a single window with PyVista.
import pyvista as pv

cavities = pv.PolyData()
for c in model.cavities:
    cavities += c.surface
cavities.plot(show_edges=True)
