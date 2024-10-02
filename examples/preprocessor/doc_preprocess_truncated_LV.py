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

Create a truncated ellipsoid model
----------------------------------
This example shows you how to build a basic ellipsoidal model
from primitive shapes. Shape based on
Land et al (2015): https://doi.org/10.1098/rspa.2015.0641
"""


###############################################################################
# Example setup
# -------------
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory and generated model
import os

import numpy as np
import pyvista as pv

import ansys.heart.core.models as models

# sphinx_gallery_start_ignore
# sphinx_gallery_thumbnail_path = '_static/images/truncated_LV_mesh.png'
# sphinx_gallery_end_ignore


###############################################################################
# Create a truncated ellipsoid using pyvista
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
workdir = r"simplified-geometries\truncated_lv"

# create the surfaces of a truncated LV model
ellips_epi: pv.PolyData = pv.ParametricEllipsoid(xradius=10, yradius=10, zradius=20)
ellips_endo: pv.PolyData = pv.ParametricEllipsoid(xradius=7, yradius=7, zradius=17)

# clip ellips at z=5
z_truncate = 5  # z-coordinate to truncate at
ellips_endo = ellips_endo.clip(normal="z", origin=[0, 0, z_truncate])
ellips_epi = ellips_epi.clip(normal="z", origin=[0, 0, z_truncate])

# compute x and y radius to create a closing disc.
endo_bounds = ellips_endo.extract_feature_edges().bounds
epi_bounds = ellips_epi.extract_feature_edges().bounds

base: pv.PolyData = pv.Disc(
    center=(0, 0, z_truncate), inner=endo_bounds[1], outer=epi_bounds[1], c_res=200
).triangulate()

# add "surface-id" to cell data
base.cell_data["surface-id"] = 3
ellips_endo.cell_data["surface-id"] = 1
ellips_epi.cell_data["surface-id"] = 2

# combine into single poly data object.
heart: pv.PolyData = ellips_endo + ellips_epi + base
heart.plot(show_edges=True)

###############################################################################
# .. image:: /_static/images/truncated_LV_mesh_input.png
#   :width: 400pt
#   :align: center

###############################################################################
# Convert the input to a HeartModel
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# construct part definition dictionary
part_definitions = {
    "Left ventricle": {
        "id": 1,
        "enclosed_by_boundaries": {
            "left-ventricle-endocardium": 1,
            "left-ventricle-epicardium": 2,
            "interface_left-ventricle-myocardium_mitral-valve": 3,
        },
    }
}

# use the combined polydata `heart` as input, where "surface-id" identifies each
# of the relevant regions.
# part definitions is used to map the remeshed model to the HeartModel parts/boundaries
info = models.ModelInfo(
    input=heart,
    scalar="surface-id",
    part_definitions=part_definitions,
    work_directory=workdir,
    mesh_size=0.5,
)
info.path_to_model = os.path.join(info.workdir, "heart_model.pickle")

# initialize left-ventricular heart model
model = models.LeftVentricle(info)

# clean working directory
model.info.clean_workdir([".stl", ".msh.h5", ".pickle"])

# load input model
model.load_input()

###############################################################################
# Remesh the surfaces and volume
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# .. note::
#    The individual surfaces in the combined PolyData object are
#    unconnected. Using the wrapper automatically fixes any small gaps
#    and ensures proper connectivity.

# remesh the model using wrapping
model.mesh_volume(use_wrapper=True)

# assign axis of model manually.
model.l4cv_axis = {"center": base.center, "normal": np.array([1, 0, 0])}
model.l2cv_axis = {"center": base.center, "normal": np.array([0, 1, 0])}
model.short_axis = {"center": base.center, "normal": np.array([0, 0, 1])}

# update the model
model._sync_input_parts_to_model_parts()

model._assign_elements_to_parts()
model._assign_surfaces_to_parts()

model._validate_parts()
model._validate_surfaces()

model._add_nodal_areas()
model._add_surface_normals()

model._assign_cavities_to_parts()
model._update_cap_names()
model._validate_cap_names()
model._extract_apex()

model.compute_left_ventricle_aha17()

# dump the model to disk
model.dump_model()

# plot the clipped volume mesh.
model.mesh.clip(crinkle=True).plot(show_edges=True)

###############################################################################
# .. image:: /_static/images/truncated_LV_mesh.png
#   :width: 400pt
#   :align: center

# sphinx_gallery_start_ignore
# Generate static images for docs.
#
from pathlib import Path

docs_images_folder = Path(Path(__file__).resolve().parents[2], "doc", "source", "_static", "images")

filename = Path(docs_images_folder, "truncated_LV_mesh_input.png")
# model.plot_surfaces()
plotter = pv.Plotter(off_screen=True)
# plotter.view_zx()
# plotter.add_mesh(model.mesh, show_edges=False, color="white")
# plotter.add_mesh(model.mesh.clip(), show_edges=True, color="w")
plotter.add_mesh(heart, show_edges=True)
# plotter.view_yz()
# plotter.add_legend()
# plotter.camera.roll = -60
plotter.screenshot(filename)


# Full mesh
filename = Path(docs_images_folder, "truncated_LV_mesh.png")
# model.plot_surfaces()
plotter = pv.Plotter(off_screen=True)
# plotter.view_zx()
# plotter.add_mesh(model.mesh, show_edges=False, color="white")
plotter.add_mesh(model.mesh.clip(crinkle=True), show_edges=True, color="w")
# plotter.view_yz()
# plotter.add_legend()
# plotter.camera.roll = -60
plotter.screenshot(filename)

# sphinx_gallery_end_ignore
