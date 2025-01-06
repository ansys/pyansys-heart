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

Stimulation definition example
---------------------------------
This example shows you how to define an EP stimulation. It demonstrates how you
can load a pre-computed heart model, define a stimulation region based on a sphere
centered on the apex, and a sphere centered on a point chosen in Universal
Ventricular Coordinates (UVC).
"""

###############################################################################
# Example setup
# -------------
# Loading required modules and heart model.
#
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths.

# sphinx_gallery_start_ignore
# Note that we need to put the thumbnail here to avoid weird rendering in the html page.
# sphinx_gallery_thumbnail_path = '_static/images/stimulation.png'
# sphinx_gallery_end_ignore

import os

import numpy as np
import pyvista as pv

import ansys.heart.core.models as models
from ansys.heart.simulator.settings.settings import SimulationSettings, Stimulation
from ansys.heart.simulator.simulator import DynaSettings, EPSimulator

# accept dpf license aggrement
# https://dpf.docs.pyansys.com/version/stable/getting_started/licensing.html#ref-licensing
os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"

# sphinx_gallery_start_ignore
# Overwrite with env variables: for testing purposes only. May be removed by user.

# sphinx_gallery_end_ignore
workdir = os.path.join(
    "pyansys-heart", "downloads", "Strocchi2020", "01", "FourChamber", "teststim"
)
path_to_model = os.path.join(workdir, "heart_model.vtu")


# load your four chamber heart model with uvcs (see preprocessor examples to create
# a heart model from scratch)
model: models.FourChamber = models.FourChamber(working_directory=workdir)
model.load_model_from_mesh(path_to_model)

###############################################################################
# Define stimulation at the apex
# ------------------------------
# Select points inside sphere centered at the left apex.
apex_left = model.left_ventricle.apex_points[0].xyz
sphere = pv.Sphere(center=(apex_left), radius=2)
newdata = model.mesh.select_enclosed_points(sphere)
node_ids = np.where(newdata.point_data["SelectedPoints"] == 1)[0]
apex_stim_points = model.mesh.points[node_ids, :]

pl = pv.Plotter()
pl.add_points(apex_stim_points, color="red")
pl.add_mesh(model.mesh, color="lightgrey", opacity=0.2)
pl.show()

# Define stimulation and introduce it as simulation settings
stim_apex = Stimulation(node_ids=list(node_ids), t_start=0, period=800, duration=20, amplitude=50)
settings = SimulationSettings()
settings.load_defaults()
settings.electrophysiology.stimulation = {"stim_apex": stim_apex}


# Define auxiliary function to find a point in the model based on its UVC coordinates
def get_point_from_uvc(
    model: models.HeartModel, apicobasal: float, transmural: float, rotational: float
):
    point_coords = np.array([apicobasal, transmural, rotational])
    diffs = (
        np.transpose(
            np.vstack(
                (
                    model.mesh.point_data["apico-basal"],
                    model.mesh.point_data["transmural"],
                    model.mesh.point_data["rotational"],
                )
            )
        )
        - point_coords
    )

    norms = np.linalg.norm(diffs, axis=1)
    norms[np.isnan(norms)] = 1000
    point_id = np.argmin(norms)
    return point_id


###############################################################################
# Define stimulation based on UVC
# -------------------------------
# Select points inside sphere centered at a chosen point based on UVC coordinates
# (if the model has UVC).
if (
    ("transmural" in model.mesh.point_data.keys())
    and ("apico-basal" in model.mesh.point_data.keys())
    and ("rotational" in model.mesh.point_data.keys())
):
    uvc_point_id = get_point_from_uvc(model, apicobasal=0.7, transmural=0, rotational=1)
    uvc_stimpoint = model.mesh.points[uvc_point_id, :]
    sphere = pv.Sphere(center=(uvc_stimpoint), radius=2)
    newdata = model.mesh.select_enclosed_points(sphere)
    node_ids = np.where(newdata.point_data["SelectedPoints"] == 1)[0]
    uvc_stim_points = model.mesh.points[node_ids, :]

    pl = pv.Plotter()
    pl.add_points(apex_stim_points, color="red")
    pl.add_points(uvc_stim_points, color="blue")
    pl.add_mesh(model.mesh, color="lightgrey", opacity=0.2)
    pl.show()

    stim_uvc = Stimulation(
        node_ids=list(node_ids), t_start=0, period=800, duration=20, amplitude=50
    )
    settings.electrophysiology.stimulation["stim_uvc"] = stim_uvc

# specify LS-DYNA path
lsdyna_path = r"ls-dyna_msmpi.exe"


# instantaiate dyna settings of choice
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path, dynatype="intelmpi", num_cpus=4, platform="wsl"
)

simulator = EPSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "simulation-EP"),
)
simulator.settings = settings

simulator.simulate()
