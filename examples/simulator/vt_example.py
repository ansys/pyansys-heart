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

Run a full-heart electrophysiology simulation
---------------------------------------------
This example shows how to consume a full-heart model and set it up for the
main electrophysiology simulation. It loads a pre-computed heart model
and computes the fiber orientation, Purkinje network, and conduction system. It
then simulates the electrophysiology.
"""

###############################################################################
# .. warning::
#    When using a standalone version of the DPF Server, you must accept the `license terms
#    <https://dpf.docs.pyansys.com/version/stable/getting_started/licensing.html>`_. To
#    accept these terms, you can set this environment variable:
#
#    .. code-block:: python
#
#        import os
#        os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"

###############################################################################
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory, heart model, and LS-DYNA executable file.

import os
from pathlib import Path

import numpy as np
from pint import Quantity
import pyvista as pv

from ansys.health.heart.examples import get_preprocessed_fullheart
import ansys.health.heart.models as models
from ansys.health.heart.settings.material.ep_material import CellModel, EPMaterial
from ansys.health.heart.settings.settings import Stimulation
from ansys.health.heart.simulator import DynaSettings, EPSimulator

# Set the working directory and path to the model. This example assumes that there is a

workdir = Path.home() / "pyansys-heart" / "downloads" / "Rodero2021" / "01" / "FullHeart-VT"
path_to_model, path_to_partinfo, _ = get_preprocessed_fullheart(resolution="2.0mm")

###############################################################################
# Load the full-heart model
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the full-heart model.
model: models.FullHeart = models.HeartModel.load_model(
    path_to_model, path_to_partinfo, working_directory=workdir
)

# Save the model.
model.mesh.save(os.path.join(model.workdir, "simulation_model.vtu"))

###############################################################################
# Instantiate the simulator object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# instantiate the simulator and settings appropriately.

# specify LS-DYNA path (last tested working versions is intelmpi-linux-DEV-106117)
# lsdyna_path = r"ls-dyna_msmpi.exe"
lsdyna_path = r"lsdyna.exe"

# instantaiate dyna settings of choice
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path, dynatype="intelmpi", num_cpus=6, platform="windows"
)

# instantiate simulator. Change options where necessary.
simulator = EPSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "simulation-EP"),
)

###############################################################################
# Load simulation settings
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Here we load the default settings.

simulator.settings.load_defaults()

###############################################################################
# Compute the fiber orientation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute fiber orientation and plot the computed fibers on the entire model.

###############################################################################
# .. warning::
#    Atrial fiber orientation is approximated by apex-base direction in this model

# compute ventricular fibers
simulator.compute_fibers()

# compute atrial fibers
simulator.model.right_atrium.active = True
simulator.model.left_atrium.active = True
simulator.model.right_atrium.fiber = True
simulator.model.left_atrium.fiber = True
simulator.compute_left_atrial_fiber()
simulator.compute_right_atrial_fiber(appendage=[39, 29, 98])
simulator.model.plot_fibers(n_seed_points=2000)

###############################################################################
# .. image:: /_static/images/fibers.png
#   :width: 400pt
#   :align: center

###############################################################################
# Build scar and slow conduction channel geometry
# ~~~~~~~~~~~~~~~~~~~~~~~~~


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


simulator.compute_uhc()

mesh = simulator.model.mesh.copy()

# removing surface meshes
mesh.remove_cells(mesh.celltypes != 10, inplace=True)

uvc_point_id = get_point_from_uvc(model, apicobasal=0.5, transmural=0.5, rotational=np.pi)
uvc_stimpoint = mesh.points[uvc_point_id, :]
sphere = pv.Sphere(center=uvc_stimpoint, radius=20)
plane = pv.Plane(
    center=uvc_stimpoint, direction=simulator.model.l4cv_axis["normal"], i_size=60, j_size=60
)

pl = pv.Plotter()
pl.add_mesh(sphere, color="red", opacity=0.3)
pl.add_mesh(plane, color="blue", opacity=0.3)
pl.add_mesh(mesh, opacity=0.3)
pl.show()

labels = np.array([0] * mesh.n_cells)

scar_tet_ids = np.nonzero(
    mesh.compute_implicit_distance(sphere).ptc().cell_data["implicit_distance"] < 0
)

labels[scar_tet_ids] = 2

bz_tet_ids = np.nonzero(
    np.all(
        [
            np.abs(mesh.compute_implicit_distance(plane).ptc().cell_data["implicit_distance"]) < 2,
            labels == 2,
        ],
        axis=0,
    )
)

labels[bz_tet_ids] = 1

mesh["test"] = labels
mesh.plot(scalars="test")

###############################################################################
# Build scar and slow conduction channel parts
# ~~~~~~~~~~~~~~~~~~~~~~~~~

slow_ten_tusscher = CellModel.Tentusscher(
    gkr=0.0459,
    gna=5.63844,
    gcal=1.2338e-5,
    gks=0.0784,
)

epinsulator = EPMaterial.Insulator()
epbz = ep_mat_active = EPMaterial.Active(
    sigma_fiber=0.7,
    sigma_sheet=0.2,
    sigma_sheet_normal=0.2,
    beta=140,
    cm=0.01,
    cell_model=slow_ten_tusscher,
)

eplv = ep_mat_active = EPMaterial.Active(
    sigma_fiber=0.5, sigma_sheet=0.1, sigma_sheet_normal=0.1, beta=140, cm=0.01
)

simulator.model.get_part("Left ventricle").ep_material = eplv

scar = simulator.model.create_part_by_ids(scar_tet_ids, "scar")
scar.fiber = False
scar.active = False
scar.ep_material = epinsulator

boundary_zone = simulator.model.create_part_by_ids(bz_tet_ids, "boundary_zone")
boundary_zone.fiber = True
boundary_zone.active = True
boundary_zone.ep_material = epbz

###############################################################################
# Add stimulation settings
# ~~~~~~~~~~~~~~~~~~~~~~~~~

uvc_point_id = get_point_from_uvc(model, apicobasal=0, transmural=0, rotational=np.pi)
stim_tet_id = mesh.point_cell_ids(uvc_point_id)[0]
stim_point_ids = mesh.cells.reshape((-1, 5))[stim_tet_id][1:]

# TODO: add multistim

stim_uvc = Stimulation(
    node_ids=list(stim_point_ids), t_start=0, period=230, duration=5, amplitude=50
)
simulator.settings.electrophysiology.stimulation = {"stim": stim_uvc}

###############################################################################
# Change simulation settings
# ~~~~~~~~~~~~~~~~~~~~~~~~~

simulator.settings.electrophysiology.analysis.dt_d3plot = Quantity(50)
simulator.settings.electrophysiology.analysis.end_time = Quantity(1000)

###############################################################################
# Start main simulation
# ~~~~~~~~~~~~~~~~~~~~~
# Start the main EP simulation. This uses the previously computed fiber orientation
# and purkinje network to set up and run the LS-DYNA model.


# # switch to ReactionEikonal
simulator.settings.electrophysiology.analysis.solvertype = "ReactionEikonal"
simulator.simulate(folder_name="main-ep-ReactionEikonal")
