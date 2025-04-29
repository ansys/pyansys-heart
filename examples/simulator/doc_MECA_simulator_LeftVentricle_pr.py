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

Run a left ventricle mechanical simulation
------------------------------------------
This example shows how to generate a left-ventricle model and set it up for the
mechanical simulation. It loads the surfaces, creates the mesh and computes
the fiber orientation, assign the customized material and then simulates
the heart contraction.
"""

###############################################################################
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory, heart model, and LS-DYNA executable file.

import copy
import glob
import json
import os
from pathlib import Path
import time

import pyvista as pv

from ansys.health.heart.examples import get_input_leftventricle
import ansys.health.heart.models as models
from ansys.health.heart.post.auto_process import zerop_post
from ansys.health.heart.settings.material.material import ACTIVE, ANISO, ISO, Mat295
from ansys.health.heart.simulator import DynaSettings, MechanicsSimulator

# Accept the DPF license agreement.
# https://dpf.docs.pyansys.com/version/stable/getting_started/licensing.html#ref-licensing
# by setting the environment variable ``ANSYS_DPF_ACCEPT_LA`` to ``Y``.
# for instance by: os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"

path_to_surface, path_to_part_definitions = get_input_leftventricle()
surface = pv.read(path_to_surface)
surface.plot()

with open(path_to_part_definitions, "r") as f:
    part_definitions = json.load(f)
print(part_definitions)

###############################################################################
# Load the input surface
# ~~~~~~~~~~~~~~~~~~~~~~~
# Set the working directory and path to the model
workdir = Path.home() / "pyansys-heart" / "downloads" / "Rodero2021" / "01" / "LeftVentricle"
model: models.LeftVentricle = models.LeftVentricle(working_directory=workdir)

# Load the input surface.
model.load_input(surface, part_definitions, scalar="surface-id")

# mesh the volume
model.mesh_volume(use_wrapper=True, global_mesh_size=4.0, _global_wrap_size=4.0)
model._update_parts()
model.plot_mesh()

###############################################################################
# Instantiate the simulator
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Instantiate the simulator and define settings.

# Specify the LS-DYNA path. (The last tested working version is ``intelmpi-linux-DEV-106117``.)
# TODO: replace path
lsdyna_path = r"D:\wsl\lsdyna_mpp\ls-dyna_mpp_d_R16"

# Instantiate DYNA settings.
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path, dynatype="intelmpi", num_cpus=2, platform="wsl"
)

# Instantiate the simulator, modifying options as necessary.
simulator = MechanicsSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "simulation-MECA"),
)

###############################################################################
# Load simulation settings
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Load the default settings.
simulator.settings.load_defaults()

###############################################################################
# Compute fiber orientation
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute fiber orientation and plot the fibers on the entire model.

# Compute ventricular fibers.
simulator.compute_fibers(method="D-RBM")
simulator.model.plot_fibers(n_seed_points=2000)


myocardium = Mat295(
    rho=0.001,
    iso=ISO(itype=-3, beta=2, kappa=1.0, k1=0.20e-3, k2=6.55),
    aniso=ANISO(
        atype=-1,
        fibers=[ANISO.HGOFiber(k1=0.00305, k2=29.05), ANISO.HGOFiber(k1=1.25e-3, k2=36.65)],
        k1fs=0.15e-3,
        k2fs=6.28,
    ),
    active=ACTIVE(),  # Default active model
)
simulator.model.left_ventricle.meca_material = myocardium

simulator.compute_stress_free_configuration()

# Plot stress-free (zeropressure) geometry
report, stress_free_coord, guess_ed_coord = zerop_post(
    os.path.join(workdir, "simulation-MECA", "zeropressure"), model
)
zerop = copy.deepcopy(model.mesh)
zerop.points = stress_free_coord

#
plotter = pv.Plotter()
plotter.add_mesh(zerop, color="red", opacity=0.3, label="zeropressure")
plotter.add_mesh(simulator.model.mesh, color="grey", opacity=0.2, label="end-of-diastole")
plotter.add_legend()
plotter.show()
###############################################################################
# Start the main simulation
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Start the main meca simulation.

# run 1 heartbeat
simulator.simulate()


mesh_folder = os.path.join(workdir, "simulation-MECA", "main-mechanics", "post", "vtks")
mesh_files = sorted(glob.glob(f"{mesh_folder}/*.vtu"))
mesh_list = []
for mesh_file in mesh_files:
    mesh = pv.read(mesh_file)
    mesh.set_active_vectors("displacement")
    mesh_list.append(mesh)
# Create a PyVista plotter
plotter = pv.Plotter(notebook=True)

if mesh_files:
    actor = plotter.add_mesh(mesh_list[0])

plotter.show(interactive_update=True)


# Function to animate meshes
def animate():
    for m in mesh_list:
        plotter.clear()
        plotter.add_mesh(m)
        plotter.update()
        plotter.render()
        time.sleep(0.5)  # Control animation speed


animate()
