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
This example shows how to perform the following actions:
    
- Generate a left-ventricle model from a labeled surface.
- Set up the model for mechanical simulation:

  - Generate fibers.
  - Assign material.
  - Tune boundary conditions.
  - Compute stress-free configuration.
  -  Run the simulation for one heartbeat.

- Postprocess the results:

  - Plot the stress-free configuration versus the end-diastole configuration.
  - Plot the Klotz curve.
  - Animate the simulation results.
  - Plot the PV loop.
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

import matplotlib.pyplot as plt
from pint import Quantity
import pyvista as pv

from ansys.health.heart.examples import get_input_leftventricle
import ansys.health.heart.models as models
from ansys.health.heart.post.auto_process import zerop_post
from ansys.health.heart.post.dpf_utils import ICVoutReader
from ansys.health.heart.settings.material.material import ACTIVE, ANISO, ISO, Mat295
from ansys.health.heart.simulator import DynaSettings, MechanicsSimulator

# Accept the DPF license agreement.
# https://dpf.docs.pyansys.com/version/stable/getting_started/licensing.html#ref-licensing
# by setting the environment variable ``ANSYS_DPF_ACCEPT_LA`` to ``Y``.
# for instance by: os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"

###############################################################################
# Get the model input surface and part definitions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path_to_surface, path_to_part_definitions = get_input_leftventricle()

# Left ventricle is enclosed by 4 surfaces
surface = pv.read(path_to_surface)
surface.plot()

# The part is defined from the JSON file.
with open(path_to_part_definitions, "r") as f:
    part_definitions = json.load(f)
print(part_definitions)

###############################################################################
# Load the input surface
# ~~~~~~~~~~~~~~~~~~~~~~~
# Set the working directory
workdir = Path.home() / "pyansys-heart" / "downloads" / "Rodero2021" / "01" / "LeftVentricle"
model: models.LeftVentricle = models.LeftVentricle(working_directory=workdir)

# Load the input surface.
model.load_input(surface, part_definitions, scalar="surface-id")

# Mesh the volume.
model.mesh_volume(use_wrapper=True, global_mesh_size=4.0, _global_wrap_size=4.0)
###############################################################################
# .. note::
#    Mesh is over coarsed for demonstration purpose. In practice, the mesh size should be smaller.

# Update the model.
model._update_parts()

# Plot mesh.
model.plot_mesh()

###############################################################################
# Instantiate the simulator
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Instantiate the simulator and define settings.

# Specify the LS-DYNA path. (The last tested working version is ``intelmpi-linux-DEV-106117``.)
# TODO: remove before merge
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
# Compute fiber orientation and plot the fibers.
simulator.compute_fibers(method="D-RBM")
simulator.model.plot_fibers(n_seed_points=2000)


###############################################################################
# Assign the material to the left ventricle.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myocardium = Mat295(
    rho=0.001,
    iso=ISO(itype=-3, beta=2, kappa=1.0, k1=0.20e-3, k2=6.55),
    aniso=ANISO(
        atype=-1,
        fibers=[ANISO.HGOFiber(k1=0.00305, k2=29.05), ANISO.HGOFiber(k1=1.25e-3, k2=36.65)],
        k1fs=0.15e-3,
        k2fs=6.28,
    ),
    active=ACTIVE(),  # Use default active model
)
simulator.model.left_ventricle.meca_material = myocardium

################################################################################
# Compute the stress-free configuration
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###############################################################################
# .. note::
#    Input geometry is assumed to be at the end-of-diastole with a
#    pressure is set to 15 mmHg, you can modify it in the settings.

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
# the main simulation
# ~~~~~~~~~~~~~~~~~~~
###############################################################################
# .. note::
#    By defaylt, the simulation is coupled to the circulation system model which has
#    constant-preload and winkessel-type afterload
#    The simulation is set to run for 1 heartbeat by default.

# Tune the boundary conditions (springs at valves regions) which is supposed for this model.
simulator.settings.mechanics.boundary_conditions.valve["stiffness"] = Quantity(0.02, "MPa/mm")

# Start the mechanical simulation.
simulator.simulate()


###############################################################################
# Post processing
# ~~~~~~~~~~~~~~~
# TODO: remove before merge
workdir = r"C:\Users\wye\pyansys-heart\downloads\Rodero2021\01\LeftVentricle"

# Plot Pressure volume loop
# ~~~~~~~~~~~~~~~~~~~~~~~~~

# read LS-DYNA binout file
icvout = ICVoutReader(os.path.join(workdir, "simulation-MECA", "main-mechanics", "binout0000"))
pressure = icvout.get_pressure(1)
volume = icvout.get_volume(1)
plt.plot(volume / 1000, pressure * 7500, label="Left ventricle")
plt.xlabel("Volume (mL)")
plt.ylabel("Pressure (mmHg)")
plt.title("Pressure-volume loop")
plt.show()

# Plot the displacement field
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
mesh_folder = os.path.join(workdir, "simulation-MECA", "main-mechanics", "post", "vtks")
mesh_files = sorted(glob.glob(f"{mesh_folder}/*.vtu"))
mesh_list = []
for mesh_file in mesh_files:
    mesh = pv.read(mesh_file)
    mesh.set_active_vectors("displacement")
    mesh_list.append(mesh)

# Create a PyVista plotter
plotter = pv.Plotter()

if mesh_files:
    actor = plotter.add_mesh(mesh_list[0])

plotter.show(interactive_update=True)


# Function to animate meshes
def animate():
    for m in mesh_list:
        plotter.clear()
        plotter.add_mesh(m, show_edges=True, opacity=0.5, scalars=None)
        plotter.update()
        plotter.render()
        time.sleep(0.1)  # Control animation speed


animate()
