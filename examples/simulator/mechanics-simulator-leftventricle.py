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
This example shows you how to perform the following actions:

- Generate a left-ventricle model from a labeled surface.
- Set up the model for mechanical simulation:

  - Generate fibers.
  - Assign material.
  - Tune boundary conditions.
  - Compute stress-free configuration.
  - Run the simulation for one heartbeat.

- Postprocess the results:

  - Plot the stress-free configuration versus the end-diastole configuration.
  - Plot the Klotz curve.
  - Animate the simulation results.
  - Plot the PV loop.
"""

###############################################################################
# .. warning::
#    When using a standalone version of the DPF Server you will need to accept the license terms
#    https://dpf.docs.pyansys.com/version/stable/getting_started/licensing.html. The license terms
#    in the license agreement can be accepted by setting the environment variable:
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

###############################################################################
# Get the model input surface and part definitions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path_to_surface, path_to_part_definitions = get_input_leftventricle()

# Plot the input surface.
surface = pv.read(path_to_surface)
surface.plot()

# Load and print the parts defined in the JSON file.
with open(path_to_part_definitions, "r") as f:
    part_definitions = json.load(f)

# Print the part definitions using indentation for better readability.
print(json.dumps(part_definitions, indent=4))

###############################################################################
# Load the input surface and mesh the model.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set the working directory.
workdir = Path.home() / "pyansys-heart" / "downloads" / "Rodero2021" / "01" / "LeftVentricle"

# Initialize a left ventricle model.
model: models.LeftVentricle = models.LeftVentricle(working_directory=workdir)

# Load the input surface.
model.load_input(surface, part_definitions, scalar="surface-id")

# Mesh the volume.
model.mesh_volume(use_wrapper=True, global_mesh_size=4.0, _global_wrap_size=4.0)

# Update the model.
model._update_parts()

# Plot the mesh.
model.plot_mesh()

###############################################################################
# .. note::
#    This mesh is very coarse for and for demonstration purposes. In practice, the mesh size should
#    be smaller.

###############################################################################
# Instantiate the simulator
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Instantiate the simulator and define settings.

# Specify the LS-DYNA path.
lsdyna_path = r"ls-dyna_mpp"

# Instantiate LS-DYNA settings.
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path, dynatype="intelmpi", num_cpus=2, platform="windows"
)

# Instantiate the simulator, modifying options as necessary.
simulator = MechanicsSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "simulation-mechanics"),
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

# Plot the fiber orientation by streamlines.
simulator.model.plot_fibers(n_seed_points=2000)


###############################################################################
# Assign a material to the left ventricle.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define a material for the myocardium.
myocardium = Mat295(
    rho=0.001,
    iso=ISO(itype=-3, beta=2, kappa=1.0, k1=0.20e-3, k2=6.55),
    aniso=ANISO(
        atype=-1,
        fibers=[ANISO.HGOFiber(k1=0.00305, k2=29.05), ANISO.HGOFiber(k1=1.25e-3, k2=36.65)],
        k1fs=0.15e-3,
        k2fs=6.28,
    ),
    active=ACTIVE(),  # This uses the default active model
)

# Assign the material to the left ventricle.
simulator.model.left_ventricle.meca_material = myocardium

################################################################################
# Compute the stress-free configuration
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compute the stress-free configuration and estimate initial stress.
simulator.compute_stress_free_configuration()

# Plot stress-free geometry.
report, stress_free_coord, guess_ed_coord = zerop_post(
    os.path.join(workdir, "simulation-mechanics", "zeropressure"), model
)
zerop = copy.deepcopy(model.mesh)

# Update the points of the mesh with the stress-free coordinates.
zerop.points = stress_free_coord

# Plot the original and stress-free meshes.
plotter = pv.Plotter()
plotter.add_mesh(zerop, color="red", opacity=0.3, label="stress-free shape")
plotter.add_mesh(simulator.model.mesh, color="grey", opacity=0.2, label="end-of-diastole")
plotter.add_legend()
plotter.show()

# Print the stress free report in readable format.
print(json.dumps(report, indent=4))

###############################################################################
# .. note::
#    The input geometry is assumed to be at the end-of-diastole with a
#    ventricular pressure of 15 mmHg. You can modify this pressure in the settings.

###############################################################################
# Run the main simulation
# ~~~~~~~~~~~~~~~~~~~~~~~

# Tune the boundary conditions at the valve regions: in this model there are
# springs in radial and normal direction of the valve region which constrain the model.
simulator.settings.mechanics.boundary_conditions.valve["stiffness"] = Quantity(0.02, "MPa/mm")

# Start the mechanical simulation.
simulator.simulate()

###############################################################################
# .. note::
#    By default, the simulation is coupled to an open-loop circulation model with
#    a constant preload and a Windkessel-type afterload. Moreover, the simulation
#    is set to run for a single heartbeat by default.

###############################################################################
# Postprocessing
# ~~~~~~~~~~~~~~

# Plot pressure-volume loop
# ~~~~~~~~~~~~~~~~~~~~~~~~~

# Read the LS-DYNA binout file that contains cavity volumes and pressure.
icvout = ICVoutReader(os.path.join(workdir, "simulation-mechanics", "main-mechanics", "binout0000"))
pressure = icvout.get_pressure(1)
volume = icvout.get_volume(1)

# Convert to mL and mmHg and plot the pressure and volume curves.
plt.plot(volume / 1000, pressure * 7500, label="Left ventricle")
plt.xlabel("Volume (mL)")
plt.ylabel("Pressure (mmHg)")
plt.title("Pressure-volume loop")
plt.show()

# Plot the displacement field
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
mesh_folder = os.path.join(workdir, "simulation-mechanics", "main-mechanics", "post", "vtks")
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


# Animate the meshes over time
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Use a function for animation.
def animate():
    for m in mesh_list:
        plotter.clear()
        plotter.add_mesh(m, show_edges=True, opacity=0.5, scalars=None)
        plotter.update()
        plotter.render()
        time.sleep(0.1)  # Control animation speed


# Animate the model.
animate()
