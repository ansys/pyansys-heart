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
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
What does Mrs Jones heart look like?
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
"""

# -----------------------------------------------------------------------------------
# 1. Download patient geometry
# -----------------------------------------------------------------------------------
import glob

from ansys.heart.core.helpers.downloader import (download_case_from_zenodo,
                                                 unpack_case)

# Download the tar file of Rodero2021 from the Zenodo database.
tar_file = download_case_from_zenodo("Rodero2021", 5, "downloads")
# Unpack the tar file
unpack_case(tar_file)
# list all files
glob.glob("downloads" + "/**/*.*", recursive=True)


# -----------------------------------------------------------------------------------
# 2. Process patient geometry
# -----------------------------------------------------------------------------------
import json
import os
from pathlib import Path

import ansys.heart.core.models as models
# Use Fluent 24.1 for meshing.
import ansys.heart.preprocessor.mesher as mesher
from ansys.heart.core.helpers.general import clean_directory
from ansys.heart.preprocessor.database_preprocessor import get_compatible_input

mesher._fluent_version = "24.1"
# specify necessary paths.
case_file = r"D:\REPOS\pyheart\pyansys-heart\downloads\Rodero2021\05\05.vtk"
workdir = os.path.join(os.path.dirname(case_file), "FullHeart")

if not os.path.isdir(workdir):
    os.makedirs(workdir)

path_to_model = os.path.join(workdir, "heart_model.pickle")
path_to_input = os.path.join(workdir, "input_model.vtp")
path_to_part_definitions = os.path.join(workdir, "part_definitions.json")

###############################################################################
# Convert the .vtk file into compatible input
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
input_geom, part_definitions = get_compatible_input(
    case_file, model_type="FullHeart", database="Rodero2021"
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
else:
    clean_directory(workdir, [".stl", ".msh.h5", ".pickle"])

###############################################################################
# Initialize the heart model with info
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize the desired heart model with info.

# initialize full heart model
model = models.FullHeart(working_directory=workdir)

# load input model generated in an earlier step.
model.load_input(input_geom, part_definitions, "surface-id")

# mesh the volume of all structural parts.
model.mesh_volume(use_wrapper=True, global_mesh_size=1.5)

# update the model and extract the required (anatomical) features
model._update_parts()

# Optionally save the simulation mesh as a vtk object for "offline" inspection
model.mesh.save(os.path.join(model.workdir, "simulation-mesh.vtu"))
model.save_model(os.path.join(model.workdir, "heart_model.vtu"))

# print some info about the processed model.
print(model)

# clean the working directory
clean_directory(workdir, extensions_to_remove=[".stl", ".vtk", ".msh.h5"])

# print part names
print(model.part_names)

# -----------------------------------------------------------------------------------
# 3. Visualize patient geometry
# -----------------------------------------------------------------------------------
print(f"Volume of LV cavity: {model.left_ventricle.cavity.volume} mm^3")
print(f"Volume of LV cavity: {model.left_atrium.cavity.volume} mm^3")
# plot the remeshed model
model.plot_mesh(show_edges=False)


"""
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
How are the fibers organized in this heart?
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
"""

import os

import ansys.heart.core.models as models
from ansys.heart.simulator.settings.material.ep_material import EPMaterial
from ansys.heart.simulator.settings.material.material import NeoHookean
from ansys.heart.simulator.simulator import (DynaSettings,
                                             EPMechanicsSimulator, EPSimulator)
from pint import Quantity

os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"

# specify necessary paths.
# Note that we need to cast the paths to strings to facilitate serialization.
case_file = os.path.join("pyansys-heart", "downloads", "Rodero2021", "05", "01.vtk")
workdir = r"D:\REPOS\pyheart\pyansys-heart\downloads\Rodero2021\05\FullHeart"

path_to_model = os.path.join(workdir, "heart_model.vtu")
# -----------------------------------------------------------------------------------
# 1. Load model
# -----------------------------------------------------------------------------------

# instantiate a four chamber model
model: models.FullHeart = models.FullHeart(working_directory=workdir)
model.load_model_from_mesh(path_to_model, path_to_model.replace(".vtu", ".partinfo.json"))


# -----------------------------------------------------------------------------------
# 2. Create a simulator
# -----------------------------------------------------------------------------------

# instantaiate dyna settings of choice
lsdyna_path = r"D:\dynaexe\mppdyna"  # tested with DEV-111820
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path, dynatype="intelmpi", platform="wsl", num_cpus=6
)


# instantiate simulator object

simulator = EPMechanicsSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "ep-mechanics-eikonal"),
)

# load default simulation settings
simulator.settings.load_defaults()

# -----------------------------------------------------------------------------------
# 3. Compute fiber orientation
# -----------------------------------------------------------------------------------
simulator.compute_fibers(method="D-RBM")
simulator.compute_left_atrial_fiber()
simulator.compute_right_atrial_fiber(appendage=[30, 32, 86])



# switch atria to active
simulator.model.left_atrium.fiber = True
simulator.model.left_atrium.active = True

simulator.model.right_atrium.fiber = True
simulator.model.right_atrium.active = True

# Optionally, we can create more anatomical details.
# Sometimes, it's in favor of convergence rate of mechanical solve

# Extract elements around atrial caps and assign as a passive material
ring = simulator.model.create_atrial_stiff_ring(radius=5)
# material is stiff and value is arbitrarily chosen
ring.meca_material = NeoHookean(rho=0.001, c10=0.1, nu=0.499)
# assign default EP material as for atrial
ring.ep_material = EPMaterial.Active()

# Compute universal coordinates:
simulator.compute_uhc()

# Extract elements around atrialvenricular valves and assign as a passive material
simulator.model.create_stiff_ventricle_base(stiff_material=NeoHookean(rho=0.001, c10=0.1, nu=0.499))

# Estimate the stress-free-configuration
simulator.compute_stress_free_configuration()
# -----------------------------------------------------------------------------------
# 4. Plot fibers
# -----------------------------------------------------------------------------------
simulator.model.plot_fibers()

"""
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
Can you highlight the conduction system?
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
"""

# Compute the conduction system and plot
simulator.compute_purkinje()
simulator.compute_conduction_system()
simulator.model.plot_purkinje()



"""
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
can you show me the electrical wave propagation?
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
"""
# -----------------------------------------------------------------------------------
# 1. Create electrophysiology simulator
# -----------------------------------------------------------------------------------
simulator_ep = EPSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "simulation-EP"),
)


# -----------------------------------------------------------------------------------
# 2. Apply settings and run simulation
# -----------------------------------------------------------------------------------
simulator_ep.model.get_part('Aorta').ep_material= EPMaterial.Insulator()
simulator_ep.model.get_part('Pulmonary artery').ep_material= EPMaterial.Insulator()
simulator_ep.settings.electrophysiology.analysis.solvertype = "ReactionEikonal"
simulator_ep.simulate(folder_name="main-ep-ReactionEikonal")


# -----------------------------------------------------------------------------------
# 3. Animation of electrical wave propagation to be performed in omniverse
# -----------------------------------------------------------------------------------


"""
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
How does it look like beating now?
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
"""
# -----------------------------------------------------------------------------------
# 1. Adjust settings and run mechanical simulation
# -----------------------------------------------------------------------------------
simulator.settings.mechanics.analysis.end_time = Quantity(800, "ms")
simulator.settings.mechanics.analysis.dt_d3plot = Quantity(10, "ms")
simulator.settings.electrophysiology.analysis.solvertype="ReactionEikonal"
simulator.dyna_settings.num_cpus = 10
simulator.model.get_part('Aorta').ep_material= EPMaterial.Insulator()
simulator.model.get_part('Pulmonary artery').ep_material= EPMaterial.Insulator()
simulator.simulate()

# -----------------------------------------------------------------------------------
# 2. Animation of beating heart + electrical wave to be performed in omniverse
# -----------------------------------------------------------------------------------

"""
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
How can I know if Mrs Jones has a healthy heart?
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
"""
# -----------------------------------------------------------------------------------
# 1. Read results
# -----------------------------------------------------------------------------------

import pathlib

import ansys.heart.core.models as models
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from ansys.heart.postprocessor.system_model_post import SystemModelPost

model: models.LeftVentricle = models.HeartModel.load_model(path_to_model)

# set simulation path
meca_folder = pathlib.Path(r"D:\REPOS\pyheart\pyansys-heart\downloads\Rodero2021\01\FullHeart\ep-mechanics-eikonal\ep_meca")

# -----------------------------------------------------------------------------------
# 2. Compute pressure volume diagram and plot
# -----------------------------------------------------------------------------------
# Pressure-volume loop figure is an important metric for heart function
system = SystemModelPost(meca_folder)
fig = system.plot_pv_loop()
plt.show()


"""
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
"What would it look like if Mrs Jones has a bundle branch block?"
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
"""
import copy

# -----------------------------------------------------------------------------------
# 1. Create bundle branch block
# -----------------------------------------------------------------------------------
# get bundle branch beamnetwork
simulator_lbbb: EPMechanicsSimulator = copy.deepcopy(simulator)
for beamnet in simulator_lbbb.model.beam_network:
    if 'Left bundle branch' in beamnet.name:
        lbb = beamnet
# Insulate the left bundle branch electrically
lbb.ep_material=EPMaterial.Insulator()

# -----------------------------------------------------------------------------------
# 2. Run simulation
# -----------------------------------------------------------------------------------
simulator_lbbb.simulate()

# ----------------------------------------------------------------------------------------
# 3. animations of healthy and left bundle branch block displayed side by side in Omniverse
# -----------------------------------------------------------------------------------------

"""
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
How can you tell there is a problem?
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
"""
# ----------------------------------------------------------------------------------------
# Animation in Omniverse highlighting a delay between left and right ventricles activation, 
# compared to 'healthy' case
# -----------------------------------------------------------------------------------------
"""
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
What if Mrs Jones had a pacemaker?
------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
"""
from ansys.heart.simulator.settings.settings import Stimulation

# -----------------------------------------------------------------------------------
# 1. Create a cardiac stimulation to simulate pacemaker effect
# -----------------------------------------------------------------------------------

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
simulator_lbbb.settings.electrophysiology.stimulation["stim_apex"] = stim_apex


# -----------------------------------------------------------------------------------
# 2. Run simulation
# -----------------------------------------------------------------------------------

simulator_lbbb.simulate()

# -----------------------------------------------------------------------------------
# 3. Rendering in Omniverse
# -----------------------------------------------------------------------------------
