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

In order to obtain valid ECGs, a 1.5 mm mesh is required :
model.mesh_volume(use_wrapper=True, global_mesh_size=2.0, _global_wrap_size=2.0)
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

import json
import os
from pathlib import Path

import numpy as np
import pyvista as pv

import ansys.health.heart.models as models
from ansys.health.heart.objects import SurfaceMesh
from ansys.health.heart.pre.conduction_path import ConductionPath, ConductionPathType
import ansys.health.heart.settings.material.cell_models as cell_mod
import ansys.health.heart.settings.material.ep_material as mat
from ansys.health.heart.simulator import DynaSettings, EPSimulator

os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"

os.environ["ANSYSLMD_LICENSE_FILE"] = "1055@SERVER"

############################### à modifier  ###################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Instantiate the simulator and define settings.

# .. note::
#    The ``DynaSettings`` object supports several LS-DYNA versions and platforms,
#    including ``smp``, ``intempi``, ``msmpi``, ``windows``, ``linux``, and ``wsl``.
#    Choose the one that works for your setup.
# Specify the LS-DYNA path. (The last tested working version is ``intelmpi-linux-DEV-106117``.)
ls_path = Path("C:\\") / "Users" / "mbouchai" / "pyansys-heart" / "ansys-dev-env" / "exe_07_07"
lsdyna_path = os.path.join(
    ls_path, r"ls-dyna_mpp_d_DEV_122687-g00479ad686_winx64_ifort190_sse2_msmpi.exe"
)

# Instantiate LS-DYNA settings.
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path, dynatype="msmpi", num_cpus=18, platform="windows"
)

# Plot option for all scripts
plot = True
compute_fib = False
compute_purk = False
# A .json file containing the information required to run the script.
# Rodero 01, 03, and 04 are available. To customize a new model, you
# must either add the information to the .json file using the same format,
#  or add the customization directly in this script.

# path of data model
ecg_data_path = (
    Path("D:\\")
    / "env_pyheart_meca"
    / "pyansys-heart"
    / "src"
    / "ansys"
    / "health"
    / "heart"
    / "data_examples"
)
rodero_file_info = os.path.join(
    ecg_data_path,
    "info_rodero_conduction_system_ecg.json",
)

# choice which model you want
num_rodero = 1

# Set the working directory and path to the model. This example assumes that there is a
workdir = (
    Path("D:\\")
    / "env_pyheart_meca"
    / "env_pyheart_meca"
    / "mes_modeles"
    / "EP"
    / "Rodero2021"
    / ("0" + str(num_rodero))
    / "FullHeart"
)

simulation_directory = os.path.join(workdir, "simulation-EP")

##############################################################################
# Import rodero model parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~
with open(rodero_file_info, "r") as f:
    part_info: dict = json.load(f)
info_rodero = part_info["Rodero_" + str(num_rodero)]

###############################################################################
# Load the full-heart model
# ~~~~~~~~~~~~~~~~~~~~~~~~
path_to_model = os.path.join(workdir, "heart_model.vtu")
path_to_partinfo = os.path.join(workdir, "heart_model.partinfo.json")

model: models.FourChamber = models.HeartModel.load_model(
    path_to_model, path_to_partinfo, working_directory=workdir
)

# Instantiate the simulator, modifying options as necessary.
simulator = EPSimulator(
    model=model, dyna_settings=dyna_settings, simulation_directory=simulation_directory
)

# Save the model.
model.mesh.save(os.path.join(model.workdir, "simulation_model.vtu"))

new_electrodes = simulator.model.set_electrodes(angle1=10, angle2=-10)

if plot:
    heart_model = pv.read(str(workdir / "heart_model.vtu"))
    plotter = pv.Plotter()
    point_cloud_1 = pv.PolyData(new_electrodes)
    plotter.add_mesh(point_cloud_1, color="blue", point_size=10, render_points_as_spheres=True)
    plotter.add_mesh(heart_model, color="red", show_edges=False, opacity=0.8)
    plotter.show()


###############################################################################
# Load simulation settings
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Load the default settings.
simulator.settings.load_defaults()
simulator.settings.electrophysiology.analysis.solvertype = "ReactionEikonal"

###############################################################-################
# Compute fiber orientation
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute fiber orientation and plot the fibers on the entire model.
if compute_fib:
    # Import the appendage landmarks.
    from ansys.health.heart.pre.database_utils import right_atrium_appendage_landmarks

    # Get the right atrium appendage landmark of the first case of Rodero2021.
    right_atrium_appendage_coordinates = right_atrium_appendage_landmarks.get("Rodero2021").get(
        num_rodero
    )

    # Compute ventricular fibers.
    simulator.compute_fibers()

    # Compute atrial fibers.
    simulator.model.right_atrium.active = True
    simulator.model.left_atrium.active = True
    simulator.model.right_atrium.fiber = True
    simulator.model.left_atrium.fiber = True
    simulator.compute_left_atrial_fiber()
    simulator.compute_right_atrial_fiber(appendage=right_atrium_appendage_coordinates)
    if plot:
        simulator.model.plot_fibers(n_seed_points=1000)

###############################################################################
# Compute the conduction system
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
# Compute the conduction system and Purkinje network, and then visualize the results.
# The action potential propagates faster through this system compared to the rest of the model.


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
# Function to merge the right ventricular endocardium and the right side of the Purkinje network
# so that the surface is treated as a single, complete surface during network generation.
def merge_surfaces(surf1: SurfaceMesh, surf2: SurfaceMesh, keep="first") -> SurfaceMesh:
    merged = surf1.merge(surf2, merge_points=True, tolerance=1e-8)
    new_surf = SurfaceMesh(merged)
    new_surf.nodes = np.array(new_surf.points)
    faces = np.reshape(new_surf.faces, (new_surf.n_cells, 4))[:, 1:]
    new_surf.triangles = faces
    if keep == "first":
        new_surf.id = surf1.id
        new_surf.name = surf1.name
    else:
        new_surf.id = surf2.id
        new_surf.name = surf2.name
    return new_surf


simulator.model.right_ventricle.endocardium = merge_surfaces(
    simulator.model.right_ventricle.endocardium,
    simulator.model.right_ventricle.septum,
    keep="first",
)

directory = os.path.join(simulator.root_directory, "purkinjegeneration")
if compute_purk:
    orig_num_cpus = simulator.dyna_settings.num_cpus
    simulator.dyna_settings.num_cpus = 1

    ureg = simulator.settings.purkinje.pmjtype._REGISTRY
    simulator.settings.purkinje.pmjtype = 1 * ureg.dimensionless
    simulator.settings.purkinje.edgelen = 0.5 * ureg.dimensionless
    simulator.settings.purkinje.nsplit = 6 * ureg.dimensionless
    simulator.settings.purkinje.ngen = 500 * ureg.dimensionless
    simulator.settings.purkinje.nbrinit = 6 * ureg.dimensionless

    simulator._write_purkinje_files(directory)
    input_file = os.path.join(directory, "main.k")
    simulator._run_dyna(input_file)
    simulator.dyna_settings.num_cpus = orig_num_cpus

if plot:
    # Visualize the entire conduction system
    simulator.model.plot_purkinje()

# By calling this method, stimulation is at the atrioventricular node.
# If you do not call this method, the two apex regions of the ventricles are stimulated.
###############################################################################
# Compute the conduction system
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import ansys.health.heart.models_utils as heart_model_utils

beam_list, simulator.model._landmarks = heart_model_utils.define_full_conduction_system(
    simulator.model, os.path.join(simulator.root_directory, "purkinjegeneration")
)

[
    left_purkinje,
    right_purkinje,
    sa_av,
    his_top,
    his_left,
    his_right,
    left_bundle,
    right_bundle,
] = beam_list

# model.assign_conduction_paths(beam_list)

# Compute the conduction system user add
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The coordinates of the fascicle terminations, or those that are imposed,
# come from a file called ‘Rodero info’. They can be modified directly in
# this file for other models. The correct definition of these coordinates
# is responsible for physiological activation of the cardiac electrical activity
# and allows obtaining valid ECGs. These points need to be redefined for each new model.

# Create the Bachmann bundle
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bachman_bundle = ConductionPath.create_from_keypoints(
    name=ConductionPathType.BACHMANN_BUNDLE,
    keypoints=[
        simulator.model._landmarks.sa_node.xyz,
        info_rodero["bachman_point_0"],
        info_rodero["bachman_point_1"],
        info_rodero["bachman_point_2"],
    ],
    id=9,
    base_mesh=model.left_atrium.endocardium,
    line_length=None,
    center=True,
)
bachman_bundle.add_pmj_path(list(range(1, bachman_bundle.mesh.n_points - 1, 4)))
bachman_bundle.up_path = sa_av

# The mid and post SA-AV node conduction paths are created by providing a list of keypoints.
mid_sa_av = ConductionPath.create_from_keypoints(
    name=ConductionPathType.MID_SAN_AVN,
    keypoints=[
        simulator.model._landmarks.sa_node.xyz,
        info_rodero["mid_sa_av_point_0"],
        info_rodero["mid_sa_av_point_1"],
        info_rodero["mid_sa_av_point_2"],
        simulator.model._landmarks.av_node.xyz,
    ],
    id=10,
    base_mesh=model.right_atrium.endocardium,
    line_length=None,
    center=True,
)
mid_sa_av.add_pmj_path(list(range(5, mid_sa_av.mesh.n_points - 5, 4)))
mid_sa_av.up_path = sa_av
mid_sa_av.down_path = sa_av

post_sa_av = ConductionPath.create_from_keypoints(
    name=ConductionPathType.POST_SAN_AVN,
    keypoints=[
        simulator.model._landmarks.sa_node.xyz,
        info_rodero["post_sa_av_point_0"],
        info_rodero["post_sa_av_point_1"],
        info_rodero["post_sa_av_point_2"],
        simulator.model._landmarks.av_node.xyz,
    ],
    id=11,
    base_mesh=model.right_atrium.endocardium,
    line_length=None,
    center=True,
)
post_sa_av.add_pmj_path(list(range(5, post_sa_av.mesh.n_points - 5, 4)))
post_sa_av.up_path = sa_av
post_sa_av.down_path = sa_av

###############################################################################
# Create the Fascicle
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

left_anterio_fascile = ConductionPath.create_from_keypoints(
    name=ConductionPathType.LEFT_ANTERIOR_FASCILE,
    keypoints=[
        simulator.model._landmarks.his_left_end_node.xyz,
        info_rodero["left_anterio_fascile_point"],
    ],
    id=12,
    base_mesh=model.left_ventricle.endocardium,
    connection=None,
    line_length=None,
)
left_anterio_fascile.up_path = his_left
left_anterio_fascile.down_path = left_purkinje

###############################################################################
# Create the Fascicle
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

left_posterior_fascile = ConductionPath.create_from_keypoints(
    name=ConductionPathType.LEFT_POSTERIOR_FASCICLE,
    keypoints=[
        simulator.model._landmarks.his_left_end_node.xyz,
        info_rodero["left_posterior_fascile_point"],
    ],
    id=13,
    base_mesh=model.left_ventricle.endocardium,
    connection=None,
    line_length=None,
)
left_posterior_fascile.up_path = his_left
left_posterior_fascile.down_path = left_purkinje

# Addition of new conduction fibers to the model
model.assign_conduction_paths(
    [
        left_purkinje,
        right_purkinje,
        sa_av,
        his_top,
        his_left,
        his_right,
        left_bundle,
        right_bundle,
        bachman_bundle,
        mid_sa_av,
        post_sa_av,
        left_anterio_fascile,
        left_posterior_fascile,
    ]
)


if plot:
    # Visualize the entire conduction system
    simulator.model.plot_purkinje()

###############################################################################
# Setup of custom parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Conduction path parameterization
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# atrium conduction path
mat_sa_av = mat.ActiveBeam()
mat_sa_av.sigma_fiber = 2.25
simulator.model.conduction_paths[2].ep_material = mat_sa_av
simulator.model.conduction_paths[8].ep_material = mat_sa_av
simulator.model.conduction_paths[9].ep_material = mat_sa_av
simulator.model.conduction_paths[10].ep_material = mat_sa_av

# hit_top
# The His_top bundle conductivity parameters must be implemented
# to coordinate atrial and ventricular activity according to the
# chosen action potential durations.
mat_his_top = mat.ActiveBeam()
mat_his_top.sigma_fiber = info_rodero["Cv_his_top"]
simulator.model.conduction_paths[3].ep_material = mat_his_top

# ventricle conduction path
mat_purkinje = mat.ActiveBeam()
mat_purkinje.sigma_fiber = 2.0

# his_left & his_right
simulator.model.conduction_paths[4].ep_material = mat_purkinje
simulator.model.conduction_paths[5].ep_material = mat_purkinje

# left_bundle & right_bundle
simulator.model.conduction_paths[6].ep_material = mat_purkinje
simulator.model.conduction_paths[7].ep_material = mat_purkinje
simulator.model.conduction_paths[0].ep_material = mat_purkinje
simulator.model.conduction_paths[1].ep_material = mat_purkinje

# fascicle left ventricle
# The fascicle conductivity must be set so that both fascicles activate
# the tissue at the same time. Therefore, the travel time variable allows
# synchronizing activation of the free wall with that of the septum.

# Set the CV of the left and right bundle branches
simulator.model.conduction_paths[6].ep_material.sigma_fiber = 1.786
simulator.model.conduction_paths[7].ep_material.sigma_fiber = 1.786
# Set the CV of the left and right purkinjes
simulator.model.conduction_paths[0].ep_material.sigma_fiber = 2
simulator.model.conduction_paths[1].ep_material.sigma_fiber = 2

# Tissue parameterization
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ep_vent = mat.Active2()
ep_vent.sigma_fiber = 0.7
ep_vent.sigma_sheet = 0.35
ep_vent.sigma_sheet_normal = 0.18
ep_vent.cond_sigma_fiber = 0.17
ep_vent.cond_sigma_sheet = 0.08
ep_vent.cond_sigma_sheet_normal = 0.08
simulator.model.parts[0].ep_material = ep_vent
simulator.model.parts[1].ep_material = ep_vent
simulator.model.parts[2].ep_material = ep_vent

ep_atrium = mat.Active2()
ep_atrium.cell_model = cell_mod.TentusscherAtria()
ep_atrium.sigma_fiber = 1.2
ep_atrium.sigma_sheet = 0.6
ep_atrium.sigma_sheet_normal = 0.3
ep_atrium.cond_sigma_fiber = 0.17
ep_atrium.cond_sigma_sheet = 0.08
ep_atrium.cond_sigma_sheet_normal = 0.08
simulator.model.parts[3].ep_material = ep_atrium
simulator.model.parts[4].ep_material = ep_atrium

simulator.compute_uhc()
# Configurable number of layers
simulator.settings.electrophysiology.nb_layers = 7

###############################################################################
# Start the main simulation
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Start the main electrophysiology simulation. This uses the previously computed fiber orientation
# and Purkinje network to set up and run the LS-DYNA model.

simulator.simulate(folder_name="main-ep-ReactionEikonal")
