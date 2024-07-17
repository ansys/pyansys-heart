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

Four-chamber EP-simulator example
---------------------------------
This example shows you how to consume a four-cavity heart model and
set it up for the main electropysiology simulation. This examples demonstrates how
you can load a pre-computed heart model, compute the fiber direction, compute the
purkinje network and conduction system and finally simulate the electrophysiology.
"""

###############################################################################
# Example setup
# -------------
# before computing the fiber orientation, purkinje network we need to load
# the required modules, load a heart model and set up the simulator.
#
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory, model, and ls-dyna executable.

# sphinx_gallery_start_ignore
# Note that we need to put the thumbnail here to avoid weird rendering in the html page.
# sphinx_gallery_thumbnail_path = '_static/images/purkinje.png'
# sphinx_gallery_end_ignore

import os

# set this environment variable to ensure you are using v0.2 of the model
os.environ["ANSYS_HEART_MODEL_VERSION"] = "v0.2"
from pint import Quantity
from ansys.heart.preprocessor.mesh.objects import Point
import ansys.heart.preprocessor.models.v0_2.models as models
from ansys.heart.simulator.simulator import DynaSettings, EPSimulator


# accept dpf license aggrement
# https://dpf.docs.pyansys.com/version/stable/getting_started/licensing.html#ref-licensing
os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"

# set working directory and path to model.
workdir = r"D:\REPOS\pyansys-heart\downloads\Strocchi2020\01\FourChamber"

path_to_model = os.path.join(workdir, "heart_model.pickle")

# specify LS-DYNA path (last tested working versions is intelmpi-linux-DEV-106117)
lsdyna_path = r"D:\Fortran\intelMPI\10JUL24\mppdyna_10JUL24"

# load four chamber heart model.
model: models.FourChamber = models.HeartModel.load_model(path_to_model)

# Define electrode positions and add them to model (correspond to patient 01 only)
# Positions were defined using a template torso geometry.
electrodes = [
    Point(name="V1", xyz=[-29.893285751342773, 27.112899780273438, 373.30865478515625]),
    Point(name="V2", xyz=[33.68170928955078, 30.09606170654297, 380.5427551269531]),
    Point(name="V3", xyz=[56.33562469482422, 29.499839782714844, 355.533935546875]),
    Point(name="V4", xyz=[100.25729370117188, 43.61333465576172, 331.07635498046875]),
    Point(name="V5", xyz=[140.29800415039062, 81.36004638671875, 349.69970703125]),
    Point(name="V6", xyz=[167.9899139404297, 135.89862060546875, 366.18634033203125]),
    Point(name="RA", xyz=[133.84518432617188, 101.44053649902344, 534.9176635742188]),
    Point(name="LA", xyz=[-176.06332397460938, 57.632076263427734, 509.14202880859375]),
    Point(name="RL", xyz=[203.38825799615842, 56.19020893502452, 538.5052677637375]),
    Point(name="LL", xyz=[128.9441375732422, 92.85327911376953, 173.07363891601562]),
]
import numpy as np

electrodesXYZ = np.array(
    [
        [-29.893285751342773, 27.112899780273438, 373.30865478515625],
        [33.68170928955078, 30.09606170654297, 380.5427551269531],
        [56.33562469482422, 29.499839782714844, 355.533935546875],
        [100.25729370117188, 43.61333465576172, 331.07635498046875],
        [140.29800415039062, 81.36004638671875, 349.69970703125],
        [167.9899139404297, 135.89862060546875, 366.18634033203125],
        [133.84518432617188, 101.44053649902344, 534.9176635742188],
        [-176.06332397460938, 57.632076263427734, 509.14202880859375],
        [203.38825799615842, 56.19020893502452, 538.5052677637375],
        [128.9441375732422, 92.85327911376953, 173.07363891601562],
    ]
)


elecV1 = np.array([-29.893285751342773, 27.112899780273438, 373.30865478515625])
elecV2 = np.array([33.68170928955078, 30.09606170654297, 380.5427551269531])
elecV3 = np.array([56.33562469482422, 29.499839782714844, 355.533935546875])
elecV4 = np.array([100.25729370117188, 43.61333465576172, 331.07635498046875])
elecV5 = np.array([140.29800415039062, 81.36004638671875, 349.69970703125])
elecV6 = np.array([167.9899139404297, 135.89862060546875, 366.18634033203125])
elecLA = np.array([-176.06332397460938, 57.632076263427734, 509.14202880859375])
elecRA = np.array([133.84518432617188, 101.44053649902344, 534.9176635742188])
elecRL = np.array([203.38825799615842, 56.19020893502452, 538.5052677637375])
elecLL = np.array([128.9441375732422, 92.85327911376953, 173.07363891601562])


model.electrodes = electrodes


plot_electrodes = True
if plot_electrodes:
    import pyvista as pv

    longvecRot = np.array([-4.6643, 110.254, 380.397]) - np.array([22.0392, 52.1924, 351.721])
    frontalvecRot = np.array([71.3322, 73.256, 352.189]) - np.array([-66.662, 116.456, 369.011])
    plotter = pv.Plotter()
    electrodes_cloud = pv.PolyData(electrodesXYZ)
    plotter.add_mesh(model.mesh)

    # plotter.add_mesh(pv.PolyData(elecV1), color="blue", point_size=10.0)
    # plotter.add_mesh(pv.PolyData(elecV2), color="blue", point_size=10.0)
    # plotter.add_mesh(pv.PolyData(elecV3), color="blue", point_size=10.0)
    # plotter.add_mesh(pv.PolyData(elecV4), color="blue", point_size=10.0)
    # plotter.add_mesh(pv.PolyData(elecV5), color="blue", point_size=10.0)
    # plotter.add_mesh(pv.PolyData(elecV6), color="blue", point_size=10.0)
    # plotter.add_mesh(pv.PolyData(elecLA), color="blue", point_size=10.0)
    # plotter.add_mesh(pv.PolyData(elecRA), color="green", point_size=10.0)
    # plotter.add_mesh(pv.PolyData(elecRL), color="yellow", point_size=10.0)
    # plotter.add_mesh(pv.PolyData(elecLL), color="black", point_size=10.0)
    plotter.add_mesh(electrodes_cloud, color="blue", point_size=10.0)

    electrodes_longitudinal_m10 = electrodes_cloud.rotate_vector(
        longvecRot, -10, (8.66776, 118.739, 350.017), inplace=False
    )
    electrodes_longitudinal_m5 = electrodes_cloud.rotate_vector(
        longvecRot, -5, (8.66776, 118.739, 350.017), inplace=False
    )
    electrodes_longitudinal_p5 = electrodes_cloud.rotate_vector(
        longvecRot, 5, (8.66776, 118.739, 350.017), inplace=False
    )
    electrodes_longitudinal_p10 = electrodes_cloud.rotate_vector(
        longvecRot, 10, (8.66776, 118.739, 350.017), inplace=False
    )

    electrodes_frontal_m10 = electrodes_cloud.rotate_vector(
        frontalvecRot, -10, (8.66776, 118.739, 350.017), inplace=False
    )
    electrodes_frontal_m5 = electrodes_cloud.rotate_vector(
        frontalvecRot, -5, (8.66776, 118.739, 350.017), inplace=False
    )
    electrodes_frontal_p5 = electrodes_cloud.rotate_vector(
        frontalvecRot, 5, (8.66776, 118.739, 350.017), inplace=False
    )
    electrodes_frontal_p10 = electrodes_cloud.rotate_vector(
        frontalvecRot, 10, (8.66776, 118.739, 350.017), inplace=False
    )
    # plotter.add_mesh(electrodes_cloud_minus10, color="maroon", point_size=10.0)
    plotter.show()

if not isinstance(model, models.FourChamber):
    raise TypeError("Expecting a FourChamber heart model.")

# set base working directory
model.info.workdir = str(workdir)

###############################################################################
# Instantiate the simulator object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# instantiate the simulator and settings appropriately.

# instantaiate dyna settings of choice
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path,
    dynatype="intelmpi",
    platform="wsl",
    num_cpus=1,
)

# instantiate simulator. Change options where necessary.
simulator = EPSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "simulation-EP-origin"),
)

###############################################################################
# Load simulation settings
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Here we load the default settings.

simulator.settings.load_defaults()

###############################################################################
# Compute Universal Ventricular Coordinates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The transmural coordinate is used to define the endo, mid and epi layers.

###############################################################################
simulator.compute_uhc()

###############################################################################
# Compute the fiber orientation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute fiber orientation and plot the computed fibers on the entire model.

###############################################################################
# .. warning::
#    Atrial fiber orientation is approximated by apex-base direction, the development is undergoing.

simulator.compute_fibers()
# simulator.model.plot_fibers(n_seed_points=2000)

###############################################################################
# .. image:: /_static/images/fibers.png
#   :width: 400pt
#   :align: center

###############################################################################
# Compute conduction system
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute conduction system and purkinje network and visualize.
# The action potential will propagate faster through this system
# compared to the rest of the model.
simulator.settings.purkinje.pmjtype = Quantity(4)
simulator.settings.purkinje.pmjradius = Quantity(2)

simulator.compute_purkinje()

# by calling this method, stimulation will at Atrioventricular node
# if you skip it, stimulation will at apex nodes of two ventricles
simulator.compute_conduction_system()

# simulator.model.plot_purkinje()

###############################################################################
# .. image:: /_static/images/purkinje.png
#   :width: 400pt
#   :align: center

###############################################################################
# Start main simulation
# ~~~~~~~~~~~~~~~~~~~~~
# Start the main EP simulation. This uses the previously computed fiber orientation
# and purkinje network to set up and run the LS-DYNA model using different solver
# options

# simulator.simulate()
# # The two following solves only work with LS-DYNA DEV-110013 or later
# simulator.settings.electrophysiology.analysis.solvertype = "Eikonal"
# simulator.simulate(folder_name="main-ep-Eikonal")
simulator.settings.electrophysiology.analysis.solvertype = "ReactionEikonal"
simulator.settings.electrophysiology.analysis.end_time = Quantity(500, "ms")
simulator.settings.electrophysiology.analysis.dt_d3plot = Quantity(5, "ms")

# In Purkinje fibers (2–3 m/s) than in myocardial cells (0.7fiber –0.2sheets m/s)
simulator.settings.electrophysiology.material.beam["velocity"] = Quantity(2.5, "ms")

# simulator.simulate(folder_name="main-ep-ReactionEikonal")


###############################################################################
# We can plot transmembrane potential in LS-PrePost

###############################################################################
# .. video:: ../../_static/images/ep_4cv.mp4
#   :width: 600
#   :loop:
#   :class: center
