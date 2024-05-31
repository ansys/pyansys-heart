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

import os

from ansys.heart.preprocessor.models.v0_2.models import HeartModel

# environment variable indicating what model version to use.
# if the preprocessor generated a model with v0.2, then
# this env variable needs to be set explicitly.
os.environ["ANSYS_HEART_MODEL_VERSION"] = "v0.2"

# import relevant classes
from ansys.heart.simulator.settings.settings import DynaSettings
from ansys.heart.simulator.simulator import EPSimulator

# set path to ls-dyna, model and working directory.
lsdyna_path = r"D:\Fortran\intelMPI\27MAY24\mppdyna_27MAY24"
path_to_model = (
    r"D:\REPOS\pyansys-heart\downloads\Strocchi2020\01\FourChamber_02\heart_model.pickle"
)
workdir = r"D:\REPOS\pyansys-heart\downloads\Strocchi2020\01\FourChamber_02"

# initialize settings
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path, dynatype="intelmpi", num_cpus=6, platform="wsl"
)

# load heart model
model = HeartModel.load_model(path_to_model)

# initialize simulator object
simulator = EPSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "ep-only"),
)

# load default simulation settings
simulator.settings.load_defaults()
simulator.compute_uhc()
# compute fiber orientation in the ventricles and atria
import numpy as np

values1 = simulator.model.mesh.point_data["transmural"]
values2 = simulator.model.mesh.point_data["apico-basal"]
values3 = simulator.model.mesh.point_data["rotational"]
littlecloud1 = np.logical_and(values1 >= 0, values1 < 2e-1)
littlecloud2 = np.logical_and(values2 >= 0.74, values2 < 0.76)
littlecloud3 = np.logical_and(values3 >= -2.38, values3 < -2.34)
cloud = np.logical_and(littlecloud1, littlecloud2)
cloud = np.logical_and(cloud, littlecloud3)
pointid = np.nonzero(cloud)[0][0]
# import pyvista as pv

# mycloud = pv.PolyData(simulator.model.mesh.points[pointid, :])

# plotter = pv.Plotter()
# plotter.add_mesh(mycloud, color="maroon", point_size=10.0, render_points_as_spheres=True)
# plotter.add_mesh(simulator.model.mesh, scalars="rotational", opacity=0.1)
# plotter.show()
simulator.compute_fibers()

# Compute the conduction system
simulator.compute_purkinje()
# simulator.model.plot_purkinje()
# simulator.model.dump_model(os.path.join(workdir, "heart_fib_purk.pickle"))
simulator.compute_conduction_system()
# simulator.model.plot_purkinje()


###############################################################################
# Start EP simulation
# ~~~~~~~~~~~~~~~~~~~

# store model with fibers and conduction system
# simulator.model.dump_model(os.path.join(workdir, "heart_fib_beam.pickle"))

# start main simulation
# simulator.settings.electrophysiology.analysis.solvertype = "Eikonal"
simulator.simulate()
