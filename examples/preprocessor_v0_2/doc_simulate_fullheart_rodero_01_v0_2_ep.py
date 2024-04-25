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
lsdyna_path = "mppdyna_d_sse2_linux86_64_intelmmpi_105630"
path_to_model = r"pyansys-heart\downloads\Rodero2021\01\FullHeart_v0_2\heart_model.pickle"
workdir = "pyansys-heart\downloads\Rodero2021\01\FullHeart_v0_2"

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

# compute fiber orientation in the ventricles and atria
simulator.compute_fibers()

# Compute the conduction system
simulator.compute_purkinje()
# simulator.model.plot_purkinje()
simulator.compute_conduction_system()
simulator.model.plot_purkinje()


###############################################################################
# Start EP simulation
# ~~~~~~~~~~~~~~~~~~~

# store model with fibers and conduction system
simulator.model.dump_model(os.path.join(workdir, "heart_fib_beam.pickle"))

# start main simulation
simulator.simulate()
