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

"""Example to preprocess and run a FourChamber EP simulation."""

import os
import pathlib

import ansys.heart.preprocessor.models.v0_1.models as models
from ansys.heart.simulator.simulator import EPSimulator
from ansys.heart.simulator.support import run_preprocessor

if __name__ == "__main__":
    """FourChamber example.

    1. Extracts simulation mesh
    2. Use simulator class to launch an EP simulation:
            - compute fibers
            - compute purkinje network
            - launch main simulation

    Please change paths according to your workspace
    """
    # extract simulation mesh(es)
    path_to_case = os.path.join(
        pathlib.Path(__file__).parents[3], "downloads\\Strocchi2020\\01\\01.case"
    )
    workdir = os.path.join(pathlib.Path(path_to_case).parent, "FourChamber")

    path_to_model = os.path.join(workdir, "heart_model.pickle")

    use_preprocessor = True

    if use_preprocessor:
        model = run_preprocessor(
            model_type=models.FourChamber,
            database="Strocchi2020",
            path_original_mesh=path_to_case,
            work_directory=workdir,
            path_to_model=path_to_model,
            mesh_size=2.0,
        )

    # specify LS-DYNA path
    lsdyna_path = r"D:\my_path_to_ls_dyna\lsdyna_executable.exe"

    model: models.FourChamber = models.HeartModel.load_model(path_to_model)
    ## instantiate simulator, please change the dynatype accordingly
    simulator = EPSimulator(
        model=model,
        lsdynapath=lsdyna_path,
        dynatype="smp",
        num_cpus=1,
        simulation_directory=os.path.join(model.info.workdir, "simulation-EP"),
        platform="wsl",
    )
    # load default settings
    simulator.settings.load_defaults()
    # compute the fiber orientation
    simulator.compute_fibers()
    # visualize computed fibers
    simulator.model.plot_fibers(n_seed_points=2000)
    # compute purkinje network
    simulator.compute_purkinje()
    # visualize purkinje
    model.dump_model()
    model: models.FourChamber = models.HeartModel.load_model(path_to_model)
    simulator.model.plot_purkinje()

    simulator.compute_conduction_system()
    # start simulation
    simulator.simulate()
    print("done")
