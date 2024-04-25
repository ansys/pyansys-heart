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
This example creates an LSOPT project that is used to calibrate passive material parameters.

The project is saved in PassiveCalibration.lsopt, a very based setting is defined:
Radial Basis function is selected as the metamodel and (only) 3 simulation points are defined.
You can refine them in LSOPT UI.

The input parameters are 'ps1' and 'ps2' which are weights before parameter k1/k1f and k2/k2f.
They are defined in file 'parameters.k' and their ranges are defined in LSOPT file.
The inputs can be redefined by rewriting the method apply_input_parameter().

The objective function is defined by the RSM error between simulated volume and Klotz curve
at different pressure loads. The error is saved as a scalr in a file named 'result'.
The objective function can be changed by rewriting the method define_objective().

Notes
-----
    Rewriting parameter/objective definition may imply necessary changes in LSOPT file.

References
----------
    doi: 10.3389/fphys.2018.00539

    doi:10.1093/europace/euw369

    http://dx.doi.org/10.1016/j.jcp.2012.09.015

"""

import os
import subprocess

from ansys.heart.calibration.passive_calibration import PassiveCalibration

if __name__ == "__main__":
    lsdyna_path = r"D:\my_path_to_ls_dyna\lsdyna_executable.exe"
    ncpu = 2
    dynatype = "smp"  # depend on your lsdyna_path

    # make sure pickle model contains right fiber information
    # i.e. fibergeneration has been run
    mdoel_path = r"path_to_model.pickle"

    # Create a folder that contains a LSOPT project for passive calibration
    target_folder = r"my_targe_folder"
    PassiveCalibration.create_calibration_project(
        target_folder, lsdyna_path, ncpu, dynatype, mdoel_path
    )

    # Run calibration
    os.chdir(target_folder)
    subprocess.run(["my_path_to_ls_opt", "PassiveCalibration.lsopt"])
