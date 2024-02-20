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
todo
"""

import os
import subprocess

from ansys.heart.calibration.active_calibration import ActiveCalibration

if __name__ == "__main__":
    lsdyna_path = r"D:\my_path_to_ls_dyna\lsdyna_executable.exe"
    ncpu = 2
    dynatype = "smp"  # depend on your lsdyna_path

    # make sure pickle model contains right fiber information
    # i.e. fibergeneration has been run
    mdoel_path = r"path_to_model.pickle"

    # Create a folder that contains a LSOPT project for passive calibration
    target_folder = r"my_targe_folder"
    ActiveCalibration.create_calibration_project(
        target_folder, lsdyna_path, ncpu, dynatype, mdoel_path
    )

    # Run calibration
    os.chdir(target_folder)
    subprocess.run(["my_path_to_ls_opt", "ActiveCalibration.lsopt"])
