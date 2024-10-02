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

from ansys.heart.core.models import HeartModel
from ansys.heart.postprocessor.auto_process import mech_post, zerop_post

if __name__ == "__main__":
    """
    Example show how to post-process mechanical simulation results
    """
    # Get heart model
    model = HeartModel.load_model(".")
    # Define directory where zeropressure simulation is carried out
    # a folder "post" will be created under it with key simulation results (json, png, vtk...)
    # use Paraview state file post_zerop2.pvsm, you can easily visualize generated results
    dir = r""
    zerop_post(dir, model)

    # Define directory where main mechanical simulation is carried out
    # a folder "post" will be created under it with key simulation results (json, png, vtk...)
    # use Paraview state file post_main2.pvsm, you can easily visualize generated results
    dir2 = r""
    mech_post(dir2, model)
