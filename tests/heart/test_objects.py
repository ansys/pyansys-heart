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

import numpy as np
import pyvista as pv

from ansys.health.heart.objects import _convert_int64_to_int32


def test_convert_int64_to_int32():
    """Test the `_change_int64_to_int32` function."""
    # Create a sample PyVista mesh with int64 arrays
    mesh = pv.Sphere()
    mesh["cell_data"] = np.arange(0, mesh.n_cells, dtype=np.int64)  # Add int64 cell data
    mesh["point_data"] = np.arange(0, mesh.n_points, dtype=np.int64)  # Add int64 point data

    mesh["data"] = np.arange(0, mesh.n_points, dtype=np.int64)  # Add int64 point data
    mesh["data"] = np.arange(0, mesh.n_cells, dtype=np.int64)  # Add int64 point data

    # Verify initial data types are int64
    assert mesh["cell_data"].dtype == np.int64
    assert mesh["point_data"].dtype == np.int64

    # Call the function to convert int64 to int32
    _convert_int64_to_int32(mesh)

    # Verify that the data types have been converted to int32
    assert mesh["cell_data"].dtype == np.int32
    assert mesh["point_data"].dtype == np.int32
    assert mesh.cell_data["data"].dtype == np.int32
    assert mesh.point_data["data"].dtype == np.int32
