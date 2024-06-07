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

from ansys.heart.preprocessor.mesh.objects_2 import Mesh
import numpy as np
import pytest  # noqa F401
import pyvista as pv  # noqa F401

from tests.heart.conftest import download_asset, get_assets_folder  # noqa F401

skip_test = os.name != "nt"


def test_mesh_object_assign():
    """Test the mesh object with a mixed tet, triangle, line mesh."""
    # construct a mixed grid with:
    # 1x Tetrahedron
    # 4x triangles
    # 3x lines
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])

    # for each face in the tet
    tets = [4, 0, 1, 2, 3]
    triangles = [3, 0, 1, 3, 3, 1, 2, 3, 3, 2, 0, 3, 3, 0, 2, 1]
    lines = [2, 0, 1, 2, 1, 2, 2, 2, 3]
    cell_types = [pv.CellType.TETRA] + [pv.CellType.TRIANGLE] * 4 + [pv.CellType.LINE] * 3
    cells = tets + triangles + lines

    grid = pv.UnstructuredGrid(cells, cell_types, points)

    grid1 = Mesh(cells, cell_types, points)

    mask = grid.celltypes == pv.CellType.LINE
    lines = grid.cells_dict[pv.CellType.LINE]
    lines = grid.extract_cells(mask)

    grid1.cell_data["cdata_tet1"] = 1
    grid1.point_data["pdata_tet1"] = 2

    # TODO
    # add lines
    lines = pv.lines_from_points([points[0, :], [-1, 0, 0], [-2, 0, 0]])
    # create some dummy data
    lines["cdata_lines1"] = np.array([10, 11], dtype=np.int32)
    lines["cdata_lines2"] = np.array([10.5, 11.5], dtype=np.float64)
    lines["pdata_lines1"] = np.array([[100, 101], [102, 103], [104, 105]], dtype=np.float64)

    grid2 = grid1.add_mesh(lines)

    # check array names
    for name in ["cdata_lines1", "cdata_lines2", "cdata_tet1", "pdata_tet1", "pdata_lines1"]:
        assert name in grid2.array_names, f"Array with name {name} does  not exist"

    # check data types
    assert isinstance(grid2.cell_data["cdata_lines1"][0], np.int32)
    assert isinstance(grid2.cell_data["cdata_lines2"][0], np.float64)

    # check shape of array
    assert grid2["pdata_lines1"].shape[1] == lines["pdata_lines1"].shape[1]

    grid

    assert grid.__add__(lines) == (grid + lines)

    pass
