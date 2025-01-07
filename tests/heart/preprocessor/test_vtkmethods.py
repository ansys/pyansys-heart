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
import pyvista.examples as examples

from ansys.heart.core.helpers.vtkmethods import (
    are_connected,
    cell_ids_inside_enclosed_surface,
    find_cells_close_to_nodes,
)


def test_check_if_connected():
    """Test whether two pyvista objects are connected."""
    box1 = pv.Box()
    box2 = pv.Box()

    assert are_connected(box1, box2.translate((2, 0, 0)))
    assert not are_connected(box1, box2.translate((3, 0, 0)))

    return


def test_find_cells_close_to_nodes():
    """Test finding cells close to list of nodes."""
    mesh = examples.load_tetbeam()

    mesh.points[0, :]
    # expecting all nodes to be enclosed.
    assert np.all(find_cells_close_to_nodes(mesh, [0], radius=6) == np.arange(0, mesh.n_cells))
    # expecting just a single enclosed node.
    assert find_cells_close_to_nodes(mesh, [0], radius=0.1) == np.array([0])

    return


def test_cell_ids_inside_enclosed_surface():
    """Test finding cell ids inside a given surface."""
    # box = pv.Box(bounds=(-0.1, 0.1, -0.1, 0.1, -0.1, 1.0))
    mesh = examples.load_tetbeam()
    surface = mesh.extract_surface()

    assert np.all(cell_ids_inside_enclosed_surface(mesh, surface) == np.arange(0, mesh.n_cells))
