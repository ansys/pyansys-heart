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
import pytest
import pyvista as pv
import pyvista.examples as examples

from ansys.health.heart.utils.vtk_utils import (
    are_connected,
    cell_ids_inside_enclosed_surface,
    find_cells_close_to_nodes,
    find_corresponding_points,
    generate_thickness_lines,
)


def test_find_corresponding_points():
    """Test correspoing points searching."""
    sphere1 = pv.Sphere(radius=1, center=(0, 0, 0))
    sphere2 = pv.Sphere(radius=1.2, center=(0, 0, 0))
    res = find_corresponding_points(sphere1, sphere2)
    assert res.shape == (2, sphere1.GetNumberOfPoints())
    assert res[1, 0] == 0
    assert res[1, 841] == 841


def test_generate_thickness_lines():
    """Test thcikness lines generation."""
    sphere1 = pv.Sphere(radius=1, center=(0, 0, 0))
    sphere2 = pv.Sphere(radius=1.2, center=(0, 0, 0))
    lines = generate_thickness_lines(sphere1, sphere2)
    assert lines.GetNumberOfCells() == sphere1.GetNumberOfPoints()
    assert pytest.approx(0.2, rel=0.01) == lines["thickness"][0]


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

    surface = pv.Box((-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)).triangulate()

    # two points. One inside, one outside.
    points = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    points = pv.PolyData(points)
    pt_ids_inside = cell_ids_inside_enclosed_surface(points, surface)

    assert len(pt_ids_inside) == 1
    assert pt_ids_inside == [0]

    # generate random points with seed for reproducibility
    rng = np.random.default_rng(seed=42)
    points = rng.random((10000, 3)) * 2 - 1
    points = pv.PolyData(points)

    pt_ids_inside = cell_ids_inside_enclosed_surface(points, surface)
    assert len(pt_ids_inside) == 1244
