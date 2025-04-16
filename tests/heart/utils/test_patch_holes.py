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

import pytest
import pyvista as pv

from ansys.health.heart.utils.vtk_utils import get_patches_delaunay, get_patches_with_centroid


@pytest.mark.parametrize("method", (get_patches_delaunay, get_patches_with_centroid))
def test_patch_holes_001(method):
    """Test hole patcher for a Sphere with 2 holes."""
    # NOTE : should we assert volume?
    # Two holes to patch
    sphere: pv.PolyData = pv.Sphere()
    tube = pv.Tube(pointa=[0.0, 0.0, -1.0], pointb=[0.0, 0, 1], radius=0.2)
    tube = tube.triangulate()
    sphere = sphere.triangulate()

    sphere_with_hole = sphere.boolean_difference(tube)
    sphere_with_hole = sphere_with_hole.extract_cells(
        sphere_with_hole.cell_data["CellSource"] == 0
    ).extract_surface()

    patches = method(sphere_with_hole)
    assert len(patches) == 2
    assert pv.merge([sphere_with_hole] + patches).clean().is_manifold


@pytest.mark.parametrize("method", (get_patches_delaunay, get_patches_with_centroid))
def test_patch_holes_003(method):
    """Test patching a box."""
    box = pv.Box()
    box = box.triangulate()
    box = box.remove_cells([0, 1])

    patches = method(box)
    assert len(patches) == 1

    return


@pytest.mark.parametrize("method", (get_patches_delaunay, get_patches_with_centroid))
def test_patch_holes_002(method):
    """Test patch holes with single curved hole to patch."""

    sphere: pv.PolyData = pv.Sphere()
    tube = pv.Tube(pointa=[0.5, 1.0, 0.0], pointb=[0.5, -1, 0.0], radius=0.2)
    tube = tube.triangulate()
    sphere = sphere.triangulate()

    sphere_with_hole = sphere.boolean_difference(tube)
    sphere_with_hole = sphere_with_hole.extract_cells(
        sphere_with_hole.cell_data["CellSource"] == 0
    ).extract_surface()

    patches = method(sphere_with_hole)

    assert len(patches) == 1
    assert pv.merge([sphere_with_hole] + patches).clean().is_manifold


def test_boundary_edge_loops():
    from ansys.health.heart.utils.vtk_utils import get_boundary_edge_loops

    # prep data
    sphere: pv.PolyData = pv.Sphere()
    tube = pv.Tube(pointa=[0.0, 0.0, -1.0], pointb=[0.0, 0, 1], radius=0.2)
    tube = tube.triangulate()
    sphere = sphere.triangulate()

    sphere_with_hole = sphere.boolean_difference(tube)
    sphere_with_hole = sphere_with_hole.extract_cells(
        sphere_with_hole.cell_data["CellSource"] == 0
    ).extract_surface()

    edge_loops = get_boundary_edge_loops(
        sphere_with_hole, remove_open_edge_loops=True, return_types=False
    )

    assert len(edge_loops) == 2
