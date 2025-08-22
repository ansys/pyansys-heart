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

"""Test connectivity methods."""

import numpy as np
import pytest
import pyvista as pv

from ansys.health.heart.utils.connectivity import (
    edge_connectivity,
    face_tetra_connectivity,
    get_edges_from_triangles,
    get_face_type,
    get_faces_tetra,
    remove_triangle_layers_from_trimesh,
)


def _get_test_data():
    tetrahedrons = np.array([[0, 1, 2, 3], [0, 1, 2, 4]])
    return tetrahedrons


def test_get_faces_tetra():
    """Test getting the faces of a tetrahedron."""
    tetrahedron = _get_test_data()[0:1]
    expected_faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    for ii, face in enumerate(expected_faces):
        assert np.all(get_faces_tetra(tetrahedron)[:, :, ii].flatten() == face)


def test_face_tetra_connectivity():
    """Test getting the faces of a tetrahedron."""
    # 1 shared face, 8 total faces. Hence, 7 unique faces in total.
    tetrahedrons = _get_test_data()
    faces, c0, c1 = face_tetra_connectivity(tetrahedrons)
    assert faces.shape == (7, 3)
    # face 1 connected to both cells
    # remaining faces connected to either cell 0 or cell 1.
    assert np.all(c0 == np.array([0, 0, 0, 0, 1, 1, 1]))
    assert np.all(c1 == np.array([1, 0, 0, 0, 1, 1, 1]))


def test_get_face_type():
    """Test getting face types."""
    tetrahedrons = _get_test_data()
    faces, c0, c1 = face_tetra_connectivity(tetrahedrons)
    face_types = get_face_type(faces, np.array([c0, c1]).T)
    # 1: interior face, 2: boundary face
    assert np.all(face_types == [1, 2, 2, 2, 2, 2, 2])


def test_get_edges_from_triangles():
    """Test getting edges from triangles."""
    triangles = np.array([[0, 1, 2]])
    edges = get_edges_from_triangles(triangles)
    assert np.all(edges == np.array([[0, 1], [1, 2], [0, 2]]))


def test_remove_triangles_layers_from_trimesh():
    """Test removing boundary faces from a triangular mesh."""
    # square with 3x3 triangles: only 2 central triangles will remain
    # after removing a triangular layer.
    plane = pv.Plane(i_resolution=3, j_resolution=3).triangulate()
    triangles = plane.cast_to_unstructured_grid().cells_dict[pv.CellType.TRIANGLE]
    triangles_reduced = remove_triangle_layers_from_trimesh(triangles)
    assert triangles_reduced.shape == (2, 3)


@pytest.mark.parametrize(
    "edges,expected_num_groups,expected_type",
    [
        (np.array([[0, 1], [1, 2], [0, 2]]), 1, ["closed"]),
        (np.array([[0, 1], [1, 2], [0, 2], [4, 5], [5, 6]]), 2, ["closed", "open"]),
    ],
    ids=["one-closed-loop", "one-closed-one-open"],
)
def test_edge_connectivity(edges, expected_num_groups, expected_type):
    """Test edge connectivity."""

    edge_groups, edge_type = edge_connectivity(edges, return_type=True, sort_closed=True)

    assert len(edge_groups) == expected_num_groups
    assert edge_type == expected_type
