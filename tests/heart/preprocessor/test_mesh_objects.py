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

from ansys.heart.preprocessor.mesh.objects import Mesh, SurfaceMesh
import numpy as np
import pytest  # noqa F401
import pyvista as pv  # noqa F401

from tests.heart.conftest import download_asset, get_assets_folder  # noqa F401

skip_test = os.name != "nt"


def test_mesh_object_assign_001():
    """Test the mesh object with a tetrahedron."""
    # test data:
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    tetrahedron = np.array([[0, 1, 2, 3]])

    mesh = Mesh()
    mesh.tetrahedrons = tetrahedron
    mesh.nodes = points

    assert np.array_equal(mesh.points.shape, points.shape)
    assert np.all(np.isclose(mesh.points, points))
    assert np.array_equal(mesh.cells, np.append([4], tetrahedron.flatten()))


# def test_surfacemesh_object_assign_001():
#     """Test the assigning nodes/faces to the surface mesh."""
#     # test data:
#     points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
#     tris = np.array([[0, 1, 2]])

#     surface = SurfaceMesh()
#     surface.triangles = tris
#     surface.nodes = points

#     polydata_pv = surface._to_pyvista_object()

#     assert np.array_equal(polydata_pv.points.shape, points.shape)
#     assert np.all(np.isclose(polydata_pv.points, points))
#     assert np.array_equal(polydata_pv.faces, np.append([3], tris.flatten()))


def test_surfacemesh_object_assign_001():
    """Test the assigning of cell and point data to SurfaceMesh."""
    # test data:
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    tris = np.array([[0, 1, 2], [0, 1, 3]])
    cell_data = np.array([1.0, 2.0])
    point_data = np.array([1.0, 2.0, 3.0, 4.0])

    surface = SurfaceMesh(name="test")
    surface.triangles = tris
    surface.nodes = points
    surface.cell_data["cells"] = cell_data
    surface.point_data["points"] = point_data

    assert np.array_equal(surface.nodes.shape, points.shape)
    assert np.all(np.isclose(surface.points, points))
    assert np.array_equal(surface.faces, np.hstack([[[3], [3]], tris]).flatten())
    assert np.array_equal(surface.triangles, tris)
    assert np.all(surface.point_data["points"] == point_data)
    assert np.all(surface.cell_data["cells"] == cell_data)


def test_mesh_object_read_001():
    """Test if strocchi's case is properly read."""
    case_file = download_asset(database="Strocchi2020", casenumber=1)
    mesh = Mesh()
    mesh.read_mesh_file(case_file)

    assert mesh.n_cells == 2349414
    assert mesh.n_points == 481066
    assert "tags" in mesh.array_names
    return


# @pytest.mark.skipif(skip_test)
# @pytest.mark.xfail(reason="Only supporting Strocchi for tests.")
def test_mesh_object_read_002():
    """Test if rodero's case is properly read."""
    case_file = download_asset(database="Rodero2021", casenumber=1)
    mesh = Mesh()
    mesh.read_mesh_file_rodero2021(case_file)
    assert mesh.n_cells == 2477866
    assert mesh.n_points == 493151
    assert "tags" in mesh.array_names

    return


def test_add_nodes_mesh():
    """Test adding points/nodes to the Mesh object."""
    # test data:
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    tetrahedron = np.array([[0, 1, 2, 3]])

    mesh = Mesh()
    mesh.tetrahedrons = tetrahedron
    mesh.nodes = points
    mesh.point_data["data-scalar"] = np.ones(mesh.n_points, dtype=float)
    mesh.point_data["data-vector"] = np.ones((mesh.n_points, 3), dtype=float)

    # test adding nodes
    mesh.nodes = np.vstack([points, [0, 0.5, 0.5]])
    assert mesh.point_data["data-scalar"].shape[0] == mesh.nodes.shape[0]
    assert mesh.point_data["data-vector"].shape[0] == mesh.nodes.shape[0]

    # test assigning same number of nodes
    mesh.nodes = mesh.nodes * 1e-3

    assert mesh.point_data["data-scalar"].shape[0] == mesh.nodes.shape[0]
    assert mesh.point_data["data-vector"].shape[0] == mesh.nodes.shape[0]


def test_add_nodes_surface_mesh():
    """Test adding points/nodes to the SurfaceMesh object."""
    # test data:
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    triangles = np.array([[0, 1, 2]])

    mesh = SurfaceMesh("triangle", triangles=triangles, nodes=points)
    mesh.point_data["data-scalar"] = np.ones(mesh.n_points, dtype=float)
    mesh.point_data["data-vector"] = np.ones((mesh.n_points, 3), dtype=float)

    # test adding nodes
    mesh.nodes = np.vstack([points, [0, 0.5, 0.5]])
    assert mesh.point_data["data-scalar"].shape[0] == mesh.nodes.shape[0]
    assert mesh.point_data["data-vector"].shape[0] == mesh.nodes.shape[0]

    # test assigning same number of nodes
    mesh.nodes = mesh.nodes * 1e-3

    assert mesh.point_data["data-scalar"].shape[0] == mesh.nodes.shape[0]
    assert mesh.point_data["data-vector"].shape[0] == mesh.nodes.shape[0]
