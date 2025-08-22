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

import unittest.mock as mock

import pytest
import pyvista as pv
from pyvista.examples import examples

from ansys.health.heart.pre.database_utils import (
    _get_interface_surfaces,
    _get_original_labels,
    _read_input_mesh,
    _Rodero2021_labels,
    _smooth_boundary_edges,
    _Strocchi2020_labels,
)


def test_read_input_mesh_strocchi():
    mesh = examples.load_tetbeam()

    with mock.patch("pyvista.read", return_value=mesh) as mock_read:
        mesh1 = _read_input_mesh(mesh, database="Strocchi2020")
        mock_read.assert_called_once()

        assert mesh1 == mesh

    multi_block = pv.MultiBlock([mesh])
    with mock.patch("pyvista.read", return_value=multi_block) as mock_read:
        mesh1 = _read_input_mesh(mesh, database="Strocchi2020")
        mock_read.assert_called_once()

        assert mesh1 == multi_block[0]

    with mock.patch("pyvista.read", return_value=[]) as mock_read:
        with pytest.raises(Exception):
            mesh1 = _read_input_mesh(mesh, database="Strocchi2020")

        mock_read.assert_called_once()


def test_read_input_mesh_rodero():
    mesh = examples.load_tetbeam()

    with mock.patch("pyvista.read", return_value=mesh) as mock_read:
        with pytest.raises(KeyError):
            mesh1 = _read_input_mesh(mesh, database="Rodero2021")

        mock_read.reset_mock()
        mesh.cell_data["ID"] = 1
        mesh1 = _read_input_mesh(mesh, database="Rodero2021")
        mock_read.assert_called_once()

        assert mesh1 == mesh


def test_get_interface_surface():
    # create test data (two connected tetbeams)
    mesh1 = examples.load_tetbeam()
    mesh2: pv.UnstructuredGrid = examples.load_tetbeam().translate((1, 0, 0))
    mesh1.cell_data["tags"] = 1
    mesh2.cell_data["tags"] = 2
    mesh = pv.merge([mesh1, mesh2])
    label_to_tag = {"mesh1": 1, "mesh2": 2}
    tag_to_label = {v: k for k, v in label_to_tag.items()}

    expected_labels = label_to_tag
    expected_labels["interface_mesh1_mesh2"] = 4

    interfaces, labels = _get_interface_surfaces(mesh, label_to_tag, tag_to_label)

    assert labels == expected_labels
    assert len(interfaces) == 1
    # interface cells all expected to have x coordinate 1.
    assert all(interfaces[0].clean().points[:, 0] == 1)


@pytest.mark.parametrize(
    "database,expected",
    (["Strocchi2020", _Strocchi2020_labels], ["Rodero2021", _Rodero2021_labels], ["unknown", None]),
)
def test_get_original_labels(database, expected):
    assert _get_original_labels(database) == expected


def test_smooth_boundary_edges():
    # test data
    poly = pv.Polygon(radius=10, n_sides=10)
    poly = poly.triangulate()
    poly.cell_data["surface-id"] = 1
    poly = poly.rotate_x(90)
    id_to_label = {1: "polygon"}
    # will return nan points if all points are in X-Y plane
    smoothed_plane, edges = _smooth_boundary_edges(poly, id_to_label, sub_label_to_smooth="polygon")
    assert len(edges) == 1
