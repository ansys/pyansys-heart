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

import numpy as np
import pytest
import pyvista as pv

import ansys.health.heart.models as models
from ansys.health.heart.objects import Mesh, Part, SurfaceMesh
from tests.heart.writer.test_dynawriter import _get_mock_conduction_system


@pytest.fixture
def _mock_input():
    """Mock biventricular model."""
    mesh = pv.examples.load_tetbeam()
    mesh.points = mesh.points * 1e1
    mesh.cell_data["_volume-id"] = 1
    fibers = np.repeat([1.0, 0.0, 0.0], mesh.n_cells, axis=0).reshape((3, mesh.n_cells)).T
    mesh.cell_data["fiber"] = fibers

    mock_biventricle: models.BiVentricle = models.BiVentricle()
    mock_biventricle.mesh = Mesh(mesh)

    # add surfaces to parts.
    ii = 1
    for part in mock_biventricle.parts:
        for surface in part.surfaces:
            mock_biventricle.mesh.add_surface(pv.Disc(), id=ii, name=surface.name)
            surface.id = ii
            super(SurfaceMesh, surface).__init__(mock_biventricle.mesh.get_surface(surface.id))
            ii += 1

    mock_biventricle.mesh.add_surface(pv.Disc(), id=1, name="valve")

    # add mock purkinje data.
    mock_biventricle.conduction_system = _get_mock_conduction_system()

    with mock.patch("pyvista.Plotter.show") as mock_show:
        yield mock_biventricle, mock_show


def test_heart_model_plot_mesh(_mock_input):
    """Test plotting the mesh."""
    # with mock.patch("pyvista.Plotter.show") as mock_show:
    mock_biventricle, mock_show = _mock_input
    mock_biventricle.plot_mesh()
    mock_show.assert_called_once()


def test_heart_model_plot_part(_mock_input):
    """Test plotting a part."""
    mock_biventricle, mock_show = _mock_input
    mock_part = mock.Mock(Part)
    mock_part.element_ids = np.array([0, 1, 2, 3, 4])
    mock_biventricle.plot_part(mock_part)
    mock_show.assert_called_once()


def test_heart_model_plot_fibers(_mock_input):
    """Test plotting of fibers."""
    mock_biventricle, mock_show = _mock_input
    assert isinstance(mock_biventricle.plot_fibers(), pv.Plotter)
    mock_biventricle.mesh.points = mock_biventricle.mesh.points / 1e1
    assert not mock_biventricle.plot_fibers()
    mock_show.assert_called_once()


def test_heart_model_plot_surfaces(_mock_input):
    """Test plotting the surfaces of the model."""
    mock_biventricle, mock_show = _mock_input

    if not isinstance(mock_biventricle, models.BiVentricle):
        assert isinstance(mock_biventricle, models.BiVentricle)

    mock_biventricle.plot_surfaces()
    mock_show.assert_called_once()


def test_heart_model_plot_purkinje(_mock_input):
    """Test plotting purkinje."""
    mock_biventricle, mock_show = _mock_input

    if not isinstance(mock_biventricle, models.BiVentricle):
        assert isinstance(mock_biventricle, models.BiVentricle)

    mock_biventricle.plot_purkinje()

    mock_show.assert_called_once()
    mock_show.reset_mock()

    # remove beam network.
    mock_biventricle.conduction_system = None
    mock_biventricle.plot_purkinje()
    mock_show.assert_not_called()
