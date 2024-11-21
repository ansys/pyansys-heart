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

if os.getenv("GITHUB_ACTIONS"):
    is_gh_action = True
else:
    is_gh_action = False

import unittest.mock as mock

import numpy as np
import pytest
import pyvista as pv

import ansys.heart.core.models as models
from ansys.heart.core.objects import Part


@pytest.fixture
def _mock_input(mocker):
    """Mock biventricular model."""
    mesh = pv.examples.load_tetbeam()
    mesh.points = mesh.points * 1e1
    mesh.cell_data["part-id"] = 1
    fibers = np.repeat([1.0, 0.0, 0.0], mesh.n_cells, axis=0).reshape((3, mesh.n_cells)).T
    mesh.cell_data["fiber"] = fibers
    mock_biventricle: models.BiVentricle = models.BiVentricle()
    mock_biventricle.mesh = mesh
    mock_show = mocker.patch("pyvista.Plotter.show")
    return mock_biventricle, mock_show


def test_heart_model_plot_mesh(_mock_input):
    """Test plotting the mesh."""
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
