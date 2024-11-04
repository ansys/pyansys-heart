import os

if os.getenv("GITHUB_ACTIONS"):
    is_gh_action = True
else:
    is_gh_action = False

import numpy as np
import pytest
import unittest.mock as mock

import ansys.heart.core.models as models
import pyvista as pv
from ansys.heart.core.objects import Part
pv.OFF_SCREEN = True


@pytest.fixture
def _mock_biventricle():
    """Mock biventricular model."""
    mesh = pv.examples.load_tetbeam()
    mesh.points = mesh.points*1e1
    mesh.cell_data["part-id"] = 1    
    fibers = np.repeat([1.0, 0.0, 0.0], mesh.n_cells, axis = 0).reshape((3, mesh.n_cells)).T
    mesh.cell_data["fiber"] = fibers
    mock_biventricle: models.BiVentricle = models.BiVentricle(mock.Mock(models.ModelInfo))
    mock_biventricle.mesh = mesh
    return mock_biventricle

def test_heart_model_plot_mesh(_mock_biventricle):
    """Test plotting the mesh."""
    # mock .show or pv.OFF_SCREEN?
    _mock_biventricle.plot_mesh()
    

def test_heart_model_plot_part(_mock_biventricle):
    """Test plotting a part."""
    mock_part = mock.Mock(Part)
    mock_part.element_ids = np.array([0,1,2,3,4])
    _mock_biventricle.plot_part(mock_part)

def test_heart_model_plot_fibers(_mock_biventricle):
    """Test plotting of fibers."""
    assert isinstance(_mock_biventricle.plot_fibers(), pv.Plotter)
    _mock_biventricle.mesh.points = _mock_biventricle.mesh.points/1e1
    assert not _mock_biventricle.plot_fibers()