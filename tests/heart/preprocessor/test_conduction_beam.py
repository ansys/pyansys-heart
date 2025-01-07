"""Contains tests for conduction system."""

import numpy as np
from pyvista.examples import load_tetbeam

from ansys.heart.preprocessor.conduction_beam import ConductionSystem, _create_line, _refine_line


def test_find_path():
    mesh = load_tetbeam()
    mesh.point_data["_global-point-ids"] = np.arange(0, mesh.n_points)
    path = ConductionSystem.find_path(mesh, np.array([0, 0, 0]), np.array([1, 0, 5]), True)

    # TODO: add reasonable assertion
    assert path[0].shape[0] > 0
    assert len(path[1]) > 0


def test_create_line():
    points = _create_line(np.array([0, 0, 0]), np.array([0, 0, 1]), 0.5)
    assert np.allclose(points, np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.5], [0.0, 0.0, 1.0]]))


def test_refine_line():
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    refined_points = _refine_line(points, 0.5)
    assert np.allclose(
        refined_points, np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0]])
    )
