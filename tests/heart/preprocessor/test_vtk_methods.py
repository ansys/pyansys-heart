"""Collection of tests for vtk methods."""
from ansys.heart.preprocessor.mesh.vtkmethods import get_cells_with_scalar_value
import numpy as np
import pyvista as pv
from pyvista import examples
import pytest


def test_get_cells_with_scalar_value():
    grid = examples.load_hexbeam()
    scalars = np.arange(grid.n_cells)
    scalars[0:19] = 1
    scalars[19:] = 2
    grid.cell_data.set_scalars(scalars, "tags")

    grid_1 = get_cells_with_scalar_value(grid, 1, "tags")
    grid_2 = get_cells_with_scalar_value(grid, 2, "tags")

    assert grid_1.n_cells == 19
    assert grid_2.n_cells == 21

    grid_3 = get_cells_with_scalar_value(grid, [1], "tags")
    assert grid_3.n_cells == 19


def test_set_get_cell_point_data():
    """Test setting and getting cell"""
    sphere = pv.Sphere()
    sphere.cell_data.set_scalars(name="surface-id-1", scalars=np.ones(sphere.n_cells, dtype=int))
    assert "surface-id-1" in list(sphere.cell_data.keys())
    assert np.all(sphere.cell_data["surface-id-1"] == np.ones(sphere.n_cells, dtype=int))

    sphere.cell_data["surface-id-2"] = np.ones(sphere.n_cells, dtype=int)
    assert "surface-id-2" in list(sphere.cell_data.keys())
    assert np.all(sphere.cell_data["surface-id-2"] == np.ones(sphere.n_cells, dtype=int))


@pytest.mark.parametrize("datatype", [np.int64, np.int32, int])
def test_add_polydata(datatype):
    """Test adding polydata's with scalar data"""
    sphere1 = pv.Sphere(radius=5)
    sphere1.cell_data.set_scalars(scalars=np.ones(sphere1.n_cells, dtype=int), name="region")

    sphere2 = pv.Sphere(radius=10)
    sphere2.cell_data.set_scalars(scalars=np.zeros(sphere1.n_cells, dtype=datatype), name="region")

    sphere: pv.PolyData = sphere1 + sphere2

    assert "region" in list(sphere.cell_data.keys())
