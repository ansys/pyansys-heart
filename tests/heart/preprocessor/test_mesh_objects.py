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
@pytest.mark.xfail(reason="Only supporting Strocchi for tests.")
def test_mesh_object_read_002():
    """Test if rodero's case is properly read."""
    case_file = download_asset(database="Cristobal2021", casenumber=1)
    mesh = Mesh()
    mesh.read_mesh_file_cristobal2021(case_file)
    assert mesh.n_cells == 2477866
    assert mesh.n_points == 493151
    assert "tags" in mesh.array_names

    return
