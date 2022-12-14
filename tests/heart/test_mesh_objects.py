from ansys.heart.preprocessor.mesh.objects import Mesh, SurfaceMesh
import numpy as np
import pyvista as pv
import pytest


def test_mesh_object_assign_001():
    """Test the mesh object with a tetrahedron."""
    # test data:
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    tetrahedron = np.array([[0, 1, 2, 3]])

    mesh = Mesh()
    mesh.tetrahedrons = tetrahedron
    mesh.nodes = points

    mesh_pv = mesh._to_pyvista_object()

    assert np.array_equal(mesh_pv.points.shape, points.shape)
    assert np.all(np.isclose(mesh_pv.points, points))
    assert np.array_equal(mesh_pv.cells, np.append([4], tetrahedron.flatten()))


def test_mesh_object_assign_002():
    """Test the assigning nodes/faces to the surface mesh."""
    # test data:
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    triangles = np.array([[0, 1, 2]])

    surface = SurfaceMesh()
    surface.faces = triangles
    surface.nodes = points

    polydata_pv = surface._to_pyvista_object()

    assert np.array_equal(polydata_pv.points.shape, points.shape)
    assert np.all(np.isclose(polydata_pv.points, points))
    assert np.array_equal(polydata_pv.faces, np.append([3], triangles.flatten()))
