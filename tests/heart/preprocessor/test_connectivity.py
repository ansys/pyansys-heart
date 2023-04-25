"""Tests for the connectivity module."""

from ansys.heart.preprocessor.mesh.connectivity import get_surfaces_from_tetmesh
import numpy as np

from tests.heart.preprocessor.common import _create_simple_unstructured_grid


def test_get_surfaces_from_tetrahedrons():
    """Test getting the boundary and interface surfaces."""

    test_grid = _create_simple_unstructured_grid()

    # expected boundaries:
    reference_interface_faces = [[0, 1, 2]]
    reference_boundary_faces = [[0, 1, 3], [0, 2, 3], [1, 2, 3], [0, 1, 4], [0, 2, 4], [1, 2, 4]]

    try:
        tets = np.reshape(test_grid.cells, (test_grid.n_cells, 5))[:, 1:]
    except:
        raise ValueError("Failed to extract tetrahedrons.")

    part_ids = test_grid.cell_data["tags"]

    boundary_faces, part_ids_boundary, interface_faces, pairs = get_surfaces_from_tetmesh(
        tets, return_boundaries=True, return_interfaces=True, part_ids=part_ids
    )

    assert np.all(
        np.equal(reference_boundary_faces, boundary_faces)
    ), "Boundary faces not the same."
    assert np.all(
        np.equal(reference_interface_faces, interface_faces)
    ), "Interface faces not the same."
    assert np.all(
        np.equal([1, 1, 1, 2, 2, 2], part_ids_boundary)
    ), "Boudary part ids not as expected."
    assert np.all(pairs == [test_grid.cell_data["tags"]]), "Pairs not the same."

    pass


"""Inexhaustive list of tests. Can consider adding:
1. get_faces_tetra
2. tetra_to_faces
3. face_tetra_connectivity
4. get_face_type
5. get_edges_from_triangles
...
"""
