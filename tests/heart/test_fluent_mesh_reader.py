""" Some tests to test function of the Fluent mesh reader on a unit-cube example """
import os

from ansys.heart.preprocessor.mesh.fluenthdf5 import FluentMesh
import numpy as np
import pytest

from .conftest import get_assets_folder

FLUENT_BOX = os.path.join(get_assets_folder(), "simple_fluent_meshes", "box.msh.h5")


@pytest.fixture(autouse=True, scope="module")
def _test_mesh():
    mesh = FluentMesh()
    mesh._open_file(FLUENT_BOX)
    yield mesh
    mesh._close_file()


def test_read_nodes(_test_mesh):
    """Tests reading of the nodes of simple box"""
    mesh = _test_mesh
    assert isinstance(mesh, FluentMesh)
    mesh._read_nodes()
    expected_nodes = np.array(
        [
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.49652330, 0.50347669, 0.49652330],
        ]
    )
    assert np.allclose(mesh.nodes, expected_nodes, atol=1e-8)


def test_read_face_zones(_test_mesh):
    """Tests reading face zones. Checks face zone names and defined faces"""
    mesh = _test_mesh
    assert isinstance(mesh, FluentMesh)
    mesh._read_nodes()
    mesh._read_face_zone_info()
    mesh._read_cell_zone_info()
    mesh._read_all_faces_of_face_zones()
    expected_names = [
        "interior-43178",
        "box1-xmin",
        "box1-xmax",
        "box1-ymin",
        "box1-ymax",
        "box1-zmin",
        "box1-zmax",
    ]
    all_in_expected = True
    for face_zone in mesh.face_zones:
        if face_zone.name not in expected_names:
            all_in_expected = False
    assert all_in_expected, "One or more face zones not in list of expected face zones"

    expected_face_zones = {
        "interior-43178": [
            [33, 29, 31],
            [33, 31, 28],
            [33, 28, 29],
            [33, 29, 30],
            [33, 30, 31],
            [33, 32, 28],
            [33, 31, 32],
            [33, 30, 32],
            [33, 28, 26],
            [33, 26, 29],
            [33, 26, 30],
            [33, 32, 27],
            [33, 27, 28],
            [33, 30, 27],
            [33, 27, 26],
            [33, 30, 25],
            [33, 25, 27],
            [33, 26, 25],
        ],
        "box1-ymin": [[29, 26, 30], [30, 26, 25]],
        "box1-zmin": [[29, 30, 31], [32, 31, 30]],
        "box1-xmin": [[31, 28, 29], [28, 26, 29]],
        "box1-xmax": [[32, 30, 27], [25, 27, 30]],
        "box1-ymax": [[31, 32, 28], [28, 32, 27]],
        "box1-zmax": [[25, 26, 27], [28, 27, 26]],
    }
    for face_zone in mesh.face_zones:
        assert np.all(expected_face_zones[face_zone.name] == face_zone.faces)


def test_read_tetrahedrons(_test_mesh):
    """Tests reading of tetrahedrons on simple box with single cell zone"""
    mesh = FluentMesh()
    mesh.load_mesh(FLUENT_BOX)
    expected_cells = np.array(
        [
            [30, 27, 28, 32],
            [32, 28, 30, 29],
            [32, 30, 27, 31],
            [32, 29, 30, 31],
            [32, 27, 28, 25],
            [32, 28, 29, 25],
            [32, 31, 27, 26],
            [32, 29, 31, 26],
            [32, 27, 25, 26],
            [32, 29, 26, 24],
            [32, 25, 29, 24],
            [32, 26, 25, 24],
        ]
    )
    assert np.all(mesh.cell_zones[0].cells == expected_cells)
    # single cell zone: so all cells in mesh should be in the first cell zone.
    assert np.all(mesh.cells == mesh.cell_zones[0].cells)
