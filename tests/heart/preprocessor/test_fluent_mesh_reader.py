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

"""Some tests to test function of the Fluent mesh reader on a unit-cube example"""

import os

import numpy as np
import pytest
import pyvista as pv

from ansys.heart.core.utils.fluent_reader import _FluentMesh
from tests.heart.conftest import get_assets_folder

FLUENT_BOX = os.path.join(get_assets_folder(), "simple_fluent_meshes", "box.msh.h5")
FLUENT_BOX1 = os.path.join(get_assets_folder(), "simple_fluent_meshes", "box_no_cells.msh.h5")


@pytest.fixture(autouse=True, scope="module")
def _test_mesh():
    mesh = _FluentMesh()
    mesh._open_file(FLUENT_BOX)
    yield mesh
    mesh._close_file()


@pytest.mark.parametrize(
    "input,expected_error",
    (
        [None, ValueError],
        ["file-without-extension", FileNotFoundError],
        ["non-existing-file.msh.h5", FileNotFoundError],
    ),
)
def test_open_file(input, expected_error):
    """Tests opening a file."""
    mesh = _FluentMesh()
    if expected_error:
        with pytest.raises(expected_error):
            mesh._open_file(input)
        return


def test_load_mesh():
    """Tests loading a mesh without cell zones."""
    mesh = _FluentMesh(FLUENT_BOX1)
    mesh.load_mesh(reconstruct_tetrahedrons=False)

    assert sorted([fz.name for fz in mesh.face_zones]) == [
        "box1-xmax",
        "box1-xmin",
        "box1-ymax",
        "box1-ymin",
        "box1-zmax",
        "box1-zmin",
    ]
    return


def test_to_vtk():
    """Tests conversion to vtk object with the various options available."""
    # mesh without any cell zones
    mesh = _FluentMesh(FLUENT_BOX1)
    mesh.load_mesh(reconstruct_tetrahedrons=False)
    vtk_object = mesh._to_vtk(add_cells=False, add_faces=True)

    # just triangles in vtk object
    assert list(vtk_object.cells_dict.keys()) == [int(pv.CellType.TRIANGLE)]
    assert vtk_object.n_cells == 6 * 2  # 2 faces per face zone, and 6 face zones

    # mesh with cell zones
    mesh = _FluentMesh(FLUENT_BOX)
    mesh.load_mesh(reconstruct_tetrahedrons=True)
    vtk_object = mesh._to_vtk(add_cells=True, add_faces=False)

    assert list(vtk_object.cells_dict.keys()) == [int(pv.CellType.TETRA)]
    assert vtk_object.n_cells == 12  # 12 tetrahedrons in cell zone

    # when extracting both cell zones and face zones
    vtk_object = mesh._to_vtk(add_cells=True, add_faces=True)

    # both triangles and tetra present
    assert list(vtk_object.cells_dict.keys()) == [int(pv.CellType.TRIANGLE), int(pv.CellType.TETRA)]
    assert vtk_object.n_cells == 42  # interior cells are also generated

    return


def test_read_nodes(_test_mesh):
    """Tests reading of the nodes of simple box"""
    mesh = _test_mesh
    assert isinstance(mesh, _FluentMesh)
    mesh._read_nodes()
    mesh._remove_duplicate_nodes()
    assert mesh.nodes.shape == (9, 3)
    expected_nodes = np.array(
        [
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
            [0.4965233, 0.5034767, 0.4965233],
        ]
    )
    assert np.allclose(np.sort(mesh.nodes, axis=0), np.sort(expected_nodes, axis=0), atol=1e-8)


def test_read_face_zones(_test_mesh):
    """Tests reading face zones. Checks face zone names and defined faces"""
    mesh = _test_mesh
    assert isinstance(mesh, _FluentMesh)
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


def test_clean_mesh():
    """Test cleaning the mesh from unreferenced nodes/points."""
    mesh = _FluentMesh()
    mesh.load_mesh(FLUENT_BOX)

    # remove interior faces
    mesh.face_zones = [fz for fz in mesh.face_zones if "interior" not in fz.name]
    # remove one face zone from the box:
    mesh.face_zones = mesh.face_zones[0:-1]

    mesh.clean()
    # even though a face zone is removed, the expected shape is
    # shape is (9,3) since it is used in a cell.
    assert mesh.nodes.shape == (9, 3)

    # remove cells and cell zones
    mesh.cell_zones = []
    mesh.cells = []

    mesh.clean()

    # since cells are also removed, we now expect 8 nodes and 10 faces
    assert np.vstack([fz.faces for fz in mesh.face_zones]).shape[0] == 10
    assert mesh.nodes.shape == (8, 3)

    return


def test_read_tetrahedrons(_test_mesh):
    """Tests reading of tetrahedrons on simple box with single cell zone"""
    mesh = _FluentMesh()
    mesh.load_mesh(FLUENT_BOX)
    expected_cells = mesh._unique_map[
        np.array(
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
    ]
    assert np.all(mesh.cell_zones[0].cells == expected_cells)
    # single cell zone: so all cells in mesh should be in the first cell zone.
    assert np.all(mesh.cells == mesh.cell_zones[0].cells)
