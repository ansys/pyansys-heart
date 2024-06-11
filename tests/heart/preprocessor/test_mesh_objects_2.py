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


from ansys.heart.preprocessor.mesh.objects_2 import Mesh
import numpy as np
import pytest  # noqa F401
import pyvista as pv  # noqa F401

from tests.heart.conftest import download_asset, get_assets_folder  # noqa F401


@pytest.mark.parametrize("dtype", [float, int])
def test_pyvista_clean_grid(dtype):
    import pyvista as pv

    if dtype == int:
        pytest.xfail(reason="clean() will fail if points are of dtype int")

    points = np.array([[0, 0, 0], [-1, 0, 0], [1, 0, 0]], dtype=dtype)

    # option 1
    line_cells = [2, 0, 1]
    cell_types = [pv.CellType.LINE]
    line_ugrid = pv.UnstructuredGrid(line_cells, cell_types, points)
    assert line_ugrid.clean().n_points == 2

    # option 2 (does not work)
    line = pv.lines_from_points(points[0:-1])
    line = line.cast_to_unstructured_grid()
    line.points = points
    assert line.clean().n_points == 2

    return


def test_mesh_add_001():
    """Test adding triangles and lines to a Mesh object."""
    # NOTE:
    # create a base mesh, and add the following:
    # - a volume mesh
    # - a surface mesh
    # - a beam mesh

    # generate some dummy data
    from pyvista import examples

    # prep data based on hexbeam
    tets = examples.load_tetbeam()
    tets.clear_cell_data()
    edges = tets.extract_feature_edges()
    edges.clear_cell_data()
    edges.clear_point_data()
    triangles = tets.extract_surface()
    triangles.clear_cell_data()
    triangles.clear_point_data()

    # init the base mesh.
    base = Mesh(tets)

    # fill attribute with some values to ensure that info is kept.
    base.conn["c1"] = [[1, 2]]

    merged = base._add_mesh(triangles)
    assert isinstance(merged, Mesh)
    assert base.n_cells == triangles.n_cells + tets.n_cells

    # assert adding edges.
    base._add_mesh(edges)
    assert base.n_cells == triangles.n_cells + tets.n_cells + edges.n_cells

    # assert that attribute was kept
    assert base.conn["c1"] == [[1, 2]]
    pass


def test_mesh_add_002():
    """Test adding another mesh and its data."""
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)

    # for each face in the tet
    tets = [4, 0, 1, 2, 3]
    cell_types = [pv.CellType.TETRA]

    grid = Mesh(tets, cell_types, points)

    # add dummy data
    grid.cell_data["cdata_tet1"] = 1
    grid.point_data["pdata_tet1"] = 2

    # add set of lines that is connected at a single point
    lines = pv.lines_from_points([points[0, :], [-1, 0, 0], [-2, 0, 0]])

    # create some dummy data
    lines["cdata_lines1"] = np.array([10, 11], dtype=np.int32)
    lines["cdata_lines2"] = np.array([10.5, 11.5], dtype=np.float64)
    lines["pdata_lines1"] = np.array([[100, 101], [102, 103], [104, 105]], dtype=np.float64)

    grid1 = grid._add_mesh(lines, keep_data=True)

    assert isinstance(grid1, Mesh)

    # check array names
    all_names = grid.array_names + lines.array_names
    for name in all_names:
        assert name in grid1.array_names, f"Array with name {name} does  not exist"

    # check data types
    assert isinstance(grid1.cell_data["cdata_lines1"][0], np.int32)
    assert isinstance(grid1.cell_data["cdata_lines2"][0], np.float64)

    # check dimension 2 of array
    assert grid1["pdata_lines1"].shape[1] == lines["pdata_lines1"].shape[1]

    pass


def test_surface_add_001():
    """Test adding a surface to a Mesh."""
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    tets = [4, 0, 1, 2, 3]
    cell_types = [pv.CellType.TETRA]

    mesh = Mesh(tets, cell_types, points)
    mesh.cell_data["volume-id"] = 1

    surface = pv.Triangle([points[0, :], [-1, 0, 0], [-1, -1, 0]])

    # test adding a single surface
    mesh.add_surface(surface, id=2)
    assert "surface-id" in mesh.cell_data.keys()

    assert np.all(mesh.celltypes == [pv.CellType.TRIANGLE, pv.CellType.TETRA])
    np.testing.assert_allclose(mesh.cell_data["volume-id"], [np.nan, 1])
    np.testing.assert_allclose(mesh.cell_data["surface-id"], [2, np.nan])

    # test adding multiple surfaces simultaneously
    mesh = Mesh(tets, cell_types, points)
    mesh.cell_data["volume-id"] = 1

    points1 = np.array([points[0, :], [-1, 0, 0], [-1, -1, 0]]) + 1.0
    surface1 = pv.Triangle(points1)
    surface.cell_data["surface-id"] = 10
    surface1.cell_data["surface-id"] = 11
    surface_to_add = surface + surface1

    mesh.add_surface(surface_to_add)
    np.testing.assert_allclose(mesh.cell_data["surface-id"], [11, 10, np.nan])


def test_lines_add_001():
    """Test adding a beam to a Mesh."""
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    tets = [4, 0, 1, 2, 3]
    cell_types = [pv.CellType.TETRA]

    mesh = Mesh(tets, cell_types, points)
    mesh.cell_data["volume-id"] = 1

    line = pv.Line(points[0, :], [-1, 0, 0])

    mesh.add_lines(line, id=2)
    assert "line-id" in mesh.cell_data.keys()
    assert np.all(mesh.celltypes == [pv.CellType.LINE, pv.CellType.TETRA])
    np.testing.assert_allclose(mesh.cell_data["volume-id"], [np.nan, 1])
    np.testing.assert_allclose(mesh.cell_data["line-id"], [2, np.nan])


def test_volume_add_001():
    """Test adding a volume (hex element) to an existing mesh."""
    from pyvista import examples

    points = np.array([[0, 0, 0], [-1, 0, 0], [0, -1, 0], [0, 0, 1]], dtype=float)
    tets = [4, 0, 1, 2, 3]
    cell_types = [pv.CellType.TETRA]
    mesh = Mesh(tets, cell_types, points)
    mesh.cell_data["volume-id"] = 1

    hex = examples.cells.Hexahedron()

    mesh.add_volume(hex, id=2)
    assert np.allclose(mesh.cell_data["volume-id"], [2, 1])


def test_mesh_remove_001():
    """Remove part of the mesh."""
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0]], dtype=float)
    # for each face in the tet
    tets = [4, 0, 1, 2, 3, 4, 0, 1, 2, 4]
    cell_types = [pv.CellType.TETRA] * 2

    grid = Mesh(tets, cell_types, points)
    grid.cell_data["test_data"] = 1

    grid2 = grid.remove_cells(0, inplace=False)
    assert isinstance(grid2, Mesh)
    assert grid2.n_cells == 1
    grid.remove_cells(0, inplace=True)
    assert isinstance(grid2, Mesh)
    assert grid.n_cells == 1


def test_mesh_clean_001():
    """Test cleaning method."""
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    tets = [4, 0, 1, 2, 3]
    cell_types = [pv.CellType.TETRA]
    grid = Mesh(tets, cell_types, points)
    # modify points
    grid.points = np.vstack((points, points))

    assert isinstance(grid, Mesh)
    grid1 = grid.clean()
    assert grid.n_points == points.shape[0] * 2
    assert grid1.n_points == points.shape[0]
    pass


def test_mesh_object_properties():
    """Test the properties of Mesh."""
    # NOTE: clean this up
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [-2, 0, 0]])

    # for each face in the tet
    tets = [4, 0, 1, 2, 3]  # 1 tet
    triangles = [3, 0, 1, 3, 3, 1, 2, 3, 3, 2, 0, 3, 3, 0, 2, 4, 3, 0, 4, 5]  # 5 triangles
    lines = [2, 0, 1, 2, 1, 2, 2, 2, 3, 2, 4, 5]  # 4 lines
    cell_types = [pv.CellType.TETRA] + [pv.CellType.TRIANGLE] * 5 + [pv.CellType.LINE] * 4
    cells = tets + triangles + lines
    mesh = Mesh(cells, cell_types, points)
    mesh.celltypes

    #
    data = np.ones(mesh.n_cells, dtype=int) * -1
    data[mesh.celltypes == pv.CellType.TETRA] = 1
    mesh.cell_data["volume-id"] = data
    assert mesh.volume_ids == [1]

    data = np.ones(mesh.n_cells, dtype=int) * -1
    data[mesh.celltypes == pv.CellType.LINE] = 1
    data[-2:] = 2
    mesh.cell_data["line-id"] = data
    assert np.all(mesh.line_ids == [1, 2])

    data = np.ones(mesh.n_cells, dtype=int) * -1
    data[mesh.celltypes == pv.CellType.TRIANGLE] = 1
    data[1] = 2
    data[2:4:] = 3
    mesh.cell_data["surface-id"] = data
    assert np.all(mesh.surface_ids == [1, 2, 3])

    pass
