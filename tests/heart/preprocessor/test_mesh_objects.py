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
import tempfile

import numpy as np
import pytest
import pyvista as pv
from pyvista import examples

from ansys.heart.preprocessor.mesh.objects import Mesh, SurfaceMesh


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
    base.boundaries = [[1, 2]]

    merged = base._add_mesh(triangles)
    assert isinstance(merged, Mesh)
    assert base.n_cells == triangles.n_cells + tets.n_cells

    # assert adding edges.
    base._add_mesh(edges)
    assert base.n_cells == triangles.n_cells + tets.n_cells + edges.n_cells

    # assert that attribute was kept
    assert base.boundaries == [[1, 2]]
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
    mesh.cell_data["_volume-id"] = 1

    surface = pv.Triangle([points[0, :], [-1, 0, 0], [-1, -1, 0]])

    assert mesh.add_surface(surface) == None
    assert mesh.add_surface(surface, float(1)) == None

    # test adding a single surface
    mesh.add_surface(surface, id=2)
    assert "_surface-id" in mesh.cell_data.keys()

    assert np.all(mesh.celltypes == [pv.CellType.TRIANGLE, pv.CellType.TETRA])
    np.testing.assert_allclose(mesh.cell_data["_volume-id"], [np.nan, 1])
    np.testing.assert_allclose(mesh.cell_data["_surface-id"], [2, np.nan])

    # test adding multiple surfaces simultaneously with celldata _surface-id
    mesh = Mesh(tets, cell_types, points)
    mesh.cell_data["_volume-id"] = 1

    points1 = np.array([points[0, :], [-1, 0, 0], [-1, -1, 0]]) + 1.0
    surface1 = pv.Triangle(points1)
    surface.cell_data["_surface-id"] = 10
    surface1.cell_data["_surface-id"] = 11
    surface_to_add = surface + surface1

    mesh.add_surface(surface_to_add)
    np.testing.assert_allclose(mesh.cell_data["_surface-id"], [11, 10, np.nan])


def test_lines_add_001():
    """Test adding a beam to a Mesh."""
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    tets = [4, 0, 1, 2, 3]
    cell_types = [pv.CellType.TETRA]

    mesh = Mesh(tets, cell_types, points)
    mesh.cell_data["_volume-id"] = 1

    line = pv.Line(points[0, :], [-1, 0, 0])

    assert mesh.add_lines(line) == None
    assert mesh.add_lines(line, float(1)) == None

    mesh.add_lines(line, id=2)
    assert "_line-id" in mesh.cell_data.keys()
    assert np.all(mesh.celltypes == [pv.CellType.LINE, pv.CellType.TETRA])
    np.testing.assert_allclose(mesh.cell_data["_volume-id"], [np.nan, 1])
    np.testing.assert_allclose(mesh.cell_data["_line-id"], [2, np.nan])


def test_volume_add_001():
    """Test adding a volume (hex element) to an existing mesh."""
    from pyvista import examples

    points = np.array([[0, 0, 0], [-1, 0, 0], [0, -1, 0], [0, 0, 1]], dtype=float)
    tets = [4, 0, 1, 2, 3]
    cell_types = [pv.CellType.TETRA]
    mesh = Mesh(tets, cell_types, points)
    mesh.cell_data["_volume-id"] = 1
    hex = examples.cells.Hexahedron()

    assert mesh.add_volume(hex) == None
    assert mesh.add_volume(hex, float(1)) == None

    mesh.add_volume(hex, id=2)
    assert np.allclose(mesh.cell_data["_volume-id"], [2, 1])


def test_get_submesh_001():
    """Test getting a submesh of Mesh."""
    # generate some dummy data
    from pyvista import examples

    # prep data based on hexbeam
    tets = examples.load_tetbeam()
    tets.clear_cell_data()
    tets.cell_data["_volume-id"] = 1
    lines = tets.extract_feature_edges()
    lines.clear_cell_data()
    lines.clear_point_data()
    triangles = tets.extract_surface()
    triangles.clear_cell_data()
    triangles.clear_point_data()
    triangles = triangles.remove_cells(np.arange(0, 20))

    mesh = Mesh()
    mesh.add_volume(tets, 1)
    mesh.add_surface(triangles, 10)
    mesh.add_lines(lines, 100)

    # test get surfaces
    triangles1 = mesh.get_surface(10)
    assert triangles.n_cells == triangles1.n_cells
    assert np.allclose(triangles.points, triangles1.points)

    # test get lines
    lines1 = mesh.get_lines(100)
    assert lines.n_cells == lines1.n_cells
    assert lines.n_points == lines1.n_points

    mesh.add_surface(triangles, 11)
    triangles1 = mesh.get_surface([10, 11])
    assert triangles1.n_cells == triangles.n_cells * 2

    # test get volume
    mesh = Mesh()
    tets.cell_data["_volume-id"] = 1
    tets.cell_data["_volume-id"][0:20] = 2
    mesh.add_volume(tets)
    assert np.allclose(mesh.volume_ids, [1, 2])

    volume1 = mesh.get_volume(2)
    assert volume1.n_cells == 20


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


def test_mesh_remove_002():
    """Remove surface based on id."""
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    tets = [4, 0, 1, 2, 3]
    cell_types = [pv.CellType.TETRA]

    mesh = Mesh(tets, cell_types, points)
    mesh.cell_data["_volume-id"] = 1

    surface = pv.Triangle([points[0, :], [-1, 0, 0], [-1, -1, 0]])

    # add a single surface
    mesh.add_surface(surface, id=2)
    # remove this surface
    mesh.remove_surface(sid=2)
    assert mesh.n_cells == 1


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
    mesh.cell_data["_volume-id"] = data
    assert mesh.volume_ids == [1]

    data = np.ones(mesh.n_cells, dtype=int) * -1
    data[mesh.celltypes == pv.CellType.LINE] = 1
    data[-2:] = 2
    mesh.cell_data["_line-id"] = data
    assert np.all(mesh.line_ids == [1, 2])

    data = np.ones(mesh.n_cells, dtype=int) * -1
    data[mesh.celltypes == pv.CellType.TRIANGLE] = 1
    data[1] = 2
    data[2:4:] = 3
    mesh.cell_data["_surface-id"] = data
    assert np.all(mesh.surface_ids == [1, 2, 3])

    pass


def test_add_nodes_mesh():
    """Test adding points/nodes to the Mesh object."""
    # test data:
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    tetrahedron = np.array([[0, 1, 2, 3]])

    mesh = Mesh()
    mesh.tetrahedrons = tetrahedron
    mesh.nodes = points
    mesh.point_data["data-scalar"] = np.ones(mesh.n_points, dtype=float)
    mesh.point_data["data-vector"] = np.ones((mesh.n_points, 3), dtype=float)

    # test adding nodes
    mesh.nodes = np.vstack([points, [0, 0.5, 0.5]])
    assert mesh.point_data["data-scalar"].shape[0] == mesh.nodes.shape[0]
    assert mesh.point_data["data-vector"].shape[0] == mesh.nodes.shape[0]

    # test assigning same number of nodes
    mesh.nodes = mesh.nodes * 1e-3

    assert mesh.point_data["data-scalar"].shape[0] == mesh.nodes.shape[0]
    assert mesh.point_data["data-vector"].shape[0] == mesh.nodes.shape[0]


def test_surface_mesh_init():
    """Test whether SurfaceMesh __init__'s behaves similar to pv.PolyData init"""

    sphere = pv.Sphere()

    ## init without var_inp
    surface = SurfaceMesh(name="test_name")
    assert surface.name == "test_name"

    # init with number of points
    sphere1 = SurfaceMesh(sphere.points, faces=sphere.faces, name="test_name")

    assert np.allclose(sphere1.points, sphere.points)
    assert np.all(sphere1.faces == sphere.faces)
    assert sphere1.name == "test_name"

    # Init with polydata
    sphere1 = SurfaceMesh(sphere)
    assert sphere1.n_cells == sphere.n_cells
    assert sphere1.n_points == sphere.n_points

    # Init with file path
    with tempfile.TemporaryDirectory(suffix=".pyansys-heart") as tempdir:
        temp_path = os.path.join(tempdir, "sphere.vtp")
        sphere.save(temp_path)

        sphere1 = SurfaceMesh(temp_path)

    assert sphere1.n_cells == sphere.n_cells
    assert sphere1.n_points == sphere.n_points


@pytest.mark.parametrize("celltype", ["quads", "triangles"])
def test_mesh_boundaries_property(celltype):
    """Test the boundaries property of Mesh."""
    if celltype == "triangles":
        polydata = pv.Sphere()
    elif celltype == "quads":
        polydata = pv.Box(level=15)

    polydata.cell_data["_surface-id"] = 1
    polydata.cell_data["_surface-id"][0:800] = 2
    mesh = Mesh()
    assert mesh._boundaries == []
    mesh.add_surface(polydata)
    assert len(mesh._boundaries) == 2


@pytest.mark.parametrize("celltype", ["hex", "tets"])
def test_mesh_volumes_property(celltype):
    """Test the boundaries property of Mesh."""
    if celltype == "hex":
        ugrid = examples.load_hexbeam()
    elif celltype == "tets":
        ugrid = examples.load_tetbeam()

    ugrid.cell_data["_volume-id"] = 1
    ugrid.cell_data["_volume-id"][0:20] = 2
    mesh = Mesh()
    assert mesh._volumes == []
    mesh.add_volume(ugrid)
    assert len(mesh._volumes) == 2
