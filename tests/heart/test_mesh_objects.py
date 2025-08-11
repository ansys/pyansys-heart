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


import os
import tempfile
from typing import Literal

import numpy as np
import pytest
import pyvista as pv
from pyvista import examples

from ansys.health.heart.objects import CapType, Mesh, SurfaceMesh

SURFACE_TYPES = [pv.CellType.TRIANGLE, pv.CellType.QUAD]
VOLUME_TYPES = [pv.CellType.TETRA, pv.CellType.HEXAHEDRON]


def _convert_to_mesh(model: pv.UnstructuredGrid) -> Mesh:
    """Convert to model to Mesh object."""
    if isinstance(model, pv.PolyData):
        model = model.cast_to_unstructured_grid()

    model.cell_data["_volume-id"] = np.array(1, dtype=float) * np.nan
    model.cell_data["_surface-id"] = np.array(1, dtype=float) * np.nan

    mask = np.isin(model.celltypes, SURFACE_TYPES)
    model.cell_data["_surface-id"][mask] = np.array(1, dtype=float)

    mask = np.isin(model.celltypes, VOLUME_TYPES)
    model.cell_data["_volume-id"][mask] = np.array(10, dtype=float)

    return Mesh(model)


# define different beam models that can be used for testing.
def _get_beam_model(
    cell_type: Literal["tets", "tets+triangles", "triangles", "hex", "hex+quads", "quads"],
) -> pv.UnstructuredGrid | pv.PolyData:
    """Generates various beam models.

    Parameters
    ----------
    cell_type : Literal[&quot;tets&quot;, &quot;tets
        Cell type for beam model to consist of.

    Returns
    -------
    pv.UnstructuredGrid | pv.PolyData
        Beam model of defined by cells of type cell_type in UnstructuredGrid or PolyData form.
    """
    from pyvista import examples

    if cell_type == "triangles":
        testmodel = examples.load_tetbeam().extract_surface()
        return testmodel

    if "tets" in cell_type:
        testmodel = examples.load_tetbeam()
        if "triangles" in cell_type:
            testmodel = testmodel + testmodel.extract_surface()
        return testmodel

    if cell_type == "quads":
        testmodel = examples.load_hexbeam().extract_surface()
        return testmodel

    if "hex" in cell_type:
        testmodel = examples.load_hexbeam()
        if "quads" in cell_type:
            testmodel = testmodel + testmodel.extract_surface()
        return testmodel


@pytest.mark.parametrize("dtype", [float, int])
def test_pyvista_clean_grid(dtype):
    import pyvista as pv

    if dtype is int:
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

    # prep data based on beam model
    tets = _get_beam_model("tets")
    tets.clear_cell_data()

    edges = tets.extract_feature_edges()
    edges.clear_cell_data()
    edges.clear_point_data()

    triangles = _get_beam_model("triangles")
    triangles.clear_cell_data()
    triangles.clear_point_data()

    # init the base mesh.
    base = Mesh(tets)

    merged = base._add_mesh(triangles)
    assert isinstance(merged, Mesh)
    assert base.n_cells == triangles.n_cells + tets.n_cells

    # assert adding edges.
    base._add_mesh(edges)
    assert base.n_cells == triangles.n_cells + tets.n_cells + edges.n_cells

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
    mesh = _convert_to_mesh(_get_beam_model("tets"))
    mesh.cell_data["_volume-id"] = float(1)
    init_n_cells = mesh.n_cells

    surface = _get_beam_model("triangles")
    # check adding without id or id as float
    assert mesh.add_surface(surface) is None
    assert mesh.add_surface(surface, float(1)) is None

    # test adding a single surface
    mesh.add_surface(surface, id=2)
    assert "_surface-id" in mesh.cell_data.keys()

    assert np.all(np.isin(mesh.celltypes, [pv.CellType.TRIANGLE, pv.CellType.TETRA]))

    expected_volume_ids = np.hstack([init_n_cells * [1], surface.n_cells * [np.nan]])
    expected_surface_ids = np.hstack([init_n_cells * [np.nan], surface.n_cells * [2]])
    np.testing.assert_allclose(mesh.cell_data["_volume-id"], expected_volume_ids)
    np.testing.assert_allclose(mesh.cell_data["_surface-id"], expected_surface_ids)

    # test adding multiple surfaces simultaneously with celldata _surface-id
    mesh = _convert_to_mesh(_get_beam_model("tets"))
    mesh.cell_data["_volume-id"] = float(1)

    triangles = _get_beam_model("triangles")
    quads = _get_beam_model("quads")
    triangles.cell_data["_surface-id"] = float(10)
    quads.cell_data["_surface-id"] = float(11)

    surface_to_add = triangles + quads

    mesh.add_surface(surface_to_add)
    np.testing.assert_allclose(np.unique(mesh.cell_data["_surface-id"]), [10, 11, np.nan])

    # check behavior when the same id is specified.
    mesh = _convert_to_mesh(_get_beam_model("tets"))

    mesh.add_surface(triangles, 10)
    assert mesh.add_surface(quads, 10) is None

    mesh.add_surface(quads, 10, overwrite_existing=True)
    # Both quads and triangles will exist:
    assert np.all(
        np.isin(
            [pv.CellType.TRIANGLE, pv.CellType.QUAD],
            mesh.get_surface(10).cast_to_unstructured_grid().celltypes,
        )
    )
    mesh.add_surface(triangles, 20, name="triangles1")
    triangles1 = mesh.get_surface_by_name("triangles1")
    assert triangles.n_cells == triangles1.n_cells
    assert triangles.n_points == triangles1.n_points
    assert isinstance(triangles1, SurfaceMesh)
    assert triangles1.name == "triangles1"


def test_lines_add_001():
    """Test adding a beam to a Mesh."""
    mesh = _convert_to_mesh(_get_beam_model("tets"))
    mesh.cell_data["_volume-id"] = float(1)

    line = pv.Line([0, 0, 0], [-1, 0, 0])

    assert mesh.add_lines(line) is None
    assert mesh.add_lines(line, float(1)) is None

    mesh.add_lines(line, id=2)

    assert "_line-id" in mesh.cell_data.keys()
    assert np.all(np.isin(mesh.celltypes, [pv.CellType.LINE, pv.CellType.TETRA]))
    np.testing.assert_allclose(np.unique(mesh.cell_data["_volume-id"]), [1, np.nan])
    np.testing.assert_allclose(np.unique(mesh.cell_data["_line-id"]), [2, np.nan])

    # test adding by name

    mesh.add_lines(line, id=3, name="lines1")
    assert mesh.line_names == ["lines1"]
    assert mesh.get_lines_by_name("lines1").n_cells == line.n_cells

    mesh.add_lines(line, id=4, name="lines2")
    assert mesh.line_names == ["lines1", "lines2"]
    assert mesh.get_lines_by_name("lines2").n_cells == line.n_cells
    assert mesh._line_id_to_name == {3: "lines1", 4: "lines2"}


def test_volume_add_001():
    """Test adding a volume (hex elements) to an existing mesh."""
    tets = _convert_to_mesh(_get_beam_model("tets"))
    tets.cell_data["_volume-id"] = float(1)
    n_tets = tets.n_cells

    hex = _get_beam_model("hex")

    assert tets.add_volume(hex) is None
    assert tets.add_volume(hex, float(1)) is None

    tets.add_volume(hex, id=2)
    expected = np.hstack(
        [
            [1] * n_tets,
            [2] * hex.n_cells,
        ]
    )
    assert np.allclose(tets.cell_data["_volume-id"], expected)

    tets.add_volume(hex, id=3, name="hex")
    hex1 = tets.get_volume_by_name("hex")
    assert hex1.n_cells == hex.n_cells
    assert hex1.n_points == hex.n_points


def test_get_submesh_001():
    """Test getting a submesh of Mesh."""
    tets = _get_beam_model("tets")
    hexs = _get_beam_model("hex")
    triangles = _get_beam_model("triangles")
    lines = triangles.extract_all_edges()
    lines.clear_cell_data()
    lines.clear_point_data()

    mesh = Mesh()
    mesh.add_volume(tets, 1)
    mesh.add_surface(triangles, 10)
    mesh.add_lines(lines, 100)

    # test get surfaces
    triangles1 = mesh.get_surface(10)
    assert triangles.n_cells == triangles1.n_cells
    assert np.allclose(triangles.points, triangles1.points)
    assert isinstance(triangles1, pv.PolyData)

    # test get lines
    lines1 = mesh.get_lines(100)
    assert lines.n_cells == lines1.n_cells
    assert lines.n_points == lines1.n_points

    # test get multiple surfaces.
    mesh.add_surface(triangles, 11)
    triangles1 = mesh.get_surface([10, 11])
    assert triangles1.n_cells == triangles.n_cells * 2

    # test get volume
    mesh = _convert_to_mesh(_get_beam_model("tets"))
    mesh.cell_data["_volume-id"] = 1
    hexs.cell_data["_volume-id"] = 2
    mesh.add_volume(hexs)
    assert np.allclose(mesh.volume_ids, [1, 2])

    volume1 = mesh.get_volume(2)
    assert volume1.n_cells == hexs.n_cells


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
    mesh = _convert_to_mesh(_get_beam_model("tets"))
    mesh.cell_data["_volume-id"] = 1
    init_n_cells = mesh.n_cells

    surface = _get_beam_model("triangles")

    # add a single surface
    mesh.add_surface(surface, id=2)
    # remove this surface
    mesh.remove_surface(sid=2)
    assert mesh.n_cells == init_n_cells


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

    # test ignoring nans in point data average works
    grid.point_data["data"] = float(1)
    grid.point_data["data"][0:4] = np.nan

    assert np.all(np.isnan(grid.clean(ignore_nans_in_point_average=False).point_data["data"]))
    assert np.allclose(grid.clean(ignore_nans_in_point_average=True).point_data["data"], 1)

    pass


def test_mesh_object_properties():
    """Test the properties of Mesh."""
    tets = _get_beam_model("tets")
    hex = _get_beam_model("hex")
    triangles = _get_beam_model("triangles")
    lines = triangles.extract_feature_edges()
    lines.clear_cell_data()
    lines.clear_point_data()

    mesh = Mesh()
    mesh.add_volume(tets, int(1))
    mesh.add_volume(hex, int(2))
    mesh.add_surface(triangles, int(10))
    mesh.add_lines(lines, int(100))

    assert np.all(
        np.isin(
            mesh.celltypes,
            [pv.CellType.TETRA, pv.CellType.HEXAHEDRON, pv.CellType.TRIANGLE, pv.CellType.LINE],
        )
    )
    assert np.all(mesh.volume_ids == [1, 2])
    assert np.all(mesh.surface_ids == [10])
    assert np.all(mesh.line_ids == [100])

    pass


def test_surface_mesh_init():
    """Test whether SurfaceMesh __init__'s behaves similar to pv.PolyData init."""

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
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
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
    assert mesh._surfaces == []
    mesh.add_surface(polydata)
    assert len(mesh._surfaces) == 2


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


def test_mesh_id_to_name():
    """Test mapping id to and from volume/surface name."""
    tets = _get_beam_model("tets")
    triangles = _get_beam_model("triangles")
    hex = _get_beam_model("hex")
    quads = _get_beam_model("quads")

    # initialize mesh and add surfaces and volumes
    mesh = Mesh()
    mesh.add_surface(triangles, int(1))
    mesh.add_surface(quads, int(2))
    mesh.add_volume(tets, int(10))
    mesh.add_volume(hex, int(11))

    mesh._surface_id_to_name = {1: "triangles", 2: "quads"}
    mesh._volume_id_to_name = {10: "tets", 11: "hex"}

    assert mesh._surface_name_to_id == {"triangles": 1, "quads": 2}
    assert mesh._volume_name_to_id == {"tets": 10, "hex": 11}

    assert mesh.validate_ids_to_name_map()

    triangles1 = mesh.get_surface_by_name("triangles")
    assert triangles.n_cells == triangles1.n_cells
    assert triangles.n_points == triangles1.n_points
    assert isinstance(triangles1, SurfaceMesh)
    assert triangles1.name == "triangles"
    assert triangles1.id == 1
    assert mesh.get_surface_by_name("silly-name") is None

    tets1 = mesh.get_volume_by_name("tets")
    assert tets.n_cells == tets1.n_cells
    assert tets.n_points == tets1.n_points
    assert mesh.get_volume_by_name("silly-name") is None

    del mesh._volume_id_to_name[10]
    assert mesh._get_unmapped_volumes() == [10]
    assert not mesh.validate_ids_to_name_map()


def test_beamsmesh_add_lines():
    """Test behavior of beamsmesh."""
    mesh = Mesh()
    lines1 = pv.Line([0, 0, 0], [-1, 0, 0])
    lines2 = pv.Line([0, 0, 0], [1, 0, 0])
    mesh.add_lines(lines1, id=1, name="lines1")
    mesh.add_lines(lines2, id=2, name="lines2")

    assert mesh.line_names == ["lines1", "lines2"]
    assert mesh._line_id_to_name == {1: "lines1", 2: "lines2"}
    assert np.allclose(mesh.get_lines(1).points, lines1.points)


def test_mesh_save_load():
    """Test saving the mesh object."""
    mesh = _convert_to_mesh(_get_beam_model("tets+triangles"))

    mesh._surface_id_to_name[1] = "triangles"
    mesh._volume_id_to_name[10] = "tets"

    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tmpdir:
        filename = os.path.join(tmpdir, "test_file.vtu")
        filename_map = filename.replace(".vtu", ".namemap.json")
        mesh.save(filename)

        assert os.path.isfile(filename)
        assert os.path.isfile(filename_map)

        # try to load the same mesh.
        mesh1 = Mesh()
        mesh1.load_mesh(filename)
        assert mesh.n_cells == mesh1.n_cells
        assert mesh.n_points == mesh1.n_points
        assert mesh._surface_id_to_name == mesh1._surface_id_to_name
        assert mesh._volume_id_to_name == mesh1._volume_id_to_name

        # try to load when no json is found
        mesh1 = Mesh()
        os.remove(filename_map)
        mesh1.load_mesh(filename)
        assert mesh.n_cells == mesh1.n_cells
        assert mesh.n_points == mesh1.n_points
        assert mesh1._surface_id_to_name == {}
        assert mesh1._volume_id_to_name == {}
        assert not mesh1.validate_ids_to_name_map()

    return


def test_force_normals_inwards():
    """Test whether compute_volume enforces inwards pointing normals."""
    # by default normals of sphere are pointing outwards.
    import copy

    sphere: pv.PolyData = pv.Sphere()
    sphere = sphere.compute_normals()
    sphere1 = SurfaceMesh(copy.deepcopy(sphere))

    sphere1.force_normals_inwards()

    assert np.allclose(sphere1.cell_data["Normals"], -1 * sphere.cell_data["Normals"])
    pass


def test_cap_properties():
    """Test getting global_node_ids_edge from Cap."""
    from ansys.health.heart.objects import Cap
    from ansys.health.heart.utils.vtk_utils import get_patches_with_centroid

    half_sphere = pv.Sphere().clip(normal="y")
    patches = get_patches_with_centroid(half_sphere)
    patch_mesh = patches[0].clean()
    patch_mesh.point_data["_global-point-ids"] = np.arange(0, patch_mesh.n_points) + 10

    cap = Cap(name="aortic-valve", cap_type=CapType.AORTIC_VALVE)
    cap._mesh = patch_mesh

    assert cap.global_node_ids_edge.shape[0] == cap._mesh.n_points - 1
    assert cap.centroid.shape == (3,)
    assert np.allclose(cap.centroid, [0, 0, 0], atol=1e-7)
    assert np.allclose(
        cap.global_node_ids_edge,
        np.array(
            [
                10,
                11,
                13,
                14,
                15,
                16,
                17,
                18,
                19,
                20,
                21,
                22,
                23,
                24,
                25,
                26,
                27,
                28,
                29,
                30,
                31,
                32,
                33,
                34,
                35,
                36,
                37,
                38,
                39,
                40,
                41,
                42,
                43,
                44,
                45,
                46,
                47,
                48,
                49,
                50,
                51,
                52,
                53,
                54,
                55,
                56,
                57,
                58,
                59,
                60,
                61,
                62,
                63,
                64,
                65,
                66,
                67,
                68,
            ]
        ),
    )
