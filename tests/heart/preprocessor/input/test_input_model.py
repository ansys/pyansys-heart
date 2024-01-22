"""Test input manager"""

import glob as glob
import os
import shutil
from typing import Union

from ansys.heart.preprocessor.input import _InputModel

# from ansys.heart.preprocessor.input import _InputBoundary, _InputPart
import numpy as np
import pyvista as pv

import tests.heart.conftest as conftest

# from tests.heart.preprocessor.common import _create_simple_unstructured_grid


def _is_same_mesh(
    object1: Union[pv.UnstructuredGrid, pv.PolyData],
    object2: Union[pv.UnstructuredGrid, pv.PolyData],
) -> bool:
    """Compares two pyvista objects for mesh correspondence."""
    if type(object1) != type(object2):
        raise ValueError("Expecting both objects to be of same type.")

    assert np.all(np.isclose(object1.points, object2.points)), "Points of objects not the same"

    if isinstance(object1, pv.UnstructuredGrid):
        assert np.all(object1.cells == object2.cells), "Cells of UnstructuredGrid not the same."

    elif isinstance(object2, pv.PolyData):
        assert np.all(object1.faces == object2.faces), "Faces of PolyData not the same."

    return True


def test_add_parts_and_boundaries_001():
    """Test adding a single part enclosed by single boundary."""
    # Prep PolyData input.
    polydata = pv.Sphere()
    part_id = 1
    boundary_id = 2
    polydata.cell_data.set_scalars(np.full(polydata.n_cells, boundary_id), "surface-tags")

    # Polydata as object
    input = _InputModel(
        polydata,
        scalar="surface-tags",
        part_definitions={
            "sphere": {"id": part_id, "enclosed_by_boundaries": {"shells": boundary_id}}
        },
    )

    # check if part names, part ids, and boundaries are added correctly to the InputModel.
    assert ["sphere"] == input.part_names
    assert [part_id] == input.part_ids
    assert ["shells"] == input.boundary_names
    assert np.all([boundary_id] == input.boundary_ids)

    pass


def test_add_parts_and_boundaries_002():
    """Test adding a single part enclosed by multiple boundaries."""
    # Prep PolyData input.
    polydata = pv.Sphere()
    part_id = 1
    boundary_ids = [1, 2]
    polydata.cell_data.set_scalars(np.full(polydata.n_cells, boundary_ids[0]), "surface-tags")
    polydata.cell_data["surface-tags"][-100:] = boundary_ids[1]

    # Polydata as object
    input = _InputModel(
        polydata,
        scalar="surface-tags",
        part_definitions={
            "sphere": {
                "id": part_id,
                "enclosed_by_boundaries": {"shells1": boundary_ids[0], "shells2": boundary_ids[1]},
            }
        },
    )

    # check if part names, part ids, and boundaries are added correctly to the InputModel.
    assert ["sphere"] == input.part_names
    assert [part_id] == input.part_ids
    assert ["shells1", "shells2"] == input.boundary_names
    assert np.all(boundary_ids == input.boundary_ids)

    pass


def test_add_parts_and_boundaries_003():
    """Test adding a multiple parts enclosed by multiple boundaries."""

    # Prep PolyData input.
    polydata1 = pv.Sphere()
    part_ids = [1, 2]
    boundary_ids = [1, 2, 3, 4]
    polydata1.cell_data.set_scalars(np.full(polydata1.n_cells, boundary_ids[0]), "boundary-id")
    polydata1.cell_data["boundary-id"][-100:] = boundary_ids[1]

    # Prep PolyData input.
    polydata2 = pv.Sphere()
    polydata2.points = polydata2.points + [1, 1, 1]
    polydata2.cell_data.set_scalars(np.full(polydata2.n_cells, boundary_ids[2]), "boundary-id")
    polydata2.cell_data["boundary-id"][-100:] = boundary_ids[3]

    polydata_input = polydata1 + polydata2

    # Polydata as object
    input = _InputModel(
        polydata_input,
        part_definitions={
            "sphere1": {
                "id": part_ids[0],
                "enclosed_by_boundaries": {"shells1": boundary_ids[0], "shells2": boundary_ids[1]},
            },
            "sphere2": {
                "id": part_ids[1],
                "enclosed_by_boundaries": {"shells3": boundary_ids[2], "shells4": boundary_ids[3]},
            },
        },
    )

    # check if part names, part ids, and boundaries are added correctly to the InputModel.
    assert ["sphere1", "sphere2"] == input.part_names
    assert np.all(part_ids == input.part_ids)
    assert ["shells1", "shells2", "shells3", "shells4"] == input.boundary_names
    assert np.all(boundary_ids == input.boundary_ids)

    pass


def test_is_manifold():
    """Test the is_manifold flag in part."""
    # Prep PolyData input.
    polydata = pv.Sphere()
    polydata.cell_data.set_scalars(np.full(polydata.n_cells, 1), "boundary-id")

    input = _InputModel(
        polydata,
        part_definitions={
            "sphere1": {
                "id": 1,
                "enclosed_by_boundaries": {"shells1": 1},
            },
        },
    )

    assert input._validate_if_parts_manifold() == True

    polydata = polydata.remove_cells([1, 2])

    input = _InputModel(
        polydata,
        part_definitions={
            "sphere1": {
                "id": 1,
                "enclosed_by_boundaries": {"shells1": 1},
            },
        },
    )

    assert input._validate_if_parts_manifold() == False


def test_write_boundaries_001():
    """Test boundary writing."""
    # Prep PolyData input.
    polydata = pv.Sphere()
    polydata.cell_data.set_scalars(np.full(polydata.n_cells, 1), "boundary-id")
    polydata.cell_data["boundary-id"][0:100] = 2
    polydata.cell_data["boundary-id"][10:200] = 3

    model = _InputModel(
        polydata,
        part_definitions={
            "sphere1": {
                "id": 1,
                "enclosed_by_boundaries": {"shells1": 1, "shells2": 2},
            },
            "sphere2": {
                "id": 2,
                "enclosed_by_boundaries": {"shells3": 3, "shells4": 2},
            },
        },
    )

    write_dir = conftest.get_workdir()
    model.write_part_boundaries(write_dir, extension=".stl")

    assert os.path.isfile(os.path.join(write_dir, "shells1.stl"))
    assert os.path.isfile(os.path.join(write_dir, "shells2.stl"))
    assert os.path.isfile(os.path.join(write_dir, "shells3.stl"))
    assert not os.path.isfile(os.path.join(write_dir, "shells4.stl"))

    os.remove(os.path.join(write_dir, "shells1.stl"))
    os.remove(os.path.join(write_dir, "shells2.stl"))
    os.remove(os.path.join(write_dir, "shells3.stl"))


def test_write_to_polydata():
    """Test writing to polydata."""
    polydata = pv.Sphere()
    polydata.cell_data.set_scalars(np.full(polydata.n_cells, 1), "boundary-id")
    polydata.cell_data["boundary-id"][0:100] = 2
    polydata.cell_data["boundary-id"][10:200] = 3

    model = _InputModel(
        polydata,
        part_definitions={
            "sphere1": {
                "id": 1,
                "enclosed_by_boundaries": {"shells1": 1, "shells2": 2},
            },
            "sphere2": {
                "id": 2,
                "enclosed_by_boundaries": {"shells3": 3, "shells4": 2},
            },
        },
    )

    write_dir = os.path.join(conftest.get_workdir(), "test-write-to-poly")
    if os.path.isdir(write_dir):
        shutil.rmtree(write_dir)
    os.mkdir(write_dir)

    model.write_part_boundaries(write_dir)

    files = glob.glob(os.path.join(write_dir, "*.stl"))
    assert len(files) == 3
    assert sorted([os.path.basename(file) for file in files]) == [
        "shells1.stl",
        "shells2.stl",
        "shells3.stl",
    ]
    return


def test_input_uniqueness():
    """Test model validator."""
    polydata = pv.Sphere()
    polydata.cell_data.set_scalars(np.full(polydata.n_cells, 1), "boundary-id")
    polydata.cell_data["boundary-id"][0:100] = 2
    polydata.cell_data["boundary-id"][10:180] = 3
    polydata.cell_data["boundary-id"][180:200] = 4

    # Not unique due to same boundary id but different name.
    model = _InputModel(
        polydata,
        part_definitions={
            "sphere1": {
                "id": 1,
                "enclosed_by_boundaries": {"shells1": 1, "shells2": 2},
            },
            "sphere2": {
                "id": 2,
                "enclosed_by_boundaries": {"shells3": 3, "shells4": 2},
            },
        },
    )

    assert model._validate_uniqueness() == False

    # Not unique due to same boundary id but same name.
    model = _InputModel(
        polydata,
        part_definitions={
            "sphere1": {
                "id": 1,
                "enclosed_by_boundaries": {"shells1": 1, "shells2": 2},
            },
            "sphere2": {
                "id": 2,
                "enclosed_by_boundaries": {"shells3": 3, "shells2": 4},
            },
        },
    )

    assert model._validate_uniqueness() == False

    # Unique due to unique id and name
    model = _InputModel(
        polydata,
        part_definitions={
            "sphere1": {
                "id": 1,
                "enclosed_by_boundaries": {"shells1": 1, "shells2": 2},
            },
            "sphere2": {
                "id": 2,
                "enclosed_by_boundaries": {"shells3": 3, "shells4": 4},
            },
        },
    )

    assert model._validate_uniqueness() == True

    # Unique due to unique id and name
    model = _InputModel(
        polydata,
        part_definitions={
            "sphere1": {
                "id": 1,
                "enclosed_by_boundaries": {"shells1": 1, "shells2": 2},
            },
            "sphere2": {
                "id": 2,
                "enclosed_by_boundaries": {"shells3": 3, "shells2": 2},
            },
        },
    )

    assert model._validate_uniqueness() == True
