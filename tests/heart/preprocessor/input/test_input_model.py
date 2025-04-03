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

"""Test input manager."""

import glob as glob
import os
import tempfile
import unittest.mock as mock

import numpy as np
import pytest
import pyvista as pv

from ansys.heart.core.preprocessor.input import _InputBoundary, _InputModel, _InputPart


def _get_test_input_model() -> tuple[_InputModel, pv.PolyData, dict]:
    """A test input model using a sphere."""
    polydata = pv.Sphere()
    polydata.cell_data.set_scalars(np.full(polydata.n_cells, 1), "boundary-id")
    polydata.cell_data["boundary-id"][0:100] = 2
    polydata.cell_data["boundary-id"][10:200] = 3

    part_definitions = {
        "sphere1": {
            "id": 1,
            "enclosed_by_boundaries": {"shells1": 1, "shells2": 2},
        },
        "sphere2": {
            "id": 2,
            "enclosed_by_boundaries": {"shells3": 3, "shells4": 2},
        },
    }

    model = _InputModel(
        input=polydata,
        part_definitions=part_definitions,
        scalar="boundary-id",
    )

    return model, polydata, part_definitions


def _is_same_mesh(
    object1: pv.UnstructuredGrid | pv.PolyData,
    object2: pv.UnstructuredGrid | pv.PolyData,
) -> bool:
    """Compares two pyvista objects for mesh correspondence."""
    if type(object1) is not type(object2):
        raise ValueError("Expecting both objects to be of same type.")

    assert np.all(np.isclose(object1.points, object2.points)), "Points of objects not the same"

    if isinstance(object1, pv.UnstructuredGrid):
        assert np.all(object1.cells == object2.cells), "Cells of UnstructuredGrid not the same."

    elif isinstance(object2, pv.PolyData):
        assert np.all(object1.faces == object2.faces), "Faces of PolyData not the same."

    return True


# TODO: extend and refactor initialization.
def test_input_model_init():
    """Test initializing the input model."""

    # specify unsupported input.
    input = pv.Box().cast_to_unstructured_grid()
    with pytest.raises(NotImplementedError):
        _InputModel(input)


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
        scalar="boundary-id",
    )

    # check if part names, part ids, and boundaries are added correctly to the InputModel.
    assert ["sphere1", "sphere2"] == input.part_names
    assert np.all(part_ids == input.part_ids)
    assert ["shells1", "shells2", "shells3", "shells4"] == input.boundary_names
    assert np.all(boundary_ids == input.boundary_ids)

    pass


def test_add_parts_and_boundaries_004():
    """Test adding a single part enclosed by single boundary with multiple surface-ids."""
    # Prep PolyData input.
    polydata = pv.Sphere()
    polydata.cell_data["surface-tags"] = 2
    polydata.cell_data["surface-tags"][0:20] = 3

    # Polydata as object
    input = _InputModel(
        polydata,
        scalar="surface-tags",
        part_definitions={"sphere": {"id": 1, "enclosed_by_boundaries": {"shells": [2, 3]}}},
    )

    # check if part names, part ids, and boundaries are added and modified
    # correctly to the InputModel.
    assert ["sphere"] == input.part_names
    assert [1] == input.part_ids
    assert ["shells"] == input.boundary_names
    assert np.all([2] == input.boundary_ids)

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

    assert input._validate_if_parts_manifold()

    polydata = polydata.remove_cells([1, 2])

    input = _InputModel(
        polydata,
        part_definitions={
            "sphere1": {
                "id": 1,
                "enclosed_by_boundaries": {"shells1": 1},
            },
        },
        scalar="boundary-id",
    )

    assert not input._validate_if_parts_manifold()


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
        scalar="boundary-id",
    )

    with tempfile.TemporaryDirectory() as write_dir:
        model.write_part_boundaries(write_dir, extension=".stl")

        assert os.path.isfile(os.path.join(write_dir, "shells1.stl"))
        assert os.path.isfile(os.path.join(write_dir, "shells2.stl"))
        assert os.path.isfile(os.path.join(write_dir, "shells3.stl"))
        assert not os.path.isfile(os.path.join(write_dir, "shells4.stl"))

        os.remove(os.path.join(write_dir, "shells1.stl"))
        os.remove(os.path.join(write_dir, "shells2.stl"))
        os.remove(os.path.join(write_dir, "shells3.stl"))

    return


def test_write_to_polydata():
    """Test writing to polydata."""

    model = _get_test_input_model()[0]

    assert np.all(
        np.unique(model.as_single_polydata.cell_data["boundary-id"]) == np.array([1, 2, 3])
    )

    with tempfile.TemporaryDirectory() as write_dir:
        model.write_part_boundaries(write_dir)

        files = glob.glob(os.path.join(write_dir, "*.stl"))
        assert len(files) == 3
        assert sorted([os.path.basename(file) for file in files]) == [
            "shells1.stl",
            "shells2.stl",
            "shells3.stl",
        ]

    return


def test_validate_input():
    """Test validate input method."""
    model, _, _ = _get_test_input_model()

    # change id to id that is not present in `boundary-id` array.
    model.parts[0].boundaries[0].id = 100
    with pytest.raises(ValueError):
        model._validate_input()

    model._parts = []
    assert model._validate_input() is None

    pass


def test_input_uniqueness():
    """Test model validator."""

    model, polydata, part_definitions = _get_test_input_model()
    polydata.cell_data["boundary-id"][0:100] = 2
    polydata.cell_data["boundary-id"][10:180] = 3
    polydata.cell_data["boundary-id"][180:200] = 4

    # Not unique due to same boundary id but different name.
    model = _InputModel(
        polydata,
        part_definitions=part_definitions,
        scalar="boundary-id",
    )

    assert not model._validate_uniqueness()

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
        scalar="boundary-id",
    )

    assert not model._validate_uniqueness()

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

    assert model._validate_uniqueness()

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

    assert model._validate_uniqueness()


@mock.patch("pyvista.Plotter.show")
def test_input_model_plot(mock_plot):
    """Test plotting the input model."""
    model = _get_test_input_model()[0]
    model.plot()
    mock_plot.assert_called_once()


def test_input_part():
    """Test part."""
    boundary_name = "sphere"
    boundary = _InputBoundary(pv.Sphere(), name=boundary_name, id=1)

    # if not passed as list raise error.
    with pytest.raises(TypeError) as exception_info:
        _InputPart(name="Part1", id=1, boundaries=boundary)
        assert exception_info.value.args[0] == "Boundaries should be a list."

    part = _InputPart(name="Part1", id=1, boundaries=[boundary])

    assert part.boundary_ids == [1]

    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        part.write_boundaries(tempdir, ".stl")

        assert os.path.isfile(os.path.join(tempdir, "sphere.stl"))
