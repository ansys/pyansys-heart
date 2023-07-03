"""Test input manager"""

import os
from typing import Union

from ansys.heart.preprocessor.input import InputManager
import numpy as np
import pyvista as pv

import tests.heart.conftest as conftest
from tests.heart.preprocessor.common import _create_simple_unstructured_grid


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


def test_inputs():
    """Test the different kind of inputs. UnstructuredGrid and PolyData file/object."""
    # unstructured grid as file or object

    workdir = conftest.get_workdir()

    # UnstructuredGrid as object
    ugrid = _create_simple_unstructured_grid()

    input = InputManager(ugrid, scalar="tags")
    assert _is_same_mesh(input.input_volume, ugrid)

    # UnstructuredGrid as vtu file
    file1 = os.path.join(workdir, "ugrid.vtu")
    ugrid.save(file1)
    input = InputManager(file1, scalar="tags")
    assert _is_same_mesh(input.input_volume, ugrid)

    # Prep PolyData input.
    polydata = pv.Sphere()
    polydata.cell_data.set_scalars(np.full(polydata.n_cells, 1), "surface-tags")

    # Polydata as object
    input = InputManager(polydata, scalar="surface-tags")
    assert _is_same_mesh(input.input_boundary, polydata)

    # Polydata as vtp file
    file2 = os.path.join(workdir, "polydata.vtp")
    polydata.save(file2)
    input = InputManager(file2, scalar="surface-tags")
    assert _is_same_mesh(input.input_boundary, polydata)

    os.remove(file1)
    os.remove(file2)


def test_reorder_ids():
    """Test reordering surface/part ids."""
    from ansys.heart.preprocessor.input import (
        _get_boundary_name_to_boundary_id_map,
        _get_part_name_to_part_id_map,
    )

    reference_map = _get_part_name_to_part_id_map()
    ugrid = _create_simple_unstructured_grid()
    ugrid.rename_array("tags", "part-id")
    ugrid.cell_data["part-id"] = [1, 2]
    part_name_to_part_id = {"Right ventricle myocardium": 1, "Left ventricle myocardium": 2}

    input = InputManager(ugrid, name_to_id_map=part_name_to_part_id)

    assert np.all(np.equal(input.input_volume.cell_data["part-id"], [2, 1]))

    reference_map = _get_boundary_name_to_boundary_id_map()
    test_map = {
        "left-ventricle-endocardium": 2,
        "left-ventricle-epicardium": 4,
        "interface@left-ventricle_aortic-valve": 1,
        "interface@left-ventricle_mitral-valve": 3,
    }
    # Prep PolyData input.
    cube = pv.Cube()
    cube.cell_data.set_scalars([1, 1, 2, 2, 3, 4], "boundary-id")

    # ref ids
    ref_ids = [
        reference_map["interface@left-ventricle_aortic-valve"],
        reference_map["interface@left-ventricle_aortic-valve"],
        reference_map["left-ventricle-endocardium"],
        reference_map["left-ventricle-endocardium"],
        reference_map["interface@left-ventricle_mitral-valve"],
        reference_map["left-ventricle-epicardium"],
    ]

    input = InputManager(cube, name_to_id_map=test_map)

    assert np.all(np.equal(input.input_boundary.cell_data["boundary-id"], ref_ids))

    return


def test_reorder_part_ids_001():
    """Test reorder of part-ids given a part-name to part-id dictionary."""

    from ansys.heart.preprocessor.input import _get_part_name_to_part_id_map

    reference_map = _get_part_name_to_part_id_map()
    ugrid = _create_simple_unstructured_grid()
    ugrid.rename_array("tags", "part-id")
    ugrid.cell_data["part-id"] = [1, 2]
    input = InputManager(ugrid)

    # first test case.
    part_name_to_part_id = {"Right ventricle myocardium": 1, "Left ventricle myocardium": 2}
    input.input_volume.cell_data["part-id"] = [1, 2]
    input._reorder_part_ids(part_name_to_part_id)

    assert np.all(
        input.input_volume.cell_data["part-id"]
        == [reference_map["Right ventricle myocardium"], reference_map["Left ventricle myocardium"]]
    )

    # second test case.
    input.input_volume.cell_data["part-id"] = [1, 2]
    input._reorder_part_ids({"Aorta wall": 1, "Pulmonary artery wall": 2})

    assert np.all(
        input.input_volume.cell_data["part-id"]
        == [reference_map["Aorta wall"], reference_map["Pulmonary artery wall"]]
    )

    # third test case. With "unknown" name
    input.input_volume.cell_data["part-id"] = [1, 2]
    input._reorder_part_ids({"Aorta wall": 1, "Randomly named zone": 2})
    assert np.all(
        input.input_volume.cell_data["part-id"]
        == [reference_map["Aorta wall"], max(reference_map.values()) + 1]
    )

    pass


def test_reorder_boundary_ids_001():
    """Test reordering of boundary-ids given a surface-name to boundary-id dictionary."""
    from ansys.heart.preprocessor.input import _get_boundary_name_to_boundary_id_map

    reference_map = _get_boundary_name_to_boundary_id_map()

    test_map = {
        "left-ventricle-endocardium": 2,
        "left-ventricle-epicardium": 4,
        "interface@left-ventricle_aortic-valve": 1,
        "interface@left-ventricle_mitral-valve": 3,
    }
    # Prep PolyData input.
    cube = pv.Cube()
    cube.cell_data.set_scalars([1, 1, 2, 2, 3, 4], "boundary-id")

    # ref ids
    ref_ids = [
        reference_map["interface@left-ventricle_aortic-valve"],
        reference_map["interface@left-ventricle_aortic-valve"],
        reference_map["left-ventricle-endocardium"],
        reference_map["left-ventricle-endocardium"],
        reference_map["interface@left-ventricle_mitral-valve"],
        reference_map["left-ventricle-epicardium"],
    ]

    input = InputManager(cube)

    input._reorder_boundary_ids(test_map)
    assert np.all(input.input_boundary.cell_data["boundary-id"] == ref_ids)


def test_reorder_boundary_ids_002():
    """Test reordering of boundary-ids when unknown surface name given."""
    from ansys.heart.preprocessor.input import _get_boundary_name_to_boundary_id_map

    reference_map = _get_boundary_name_to_boundary_id_map()

    test_map = {
        "left-ventricle-endocardium": 2,
        "my-random-surface-name": 4,
        "interface@left-ventricle_aortic-valve": 1,
        "interface@left-ventricle_mitral-valve": 3,
    }
    # Prep PolyData input.
    cube = pv.Cube()
    cube.cell_data.set_scalars([1, 1, 2, 2, 3, 4], "boundary-id")

    # ref ids
    ref_ids = [
        reference_map["interface@left-ventricle_aortic-valve"],
        reference_map["interface@left-ventricle_aortic-valve"],
        reference_map["left-ventricle-endocardium"],
        reference_map["left-ventricle-endocardium"],
        reference_map["interface@left-ventricle_mitral-valve"],
        max(reference_map.values()) + 1,
    ]

    input = InputManager(cube)

    input._reorder_boundary_ids(test_map)

    assert np.all(input.input_boundary.cell_data["boundary-id"] == ref_ids)


def test_boundary_export():
    """Test the export of surfaces."""

    # Prep PolyData input.
    cube = pv.Cube()
    cube.cell_data.set_scalars([1, 1, 2, 2, 3, 4], "boundary-id")

    input = InputManager(cube)
    input.export_boundaries(".stl", folder=conftest.get_workdir())
