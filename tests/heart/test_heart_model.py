# Copyright (C) 2023 - 2026 ANSYS, Inc. and/or its affiliates.
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

"""Tests shigh-level heart model class."""

import json
import os
import tempfile
import unittest.mock as mock

import numpy as np
import pytest
import pyvista as pv
from pyvista import examples

import ansys.health.heart.models as models
from ansys.health.heart.objects import Mesh
from ansys.health.heart.parts import _PartType
from ansys.health.heart.settings.material.cell_models import Tentusscher


def test_set_workdir(monkeypatch):
    """Test setting the working directory."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        workdir = models._set_workdir(os.path.join(tempdir, "test"))
        assert workdir == os.path.join(tempdir, "test")
        assert os.path.isdir(workdir)

        # test setting workdir to current workdir
        with mock.patch("os.getcwd") as mock_getcwd:
            mock_getcwd.return_value = os.path.join(tempdir, "test1")
            workdir = models._set_workdir()
            assert workdir == os.path.join(tempdir, "test1")

        monkeypatch.setenv("PYANSYS_HEART_WORKDIR", os.path.join(tempdir, "test2"))
        workdir = models._set_workdir()
        assert workdir == os.path.join(tempdir, "test2")


def test_import_mesher_missing_dependency(monkeypatch):
    """Test that _import_mesher raises a clear error when ansys-fluent-core is missing."""
    import sys

    for key in list(sys.modules):
        if "ansys.health.heart.pre.mesher" in key or "ansys.fluent" in key:
            monkeypatch.delitem(sys.modules, key)

    monkeypatch.setitem(sys.modules, "ansys.fluent", None)
    monkeypatch.setitem(sys.modules, "ansys.fluent.core", None)

    with pytest.raises(ModuleNotFoundError, match="pip install ansys-health-heart\\[meshing\\]"):
        models._import_mesher()


def test_save_model():
    """Test dumping of model to disk."""

    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:
        model = models.BiVentricle(working_directory=workdir)

        expected_mesh_path = os.path.join(os.path.join(workdir, "test.vtu"))
        expected_info_path = os.path.join(os.path.join(workdir, "test.partinfo.json"))

        mesh_path, info_path = model.save_model(os.path.join(workdir, "test.vtu"))

        assert mesh_path == expected_mesh_path
        assert info_path == expected_info_path

        assert os.path.isfile(expected_mesh_path)
        assert os.path.isfile(expected_info_path)


def test_remove_part_updates_part_info():
    """Test that remove_part also removes the part from _part_info."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:
        model = models.BiVentricle(working_directory=workdir)

        # Build _part_info so all parts are registered.
        model._get_parts_info()
        assert "Septum" in model._part_info

        # Remove the septum part.
        model.remove_part("Septum")

        # _part_info should no longer contain the removed part.
        assert "Septum" not in model._part_info

        # Saving and reloading should not bring the part back.
        mesh_path, info_path = model.save_model(os.path.join(workdir, "test.vtu"))
        with open(info_path, "r") as f:
            saved_info = json.load(f)
        assert "Septum" not in saved_info


@pytest.mark.parametrize(
    "model_type,expected_part_names",
    [
        (models.LeftVentricle, ["Left ventricle"]),
        (models.BiVentricle, ["Left ventricle", "Right ventricle", "Septum"]),
        (
            models.FourChamber,
            [
                "Left ventricle",
                "Right ventricle",
                "Septum",
                "Left atrium",
                "Right atrium",
            ],
        ),
        (
            models.FullHeart,
            [
                "Left ventricle",
                "Right ventricle",
                "Septum",
                "Left atrium",
                "Right atrium",
                "Aorta",
                "Pulmonary artery",
            ],
        ),
    ],
)
def test_model_part_names(model_type, expected_part_names):
    """Test whether all parts exist in the model."""
    model: models.HeartModel = model_type()

    assert model.part_names == expected_part_names


def _get_test_model():
    #! Note, can modify to create something more meaningful,
    #! e.g. a sphere with inner/outer surface and caps.
    mesh = Mesh()
    mesh.add_surface(pv.Sphere(), int(1))
    mesh.add_surface(pv.Box(), int(2))
    mesh.add_surface(pv.Disc(), int(3))
    mesh.add_surface(pv.Disc().translate((1, 0, 0)), int(4))
    mesh.add_surface(pv.Sphere(radius=0.3), int(5))

    mesh.add_volume(examples.load_tetbeam(), int(10))
    mesh.add_volume(examples.load_tetbeam(), int(11))
    mesh.add_volume(examples.load_tetbeam(), int(12))

    mesh._volume_id_to_name[10] = "Left ventricle"
    mesh._volume_id_to_name[11] = "Right ventricle"
    mesh._volume_id_to_name[12] = "Septum"

    mesh._surface_id_to_name[1] = "Left ventricle endocardium"
    mesh._surface_id_to_name[2] = "Right ventricle epicardium"
    mesh._surface_id_to_name[3] = "mitral-valve"
    mesh._surface_id_to_name[4] = "aortic-valve"
    mesh._surface_id_to_name[5] = "Left ventricle cavity"

    part_info = {
        "Left ventricle": {
            "part-id": 10,
            "part-type": _PartType.VENTRICLE.value,
            "surfaces": {"Left ventricle endocardium": 1},
            "caps": {"mitral-valve": 3, "aortic-valve": 4},
            "cavity": {"Left ventricle endocardium": 1},
        },
        "Right ventricle": {
            "part-id": 11,
            "part-type": _PartType.VENTRICLE.value,
            "surfaces": {"Right ventricle epicardium": 2},
            "caps": {},
            "cavity": {},
        },
        "Septum": {
            "part-id": 12,
            "part-type": _PartType.SEPTUM.value,
            "surfaces": {},
            "caps": {},
            "cavity": {},
        },
    }

    return mesh, part_info


def test_load_model_arbitrary_part():
    """Test loading model with an arbitrary part."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tmpdir:
        mesh_path = os.path.join(tmpdir, "mesh.vtu")

        mesh, part_info = _get_test_model()

        mesh.save(mesh_path)

        part_info_path = os.path.join(tmpdir, "partinfo.json")

        # Add an arbitrary part to the part info.
        part_info["ArbitraryPart"] = part_info["Septum"]
        part_info["ArbitraryPart"]["part-id"] = 14
        part_info["ArbitraryPart"]["part-type"] = "undefined"

        with open(part_info_path, "w") as f:
            json.dump(part_info, f, indent=4)

        model = models.HeartModel.load_model(mesh_path, part_info_path, tmpdir)

        assert isinstance(model, models.BiVentricle)
        assert "ArbitraryPart" in model.part_names
        assert "arbitrarypart" in model.__dict__.keys()


def test_load_model():
    """Test loading model from a mesh file and part information file."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tmpdir:
        mesh_path = os.path.join(tmpdir, "mesh.vtu")

        mesh, part_info = _get_test_model()

        mesh.save(mesh_path)

        part_info_path = os.path.join(tmpdir, "partinfo.json")
        with open(part_info_path, "w") as f:
            json.dump(part_info, f, indent=4)

        with mock.patch("ansys.health.heart.models.HeartModel._extract_apex") as mock_extract_apex:
            with mock.patch(
                "ansys.health.heart.models.HeartModel._define_anatomy_axis"
            ) as mock_define_axis:
                model = models.HeartModel.load_model(mesh_path, part_info_path, tmpdir)

                assert isinstance(model, models.BiVentricle)
                assert model.mesh.n_cells == mesh.n_cells
                assert model.mesh.n_points == mesh.n_points

                assert model._part_info == part_info
                mock_extract_apex.assert_called_once()
                mock_define_axis.assert_called_once()

        assert model.part_names == list(part_info.keys())

        assert (
            model.left_ventricle.get_element_ids(model.mesh).shape[0]
            == examples.load_tetbeam().n_cells
        )
        assert (
            model.right_ventricle.get_element_ids(model.mesh).shape[0]
            == examples.load_tetbeam().n_cells
        )
        assert model.septum.get_element_ids(model.mesh).shape[0] == examples.load_tetbeam().n_cells

        assert model.left_ventricle.endocardium.n_cells == pv.Sphere().n_cells
        assert model.left_ventricle.endocardium.n_points == pv.Sphere().n_points

        assert model.right_ventricle.epicardium.n_cells == pv.Box().n_cells
        assert model.right_ventricle.epicardium.n_points == pv.Box().n_points

        assert model.left_ventricle.caps[0]._mesh.n_cells == pv.Disc().n_cells
        assert model.left_ventricle.caps[1]._mesh.n_cells == pv.Disc().n_cells

        assert model.left_ventricle.cavity.surface.n_cells == pv.Sphere().n_cells
        assert model.left_ventricle.cavity.surface.n_points == pv.Sphere().n_points

    return


def test_load_model_from_mesh():
    """Test loading mesh from mesh file and ID map."""
    # generate a dummy mesh.

    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tmpdir:
        mesh_path = os.path.join(tmpdir, "mesh.vtu")

        mesh, part_info = _get_test_model()

        mesh.save(mesh_path)

        part_info_path = os.path.join(tmpdir, "partinfo.json")
        with open(part_info_path, "w") as f:
            json.dump(part_info, f, indent=4)

        with mock.patch("ansys.health.heart.models.BiVentricle._extract_apex") as mock_extract_apex:
            with mock.patch(
                "ansys.health.heart.models.BiVentricle._define_anatomy_axis"
            ) as mock_define_axis:
                model = models.BiVentricle.load_model(mesh_path, part_info_path, tmpdir)
                mock_extract_apex.assert_called_once()
                mock_define_axis.assert_called_once()

        assert model.part_names == list(part_info.keys())

        assert (
            model.left_ventricle.get_element_ids(model.mesh).shape[0]
            == examples.load_tetbeam().n_cells
        )
        assert (
            model.right_ventricle.get_element_ids(model.mesh).shape[0]
            == examples.load_tetbeam().n_cells
        )
        assert model.septum.get_element_ids(model.mesh).shape[0] == examples.load_tetbeam().n_cells

        assert model.left_ventricle.endocardium.n_cells == pv.Sphere().n_cells
        assert model.left_ventricle.endocardium.n_points == pv.Sphere().n_points

        assert model.right_ventricle.epicardium.n_cells == pv.Box().n_cells
        assert model.right_ventricle.epicardium.n_points == pv.Box().n_points

        assert model.left_ventricle.caps[0]._mesh.n_cells == pv.Disc().n_cells
        assert model.left_ventricle.caps[1]._mesh.n_cells == pv.Disc().n_cells

        assert model.left_ventricle.cavity.surface.n_cells == pv.Sphere().n_cells
        assert model.left_ventricle.cavity.surface.n_points == pv.Sphere().n_points

        # TODO: Could add more asserts, e.g. for the cavity.
    return


def test_model_get_set_axes():
    """Test getting and setting the axis of a model."""
    model = models.HeartModel()
    model.mesh = pv.Sphere().cast_to_unstructured_grid()
    test_axis = {"center": [0.0, 0.0, 0.0], "normal": [1.0, 0.0, 0.0]}

    assert model.l2cv_axis is None
    assert model.l4cv_axis is None
    assert model.short_axis is None

    model.l2cv_axis = test_axis
    model.l4cv_axis = test_axis
    model.short_axis = test_axis

    assert len(set(model.l2cv_axis) - set(test_axis)) == 0
    assert len(set(model.l4cv_axis) - set(test_axis)) == 0
    assert len(set(model.short_axis) - set(test_axis)) == 0


def test_create_stiff_ventricle_base():
    """Test creating a stiff ventricle base."""
    model = models.BiVentricle()
    # create a test mesh (apico-basal coordinate defined by z-coordinate)
    mesh = examples.load_tetbeam()
    mesh.cell_data["_volume-id"] = 1
    mesh.point_data["apico-basal"] = mesh.points[:, 2] / 5
    total_num_cells = mesh.n_cells

    mesh1 = Mesh()
    mesh1.add_volume(mesh, id=1, name="Left ventricle")
    model.mesh = mesh1
    model.left_ventricle.pid = 1

    part = model.create_stiff_ventricle_base()

    part_element_ids = part.get_element_ids(model.mesh)
    assert len(part_element_ids) == 20
    assert np.all(part_element_ids == np.arange(180, 200))
    # Check whether the left ventricle part has the correct number of cells
    assert len(model.left_ventricle.get_element_ids(model.mesh)) == total_num_cells - 20


class TestNodesetCellmodel:
    """Unit tests for the HeartModel._nodeset_cellmodel property."""

    def setup_method(self):
        """Create a fresh HeartModel for each test."""
        self.model = models.HeartModel()

    def _make_valid(self, n=2):
        """Return a valid (node_sets, cell_models) tuple of length n."""
        node_sets = [np.array([0, 1, 2]) for _ in range(n)]
        cell_models = [Tentusscher() for _ in range(n)]
        return node_sets, cell_models

    # --- default state ---

    def test_default_is_none(self):
        """_nodeset_cellmodel is None on a fresh model."""
        assert self.model._nodeset_cellmodel is None

    # --- valid assignments ---

    def test_valid_assignment(self):
        """A correctly typed tuple is accepted and stored."""
        node_sets, cell_models = self._make_valid(3)
        self.model._nodeset_cellmodel = (node_sets, cell_models)
        result = self.model._nodeset_cellmodel
        assert result is not None
        assert result[0] is node_sets
        assert result[1] is cell_models

    # --- type errors: outer container ---

    def test_raises_if_not_tuple(self):
        """Passing a list instead of tuple raises TypeError."""
        node_sets, cell_models = self._make_valid()
        with pytest.raises(TypeError, match="tuple"):
            self.model._nodeset_cellmodel = [node_sets, cell_models]

    # --- type errors: node_sets element ---

    def test_raises_if_node_set_elements_not_ndarray(self):
        """node_sets containing non-ndarray entries raises TypeError."""
        _, cell_models = self._make_valid(1)
        with pytest.raises(TypeError, match="list of np.ndarray"):
            self.model._nodeset_cellmodel = ([[0, 1, 2]], cell_models)

    # --- length mismatch ---

    def test_raises_if_length_mismatch(self):
        """node_sets and cell_models of different lengths raises ValueError."""
        node_sets, _ = self._make_valid(3)
        _, cell_models = self._make_valid(2)
        with pytest.raises(ValueError, match="same length"):
            self.model._nodeset_cellmodel = (node_sets, cell_models)
