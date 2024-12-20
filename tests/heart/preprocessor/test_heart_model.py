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

"""Tests shigh-level heart model class."""

import json
import os
import tempfile

import pyvista as pv
from pyvista import examples

from ansys.heart.core.objects import Mesh, PartType

if os.getenv("GITHUB_ACTIONS"):
    is_gh_action = True
else:
    is_gh_action = False

from pathlib import Path
import unittest.mock as mock

import numpy as np
import pytest

import ansys.heart.core.models as models


def test_dump_model_001():
    """Test dumping of model to disk."""

    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:
        model = models.BiVentricle(working_directory=workdir)

        expected_path = os.path.join(model.workdir, "heart_model.pickle")

        model.dump_model()
        assert os.path.isfile(expected_path)

        expected_path = os.path.join(model.workdir, "heart_model1.pickle")
        model.dump_model(expected_path)
        assert os.path.isfile(expected_path)

        expected_path = Path(os.path.join(model.workdir, "heart_model2.pickle"))
        model.dump_model(expected_path)
        assert os.path.isfile(expected_path)


def test_model_load_001():
    """Test dumping and reading of model with data."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:
        model = models.BiVentricle(working_directory=workdir)

        path_to_model = os.path.join(model.workdir, "heart_model.pickle")
        model.left_ventricle.endocardium.triangles = np.array([[0, 1, 2]], dtype=int)
        model.left_ventricle.endocardium.nodes = np.eye(3, 3, dtype=float)

        model.dump_model()

        assert os.path.isfile(path_to_model)

        model_loaded: models.BiVentricle = models.HeartModel.load_model(path_to_model)
        assert np.array_equal(
            model_loaded.left_ventricle.endocardium.triangles,
            model.left_ventricle.endocardium.triangles,
        )
        assert np.array_equal(
            model_loaded.left_ventricle.endocardium.nodes,
            model.left_ventricle.endocardium.nodes,
        )


def test_model_load_002():
    """Test loading model from pickle."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:
        model: models.BiVentricle = models.BiVentricle(working_directory=workdir)
        model.workdir = workdir
        # populate model
        model.left_ventricle.element_ids = np.array([1, 2, 3, 4], dtype=int)
        model.right_ventricle.element_ids = np.array([11, 66, 77, 88], dtype=int)

        model.left_ventricle.endocardium.triangles = np.array([[0, 1, 2], [0, 2, 3]], dtype=int)
        model.left_ventricle.endocardium.nodes = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float
        )

        model.right_ventricle.endocardium.triangles = np.array([[0, 3, 1], [3, 2, 0]], dtype=int)
        model.right_ventricle.endocardium.nodes = (
            np.array(
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float
            )
            + 10
        )

        model.mesh.tetrahedrons = np.array([[0, 1, 2, 3]], dtype=int)
        model.mesh.nodes = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float
        )

        # dump model to disk
        path_to_heart_model = os.path.join(workdir, "heart_model.pickle")
        model.dump_model(path_to_heart_model)

        assert os.path.isfile(path_to_heart_model), "File does not exist"

        # load model
        model1 = models.HeartModel.load_model(path_to_heart_model)

        assert isinstance(model1, models.BiVentricle), "Expecting model of type BiVentricle"

        # compare contents to original
        assert np.array_equal(model1.left_ventricle.element_ids, model.left_ventricle.element_ids)
        assert np.array_equal(model1.right_ventricle.element_ids, model.right_ventricle.element_ids)

        assert np.array_equal(
            model1.left_ventricle.endocardium.triangles, model.left_ventricle.endocardium.triangles
        )
        assert np.array_equal(
            model1.right_ventricle.endocardium.triangles,
            model.right_ventricle.endocardium.triangles,
        )
        assert np.array_equal(model1.mesh.tetrahedrons, model.mesh.tetrahedrons)
        assert np.allclose(model1.mesh.nodes, model.mesh.nodes, atol=1e-8)

    pass


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


def test_load_from_mesh():
    """Test loading mesh from mesh file and id map."""
    # generate a dummy mesh.
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

    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tmpdir:
        mesh_path = os.path.join(tmpdir, "mesh.vtu")

        mesh.save(mesh_path)
        model = models.BiVentricle(working_directory=tmpdir)

        part_info = {
            "Left ventricle": {
                "part-id": 10,
                "part-type": PartType.VENTRICLE.value,
                "surfaces": {"Left ventricle endocardium": 1},
                "caps": {"mitral-valve": 3, "aortic-valve": 4},
                "cavity": {"Left ventricle endocardium": 1},
            },
            "Right ventricle": {
                "part-id": 11,
                "part-type": PartType.VENTRICLE.value,
                "surfaces": {"Right ventricle epicardium": 1},
                "caps": {},
                "cavity": {},
            },
            "Septum": {
                "part-id": 12,
                "part-type": PartType.SEPTUM.value,
                "surfaces": {},
                "caps": {},
                "cavity": {},
            },
        }

        part_info_path = os.path.join(tmpdir, "partinfo.json")
        with open(part_info_path, "w") as f:
            json.dump(part_info, f, indent=4)

        with mock.patch("ansys.heart.core.models.BiVentricle._extract_apex") as mock_extract_apex:
            with mock.patch(
                "ansys.heart.core.models.BiVentricle._define_anatomy_axis"
            ) as mock_define_axis:
                model.load_model_from_mesh(mesh_path, part_info_path)
                mock_extract_apex.assert_called_once()
                mock_define_axis.assert_called_once()

        assert model.part_names == list(part_info.keys())

        assert model.left_ventricle.element_ids.shape[0] == examples.load_tetbeam().n_cells
        assert model.right_ventricle.element_ids.shape[0] == examples.load_tetbeam().n_cells
        assert model.septum.element_ids.shape[0] == examples.load_tetbeam().n_cells

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


def test_heart_model_add_purkinje_from_file():
    """Test adding a purkinje network from a file."""
    lines = """*KEYWORD
*NODE +
              123932  0.117289119795E+03  0.244318845738E+02  0.883234881118E+01
              123933  0.116478727414E+03  0.250487507439E+02  0.795046863408E+01
              123934  0.116515199453E+03  0.236812534885E+02  0.963833959240E+01
*ELEMENT_BEAM
  560893       8  123932  123933
  560894       8  123932  123934
*END"""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        model = models.HeartModel(tempdir)
        temp_path = os.path.join(tempdir, "purkinje1.k")
        with open(temp_path, "w") as f:
            f.writelines(lines)
        model.add_purkinje_from_kfile(temp_path, "purkinje1")

        assert len(model.beam_network) == 1
        assert model.beam_network[0].pid == 8
        assert model.beam_network[0].name == "purkinje1"
        assert np.all(model.beam_network[0].edges == np.array([[0, 1], [0, 2]]))
        assert np.allclose(
            model.beam_network[0].points,
            np.array(
                [
                    [0.117289119795e03, 0.244318845738e02, 0.883234881118e01],
                    [0.116478727414e03, 0.250487507439e02, 0.795046863408e01],
                    [0.116515199453e03, 0.236812534885e02, 0.963833959240e01],
                ]
            ),
        )
