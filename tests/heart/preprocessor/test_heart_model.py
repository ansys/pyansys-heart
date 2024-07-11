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

"""Tests shigh-level heart model class. """

import json
import os
import tempfile

import yaml

if os.getenv("GITHUB_ACTIONS"):
    is_gh_action = True
else:
    is_gh_action = False

import ansys.heart.preprocessor.models as models
import numpy as np
import pytest


def _get_test_model_info() -> models.ModelInfo:
    """Get a test model info and populates it."""
    info = models.ModelInfo(
        work_directory=None,
        path_to_simulation_mesh="path-to-simulation-mesh",
        mesh_size=2.0,
        part_definitions={
            "Left ventricle": {"id": 1, "enclosed_by_boundaries": {"endocardium": 2}}
        },
    )
    models.ModelInfo()
    return info


@pytest.mark.parametrize("extension", [".json", ".yml"])
def test_model_info_dump(extension):
    """Test dumping of model info to json."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:
        info = _get_test_model_info()
        info.workdir = workdir

        path_to_model_info = os.path.join(workdir, "model_info" + extension)

        info.dump_info(path_to_model_info)

        path_to_model_info
        if extension == ".json":
            with open(path_to_model_info) as json_file:
                data = json.load(json_file)
        else:
            with open(path_to_model_info) as json_file:
                data = yaml.load(json_file, yaml.SafeLoader)

        assert info.workdir == data["workdir"], "Workdir not the same"
        assert (
            info.path_to_simulation_mesh == data["path_to_simulation_mesh"]
        ), "Path to simulation mesh not the same"
        assert info.mesh_size == data["mesh_size"], "Mesh size not the same"
        assert info.part_definitions == data["part_definitions"]
        assert info.path_to_model == data["path_to_model"]

    pass


def test_dump_model_001():
    """Test dumping of model to disk: using path in ModelInfo."""
    from pathlib import Path

    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:
        info = models.ModelInfo(work_directory=workdir)
        model = models.BiVentricle(info)

        expected_path = os.path.join(model.info.workdir, "heart_model.pickle")

        model.dump_model()
        assert os.path.isfile(expected_path)
        assert model.info.path_to_model == expected_path

        expected_path = os.path.join(model.info.workdir, "heart_model1.pickle")
        model.dump_model(expected_path)
        assert os.path.isfile(expected_path)
        assert model.info.path_to_model == expected_path

        expected_path = Path(os.path.join(model.info.workdir, "heart_model2.pickle"))
        model.dump_model(expected_path)
        assert os.path.isfile(expected_path)
        assert model.info.path_to_model == str(expected_path)


def test_model_load_001():
    """Test dumping and reading of model with data."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:

        info = models.ModelInfo()
        info.workdir = workdir
        model = models.BiVentricle(info)

        model.info.path_to_model = os.path.join(model.info.workdir, "heart_model.pickle")
        model.left_ventricle.endocardium.triangles = np.array([[0, 1, 2]], dtype=int)
        model.left_ventricle.endocardium.nodes = np.eye(3, 3, dtype=float)

        model.dump_model()

        assert os.path.isfile(model.info.path_to_model)

        model_loaded: models.BiVentricle = models.HeartModel.load_model(model.info.path_to_model)
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
        info = _get_test_model_info()
        model: models.BiVentricle = models.BiVentricle(info)
        model.info.workdir = workdir
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
    info = models.ModelInfo()
    model: models.HeartModel = model_type(info)

    assert model.part_names == expected_part_names
