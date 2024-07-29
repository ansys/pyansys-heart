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

if os.getenv("GITHUB_ACTIONS"):
    is_gh_action = True
else:
    is_gh_action = False

import numpy as np
import pytest

import ansys.heart.preprocessor.models as models
import ansys.heart.writer.dynawriter as writers


@pytest.fixture
def model_info() -> models.ModelInfo:
    """Get a test model info and populates it."""
    return models.ModelInfo(
        work_directory=None,
        path_to_simulation_mesh="path-to-simulation-mesh",
        mesh_size=2.0,
        part_definitions={
            "Left ventricle": {"id": 1, "enclosed_by_boundaries": {"endocardium": 2}}
        },
    )


def get_test_model():
    model: models.LeftVentricle = models.LeftVentricle(model_info)

    # populate model
    model.mesh.tetrahedrons = np.array([[0, 1, 2, 3]], dtype=int)
    model.mesh.nodes = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=float
    )
    model.left_ventricle.element_ids = np.array([0], dtype=int)

    model.left_ventricle.endocardium.triangles = np.array([[0, 1, 2]], dtype=int)
    model.left_ventricle.endocardium.points = model.mesh.nodes

    return model


def test_filter_bc_nodes01():
    # filter is not necessary
    model = get_test_model()
    model.left_ventricle.endocardium.triangles = np.array([[0, 1, 2]], dtype=int)
    w = writers.FiberGenerationDynaWriter(model)
    ids = w._filter_bc_nodes(model.left_ventricle.endocardium)

    # no node is removed
    assert np.allclose(ids, [0, 1, 2])


def test_filter_bc_nodes02():
    # filter is necessary
    model = get_test_model()
    model.left_ventricle.endocardium.triangles = np.array([[0, 1, 2], [0, 1, 3]], dtype=int)
    w = writers.FiberGenerationDynaWriter(model)
    ids = w._filter_bc_nodes(model.left_ventricle.endocardium)

    # one node should be removed
    # TODO @oliver, which node should be removed by your logic?
    assert np.allclose(ids, [1, 2, 3])
