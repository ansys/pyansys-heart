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

import copy

import numpy as np
import pyvista as pv

import ansys.health.heart.models as models
from ansys.health.heart.objects import Mesh, SurfaceMesh
import ansys.health.heart.writer.dynawriter as writers


def get_test_model():
    model: models.LeftVentricle = models.LeftVentricle()

    # populate model
    cells = np.array([4, 1, 2, 3, 4, 4, 5, 2, 3, 4], dtype=int)
    points = np.array(
        [
            [9.1, 9.1, 9.1],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ],
        dtype=float,
    )
    celltypes = [pv.CellType.TETRA] * 2

    model.mesh = Mesh(cells, celltypes, points)
    model.left_ventricle.element_ids = np.array([0, 1], dtype=int)

    model.left_ventricle.endocardium.triangles = np.array([[1, 2, 3]], dtype=int)
    model.left_ventricle.endocardium.points = model.mesh.points
    model.left_ventricle.endocardium.point_data["_global-point-ids"] = np.arange(
        0, model.mesh.n_points
    )

    return model


def test_filter_bc_nodes01():
    model = copy.deepcopy(get_test_model())
    # filter is not necessary
    model.left_ventricle.endocardium.triangles = np.array([[1, 2, 3]], dtype=int)
    model.left_ventricle.endocardium = SurfaceMesh(model.left_ventricle.endocardium.clean())

    w = writers.FiberGenerationDynaWriter(model)
    ids = w._filter_bc_nodes(model.left_ventricle.endocardium)

    # no node is removed
    assert np.allclose(ids, [1, 2, 3])


def test_filter_bc_nodes02():
    model = copy.deepcopy(get_test_model())
    # filter is necessary
    model.left_ventricle.endocardium.triangles = np.array(
        [[1, 2, 3], [1, 2, 4], [1, 3, 4]], dtype=int
    )
    model.left_ventricle.endocardium = SurfaceMesh(model.left_ventricle.endocardium.clean())

    w = writers.FiberGenerationDynaWriter(model)
    ids = w._filter_bc_nodes(model.left_ventricle.endocardium)

    # one node is removed, but not the node 0 since it's only attached with the issue element
    # see #656
    assert np.allclose(ids, [1, 3, 4])


def test_filter_bc_nodes03():
    model = copy.deepcopy(get_test_model())
    # filter is necessary but no individual nodes can be removed
    model.left_ventricle.endocardium.triangles = np.array(
        [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 5]], dtype=int
    )
    model.left_ventricle.endocardium = SurfaceMesh(model.left_ventricle.endocardium.clean())

    w = writers.FiberGenerationDynaWriter(model)
    ids = w._filter_bc_nodes(model.left_ventricle.endocardium)

    # all nodes in the unsolvable tetrahedron and their respective neighbors are removed
    # see #656
    assert len(ids) == 0
