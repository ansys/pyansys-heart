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

"""Contains tests for conduction system."""

import numpy as np
from pyvista.examples import load_tetbeam

from ansys.health.heart.pre.conduction_beam import (
    ConductionSystem,
    _create_line,
    _refine_line,
)


def test_find_path():
    mesh = load_tetbeam()
    mesh.point_data["_global-point-ids"] = np.arange(0, mesh.n_points)
    path = ConductionSystem.find_path(mesh, np.array([0, 0, 0]), np.array([1, 0, 5]), True)

    # TODO: add reasonable assertion
    assert path[0].shape[0] > 0
    assert len(path[1]) > 0


def test_create_line():
    points = _create_line(np.array([0, 0, 0]), np.array([0, 0, 1]), 0.5)
    assert np.allclose(points, np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.5], [0.0, 0.0, 1.0]]))


def test_refine_line():
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    refined_points = _refine_line(points, 0.5)
    assert np.allclose(
        refined_points, np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0]])
    )
