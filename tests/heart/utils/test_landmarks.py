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

import numpy as np
import pytest

import ansys.health.heart.models as models
from ansys.health.heart.utils.landmark_utils import (
    compute_aha17,
    compute_anatomy_axis,
    compute_element_cs,
)
from tests.heart.conftest import get_assets_folder


@pytest.fixture
def model() -> models.FourChamber:
    vtu_file = os.path.join(
        get_assets_folder(),
        "reference_models",
        "strocchi2020",
        "01",
        "FourChamber",
        "heart_model.vtu",
    )

    json_file = os.path.join(
        get_assets_folder(),
        "reference_models",
        "strocchi2020",
        "01",
        "FourChamber",
        "heart_model.partinfo.json",
    )

    model: models.FourChamber = models.FourChamber(working_directory=".")

    model.load_model_from_mesh(vtu_file, json_file)

    return model


def test_compute_anatomy_axis():
    mv_center = np.array([14.91736773, 138.79881432, 382.26438557])
    av_center = np.array([4.56652835, 118.0155877, 391.6282463])
    apex = np.array([70.69582127, 71.67369568, 353.56092681])

    a, b, c = compute_anatomy_axis(mv_center, av_center, apex)
    assert np.allclose(b["center"], np.array([42.8065945, 105.236255, 367.91265619]))


def test_compute_aha17(model):
    l4cv = {
        "center": np.array([30.05990578, 109.49603257, 375.8178529]),
        "normal": np.array([-0.54847906, -0.10082087, -0.83006378]),
    }
    short = {
        "center": np.array([26.07305844, 125.37379059, 376.52369382]),
        "normal": np.array([0.60711636, -0.73061828, -0.31242063]),
    }

    aha_ids = compute_aha17(model, short, l4cv)
    # solid = model.mesh.extract_cells_by_type(10)
    # solid.cell_data["aha"] = aha_ids
    # solid.save("aha.vtk")

    assert np.sum(aha_ids == 10) == 10591
    assert np.sum(aha_ids == 12) == 11220
    assert np.sum(aha_ids == 17) == 2526
    assert np.sum(np.isnan(aha_ids)) == 361818


def test_compute_element_cs(model):
    l4cv = {
        "center": np.array([30.05990578, 109.49603257, 375.8178529]),
        "normal": np.array([-0.54847906, -0.10082087, -0.83006378]),
    }
    short = {
        "center": np.array([26.07305844, 125.37379059, 376.52369382]),
        "normal": np.array([0.60711636, -0.73061828, -0.31242063]),
    }
    aha_ids = compute_aha17(model, short, l4cv)
    aha_elems = np.where(~np.isnan(aha_ids))[0]
    longitudinal, radial, circ = compute_element_cs(model, short, aha_elems)

    assert longitudinal.shape == (174298, 3)
