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
import unittest.mock as mock

import numpy as np
import pytest
import pyvista as pv

from ansys.heart.core.post.laplace_post import (
    compute_la_fiber_cs,
    compute_ra_fiber_cs,
    compute_ventricle_fiber_by_drbm,
)
from ansys.heart.simulator.settings.settings import AtrialFiber
from tests.heart.conftest import get_assets_folder

os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"


def test_compute_la_fiber_cs():
    dir = os.path.join(get_assets_folder(), "post", "la_fiber")

    setting = AtrialFiber()
    setting.set_values(
        {
            "tau_mv": 0.65,
            "tau_lpv": 0.1,
            "tau_rpv": 0.65,
        }
    )
    input_grid = pv.read(os.path.join(dir, "la_input.vtu"))
    la_endo = pv.read(os.path.join(dir, "la_endo.vtk"))

    with mock.patch(
        "ansys.heart.core.post.laplace_post.read_laplace_solution", return_value=input_grid
    ):
        res = compute_la_fiber_cs(".", setting, la_endo)
        assert np.sum(res["bundle"] == 0) == 33097
        assert np.allclose(res["e_l"][0], np.array([0.54939723, 0.77379318, -0.31528842]))

        assert pytest.approx(np.dot(res["e_l"][0], res["e_n"][0])) == 0
        assert pytest.approx(np.dot(res["e_l"][0], res["e_t"][0])) == 0


def test_compute_ra_fiber_cs():
    dir = os.path.join(get_assets_folder(), "post", "ra_fiber")

    setting = AtrialFiber()
    setting.set_values(
        {
            "tau_tv": 0.9,
            "tau_raw": 0.55,
            "tau_ct_minus": -0.18,
            "tau_ct_plus": -0.1,
            "tau_icv": 0.9,
            "tau_scv": 0.1,
            "tau_ib": 0.135,
            "tau_ras": 0.35,
        }
    )
    input_grid = pv.read(os.path.join(dir, "ra_input.vtu"))
    la_endo = pv.read(os.path.join(dir, "ra_endo.vtk"))

    with mock.patch(
        "ansys.heart.core.post.laplace_post.read_laplace_solution", return_value=input_grid
    ):
        res = compute_ra_fiber_cs(".", setting, la_endo)

        assert np.sum(res["bundle"] == 0) == 0
        assert np.sum(res["bundle"] == 10) == 15548
        assert np.allclose(res["e_l"][0], np.array([-0.3007031, 0.93551427, -0.18544755]))

        assert pytest.approx(np.dot(res["e_l"][0], res["e_n"][0])) == 0
        assert pytest.approx(np.dot(res["e_l"][0], res["e_t"][0])) == 0


def test_compute_ventricle_fiber_by_drbm():
    dir = os.path.join(get_assets_folder(), "post", "drbm")
    input_grid = pv.read(os.path.join(dir, "data.vtu"))

    with mock.patch(
        "ansys.heart.core.post.laplace_post.read_laplace_solution",
        return_value=input_grid,
    ):
        # test no outflow tract
        res = compute_ventricle_fiber_by_drbm(".")
        assert np.sum(res["label"] == 1) == 86210

        assert res["label"][0] == 1
        assert np.allclose(res["fiber"][0], np.array([0.03515216, 0.964732, 0.26087638]))
        assert res["label"][-1] == 2
        assert np.allclose(res["fiber"][-1], np.array([0.70611833, 0.32320625, -0.63002748]))

        # test with outflow tract
        res = compute_ventricle_fiber_by_drbm(
            ".",
            settings={
                "alpha_left": [-60, 60],
                "alpha_right": [90, -25],
                "alpha_ot": [90, 0],
                "beta_left": [-20, 20],
                "beta_right": [0, 20],
                "beta_ot": [0, 0],
            },
        )
        assert np.sum(res["label"] == 1) == 86210

        assert res["label"][0] == 1
        assert np.allclose(res["fiber"][0], np.array([0.03520965, 0.96476314, 0.26075345]))
        assert res["label"][-1] == 2
        assert np.allclose(res["fiber"][-1], np.array([-0.70556093, 0.21479473, 0.67531252]))
