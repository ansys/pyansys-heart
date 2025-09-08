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

from ansys.health.heart.post.laplace_post import (
    compute_la_fiber_cs,
    compute_ra_fiber_cs,
    compute_rotation_angle1,
    compute_ventricle_fiber_by_drbm,
    orthogonalization2,
)
from ansys.health.heart.settings.settings import AtrialFiber
from tests.heart.conftest import get_assets_folder


@pytest.fixture(autouse=True)
def _set_env_vars(monkeypatch):
    """Set environment variables for testing."""
    monkeypatch.setenv("ANSYS_DPF_ACCEPT_LA", "Y")


def test_compute_la_fiber_cs(_set_env_vars):
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
        "ansys.health.heart.post.laplace_post.read_laplace_solution", return_value=input_grid
    ):
        res = compute_la_fiber_cs(".", setting, la_endo)
        assert np.sum(res["bundle"] == 0) == 33097
        assert np.allclose(res["e_l"][0], np.array([0.54939723, 0.77379318, -0.31528842]))

        assert pytest.approx(np.dot(res["e_l"][0], res["e_n"][0])) == 0
        assert pytest.approx(np.dot(res["e_l"][0], res["e_t"][0])) == 0


def test_compute_ra_fiber_cs(_set_env_vars):
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
        "ansys.health.heart.post.laplace_post.read_laplace_solution", return_value=input_grid
    ):
        res = compute_ra_fiber_cs(".", setting, la_endo)

        assert np.sum(res["bundle"] == 0) == 0
        assert np.sum(res["bundle"] == 10) == 15548
        assert np.allclose(res["e_l"][0], np.array([-0.3007031, 0.93551427, -0.18544755]))

        assert pytest.approx(np.dot(res["e_l"][0], res["e_n"][0])) == 0
        assert pytest.approx(np.dot(res["e_l"][0], res["e_t"][0])) == 0


def test_compute_ventricle_fiber_by_drbm(_set_env_vars):
    dir = os.path.join(get_assets_folder(), "post", "drbm")
    input_grid = pv.read(os.path.join(dir, "data.vtu"))

    with mock.patch(
        "ansys.health.heart.post.laplace_post.read_laplace_solution",
        return_value=input_grid,
    ):
        # test no outflow tract
        res = compute_ventricle_fiber_by_drbm(".")
        assert np.sum(res["label"] == 1) == 96695

        assert res["label"][0] == 1
        assert res["label"][-1] == 2
        assert np.isclose(np.min(res["alpha"]), -60, rtol=1e-1)
        assert np.isclose(np.max(res["alpha"]), 90, rtol=1e-1)
        assert np.allclose(res["fiber"][0], np.array([-0.07979692, 0.43617534, -0.89631664]))
        assert np.allclose(res["fiber"][-1], np.array([0.33103764, 0.77190818, 0.54274473]))

        # test with outflow tract
        res = compute_ventricle_fiber_by_drbm(
            ".",
            settings={
                "alpha_left": [40, -40],
                "alpha_right": [80, -25],
                "alpha_ot": [90, 0],
                "beta_left": [-20, 20],
                "beta_right": [0, 20],
                "beta_ot": [0, 0],
            },
        )
        assert np.sum(res["label"] == 1) == 96695

        assert np.isclose(np.min(res["alpha"]), -40, rtol=1e-1)
        assert np.isclose(np.max(res["alpha"]), 80, rtol=1e-1)
        assert res["label"][0] == 1
        assert res["label"][-1] == 2

        assert np.allclose(res["fiber"][0], np.array([0.05671137, 0.61760931, -0.78443774]))
        assert np.allclose(res["fiber"][-1], np.array([-0.26596645, -0.08034016, 0.9606286]))


def test_orthogonalization2():
    """Test the orthogonalization function."""
    # Create two non-parallel vectors for 5 cells
    e_1 = np.tile(np.array([1, 0, 0]), (5, 1))
    e_2 = np.tile(np.array([1, 1, 0]), (5, 1))

    v1, v2, v3 = orthogonalization2(e_1, e_2)

    # Check shapes
    assert v1.shape == (5, 3)
    assert v2.shape == (5, 3)
    assert v3.shape == (5, 3)

    # Check normalization
    np.testing.assert_allclose(np.linalg.norm(v1, axis=1), 1)
    np.testing.assert_allclose(np.linalg.norm(v2, axis=1), 1)
    np.testing.assert_allclose(np.linalg.norm(v3, axis=1), 1)

    # Check orthogonality
    np.testing.assert_allclose(np.einsum("ij,ij->i", v1, v2), 0, atol=1e-7)
    np.testing.assert_allclose(np.einsum("ij,ij->i", v1, v3), 0, atol=1e-7)
    np.testing.assert_allclose(np.einsum("ij,ij->i", v2, v3), 0, atol=1e-7)

    # Check right-hand rule
    np.testing.assert_allclose(np.cross(v1, v2), v3, atol=1e-7)


def test_compute_rotation_angle1():
    transmural_distance = np.array([0, 0.5, 1])
    rotation = [60, -60]
    angles = compute_rotation_angle1(transmural_distance, rotation)

    # expected angles at transmural distances 0, 0.5, and 1
    expected = np.array([60, 0, -60])
    np.testing.assert_allclose(angles, expected)
