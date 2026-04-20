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

import numpy as np
import pytest

from ansys.health.heart.settings.material.material import Mat295
from ansys.health.heart.settings.material.material_factory import (
    _get_myocardium_material,
    _get_passive_material,
)


class _Q:
    """Minimal stand-in for Quantity-like objects used by factory (.m attribute)."""

    def __init__(self, value: float) -> None:
        self.m = value


def _mock_settings(
    rho: float = 1000.0,
    kappa: float = 1.0,
    k1: float = 0.1,
    k2: float = 0.2,
    k1f: float = 2.0,
    k2f: float = 1.0,
    taumax: float = 150.0,
    beat_time: float = 800.0,
    ss: float = 1.0,
    sn: float = 1.0,
) -> dict[str, dict[str, _Q | float]]:
    return {
        "isotropic": {
            "rho": _Q(rho),
            "kappa": _Q(kappa),
            "k1": _Q(k1),
            "k2": _Q(k2),
        },
        "anisotropic": {"k1f": _Q(k1f), "k2f": _Q(k2f)},
        "active": {"taumax": _Q(taumax), "beat_time": _Q(beat_time), "ss": ss, "sn": sn},
    }


def test_get_myocardium_material_returns_mat295_with_active() -> None:
    """Myocardium factory builds Mat295 with iso/aniso and active when ep_coupled=False."""
    settings = _mock_settings()
    mat = _get_myocardium_material(settings, ep_coupled=False)

    assert isinstance(mat, Mat295)
    # isotropic values
    assert pytest.approx(mat.rho, rel=1e-6) == 1000.0
    assert pytest.approx(mat.iso.kappa, rel=1e-6) == 1.0
    assert pytest.approx(mat.iso.k1, rel=1e-6) == 0.1
    assert pytest.approx(mat.iso.k2, rel=1e-6) == 0.2

    # anisotropic fibers
    assert mat.aniso is not None
    assert len(mat.aniso.fibers) == 1
    assert pytest.approx(mat.aniso.fibers[0].k1, rel=1e-6) == 2.0
    assert pytest.approx(mat.aniso.fibers[0].k2, rel=1e-6) == 1.0

    # active component present when not ep_coupled
    assert mat.active is not None
    assert pytest.approx(mat.active.model.taumax, rel=1e-6) == 150.0
    # ca2_curve should be present (ActiveCurve -> has func or similar)
    assert hasattr(mat.active, "ca2_curve")

    # additional assertions for active parameters and curve shape
    assert mat.active.ss == 1.0
    assert mat.active.sn == 1.0
    t, y = mat.active.ca2_curve.dyna_input
    assert isinstance(t, np.ndarray)
    assert isinstance(y, np.ndarray)
    assert t.shape == y.shape
    assert t.size > 1
    assert np.all(np.diff(t) > 0)


def test_get_myocardium_material_ep_coupled_disables_active() -> None:
    """Verify ep_coupled=True produces a Mat295 without an active component."""
    settings = _mock_settings()
    mat = _get_myocardium_material(settings, ep_coupled=True)

    assert isinstance(mat, Mat295)
    # active should be None when ep_coupled is True
    assert mat.active.ca2_curve is None


def test_get_myocardium_material_missing_keys_raises() -> None:
    """Factory raises when isotropic/anisotropic sections are missing."""
    # missing 'isotropic' section entirely
    bad_settings = {"anisotropic": {}, "active": {}}
    with pytest.raises(Exception):
        _get_myocardium_material(bad_settings)


# passive tests


def _mock_passive_settings(
    rho: float = 1050.0, itype: int = -1, kappa: float = 0.5, mu1: float = 0.3, alpha1: float = 0.2
) -> dict[str, _Q | int | float]:
    """Create minimal passive settings structure expected by the factory."""
    return {
        "rho": _Q(rho),
        "itype": itype,
        "kappa": _Q(kappa),
        "mu1": _Q(mu1),
        "alpha1": alpha1,
    }


def test_get_passive_material_returns_mat295() -> None:
    """Passive factory returns Mat295 with expected isotropic parameters."""
    settings = _mock_passive_settings()
    mat = _get_passive_material(settings)

    assert isinstance(mat, Mat295)
    # rho taken from Quantity-like .m
    assert pytest.approx(mat.rho, rel=1e-6) == 1050.0
    # iso fields set as expected
    assert mat.iso.itype == -1
    assert pytest.approx(mat.iso.kappa, rel=1e-6) == 0.5
    assert pytest.approx(mat.iso.mu1, rel=1e-6) == 0.3
    # alpha1 in factory is passed directly
    assert pytest.approx(mat.iso.alpha1, rel=1e-6) == 0.2
    # beta is set in factory to 2
    assert mat.iso.beta == 2


def test_get_passive_material_missing_key_raises() -> None:
    """Passive factory raises KeyError when required passive keys are absent."""
    bad = {"rho": _Q(1000.0), "itype": -3}  # missing kappa/mu1/alpha1
    with pytest.raises(KeyError):
        _get_passive_material(bad)
