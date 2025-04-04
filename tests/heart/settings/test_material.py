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

import numpy as np
import pytest

from ansys.heart.core.settings.material.curve import ActiveCurve, constant_ca2, strocchi_active
import ansys.heart.core.settings.material.material as material


class TestCa2Curve:
    @pytest.fixture
    def ca2_curve(self):
        return ActiveCurve(constant_ca2(), threshold=0.1, type="ca2")

    @pytest.fixture
    def stress_curve(self):
        return ActiveCurve(strocchi_active(), threshold=0.1, type="stress")

    def test_stress_to_ca2(self, stress_curve):
        t, v = stress_curve.dyna_input

        # test always activated after t=0
        assert v[0] < stress_curve.threshold
        assert np.all(v[1:] > stress_curve.threshold)

    def test_repeat(self, ca2_curve):
        t, v = ca2_curve.dyna_input
        assert len(t) == 501

        ca2_curve.n_beat = 10
        t, v = ca2_curve.dyna_input
        assert len(t) == 1001

    def test_check_threshold(self):
        # threshold larger than max
        with pytest.raises(ValueError) as e:
            ActiveCurve(constant_ca2(), threshold=4.36, type="ca2")
            assert str(e.value) == "Threshold must cross ca2+ curve at least once"

        # threshold lower than min
        with pytest.raises(ValueError) as e:
            ActiveCurve(constant_ca2(), threshold=-0.1, type="ca2")
            assert str(e.value) == "Threshold must cross ca2+ curve at least once"


def test_iso():
    iso = material.ISO(k1=1, k2=2)
    assert iso.itype == -3

    iso2 = material.ISO(mu1=1, alpha1=2, itype=-1)
    assert iso2.itype == -1
    assert iso2.nu == 0.499

    iso3 = material.ISO(mu1=1, alpha1=2, itype=-1, kappa=100)
    assert iso3.nu == pytest.approx(0.495, abs=0.01)


def test_aniso():
    fiber = material.ANISO.HGOFiber(k1=1, k2=2)
    sheet = material.ANISO.HGOFiber(k1=1, k2=2)
    aniso = material.ANISO(fibers=[fiber, sheet])
    assert aniso.nf == 2
    assert aniso.intype == 0
    assert fiber._theta == 0.0
    assert sheet._theta == 90.0

    aniso = material.ANISO(fibers=[fiber, sheet], k1fs=1, k2fs=2)
    assert aniso.intype == 1


def test_active():
    a = material.ACTIVE(ca2_curve=ActiveCurve(strocchi_active(), type="ca2", threshold=0.1))
    assert a.actype == 1
    assert a.acthr == 0.1
    assert a.model.l == 1.85
    assert a.ca2_curve.type == "ca2"


def test_active_couple():
    active_model = material.ActiveModel.Model3(
        ca2ion50=0.001,
        n=2,
        f=0.0,
        l=1.9,
        eta=1.45,
        sigmax=0.125,  # MPa
    )
    assert active_model.f == 0.0
    active = material.ACTIVE(
        sf=1.0, ss=0.0, sn=0.0, acthr=0.0002, model=active_model, ca2_curve=None
    )
    assert active.acthr == 0.0002
    assert active.acdir == 1


def test_mat295():
    m = material.Mat295(rho=1, iso=material.ISO(k1=1, k2=1))
    assert m.aniso is None
    assert m.active is None

    m.aniso = material.ANISO()
    assert m.aniso is not None


def test_neohookean():
    m = material.NeoHookean(rho=1, c10=1, nu=0.499)
    assert m.nu == 0.499


def test_dummy(capsys):
    m = material.MechanicalMaterialModel.DummyMaterial()
    print(m)
    captured = capsys.readouterr()
    assert captured.out == "Material is empty.\n"
