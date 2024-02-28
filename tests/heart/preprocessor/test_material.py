from ansys.heart.preprocessor.material.curve import ActiveCurve, Strocchi_active, constant_ca2
import ansys.heart.preprocessor.material.material as M
import numpy as np
import pytest


class TestCa2Curve:
    @pytest.fixture
    def ca2_curve(self):
        return ActiveCurve(constant_ca2(), threshold=0.1, type="ca2")

    @pytest.fixture
    def stress_curve(self):
        return ActiveCurve(Strocchi_active(), threshold=0.1, type="stress")

    def test_stress_to_ca2(self, stress_curve):
        t, v = stress_curve.dyna_input

        # test always activated after t=0
        assert v[0] < stress_curve.threshold
        assert np.all(v[1:] > stress_curve.threshold)

    def test_repeat(sellf, ca2_curve):
        t, v = ca2_curve.dyna_input
        assert len(t) == 501

        ca2_curve.n_beat = 10
        t, v = ca2_curve.dyna_input
        assert len(t) == 1001

    def test_check_threshold(self):
        # threshold larger than max
        with pytest.raises(ValueError) as e:
            bad = ActiveCurve(constant_ca2(), threshold=4.36, type="ca2")
            assert str(e.value) == "Threshold must cross ca2+ curve at least once"

        # threshold lower than min
        with pytest.raises(ValueError) as e:
            bad = ActiveCurve(constant_ca2(), threshold=-0.1, type="ca2")
            assert str(e.value) == "Threshold must cross ca2+ curve at least once"


def test_iso():
    iso = M.ISO(k1=1, k2=2)
    assert iso.itype == -3


def test_aniso():
    fiber = M.ANISO.HGO_Fiber(k1=1, k2=2)
    sheet = M.ANISO.HGO_Fiber(k1=1, k2=2)
    aniso = M.ANISO(fibers=[fiber, sheet])
    assert aniso.nf == 2
    assert aniso.intype == 0
    assert fiber._theta == 0.0
    assert sheet._theta == 90.0

    aniso = M.ANISO(fibers=[fiber, sheet], k1fs=1, k2fs=2)
    assert aniso.intype == 1


def test_active():
    a = M.ACTIVE(ca2_curve=ActiveCurve(Strocchi_active(), type="ca2", threshold=0.1))
    assert a.actype == 1
    assert a.acthr == 0.1


def test_mat295():
    m = M.MAT295(rho=1, iso=M.ISO())
    assert m.aniso == None
    assert m.active == None

    m.aniso = M.ANISO()
    assert m.aniso != None


def test_neohookean():
    m = M.NeoHookean(rho=1, c10=1)
    assert m.nu == 0.499


def test_dummy(capsys):
    m = M.MechaMaterialModel.DummyMaterial()
    print(m)
    captured = capsys.readouterr()
    assert captured.out == "Material is empty.\n"
