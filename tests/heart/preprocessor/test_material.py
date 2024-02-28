from ansys.heart.preprocessor.material.curve import ActiveCurve, Strocchi_active, constant_ca2
import numpy as np

# import ansys.heart.preprocessor.material.material as M
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
