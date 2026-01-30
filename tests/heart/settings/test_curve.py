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

# tests/heart/settings/test_curve_active.py
from ansys.health.heart.settings.material.curve import (
    ActiveCurve,
    constant_ca2,
    strocchi_active,
)


def test_repeat_default_and_custom_n() -> None:
    curve = ActiveCurve(func=constant_ca2(), threshold=0.1, type="ca2")
    t, v = curve.dyna_input
    # default n_beat = 5 --> 101 + 4*100 = 501
    assert len(t) == 501

    curve.n_beat = 10
    t, v = curve.dyna_input
    # 101 + 9*100 = 1001
    assert len(t) == 1001


def test_check_threshold_raises() -> None:
    # threshold larger than max
    with pytest.raises(ValueError):
        ActiveCurve(func=constant_ca2(), threshold=4.36, type="ca2")

    # threshold lower than min
    with pytest.raises(ValueError):
        ActiveCurve(func=constant_ca2(), threshold=-0.1, type="ca2")


def test_stress_to_ca2_conversion_and_threshold_reset() -> None:
    # When type == "stress", threshold is reset to 0.5e-6 inside ActiveCurve
    curve = ActiveCurve(func=strocchi_active(), type="stress")
    assert curve.threshold == pytest.approx(0.5e-6)

    t, v = curve.dyna_input
    # first value should be below threshold, rest should be above threshold
    assert v[0] < curve.threshold
    assert np.all(v[1:] > curve.threshold)


def test_invalid_stress_values_raise() -> None:
    time = np.linspace(0, 100, 5)

    # negative stress value
    stress_neg = np.array([-0.1, 0.0, 0.2, 0.3, 0.4])
    with pytest.raises(ValueError):
        ActiveCurve(func=(time, stress_neg), type="stress")

    # stress value greater than 1
    stress_gt1 = np.array([0.0, 1.2, 0.2, 0.3, 0.4])
    with pytest.raises(ValueError):
        ActiveCurve(func=(time, stress_gt1), type="stress")


def test_repeat_time_offsets() -> None:
    base_time, _ = constant_ca2()
    tb = base_time[-1]
    curve = ActiveCurve(func=constant_ca2(), threshold=0.1, type="ca2")
    t, v = curve.dyna_input

    assert t[0] == pytest.approx(0.0)
    # index 101 corresponds to the first appended point of the second beat
    assert t[101] == pytest.approx(base_time[1] + tb)
