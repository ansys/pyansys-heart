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

from matplotlib.figure import Figure
import numpy as np
import pytest

from ansys.heart.core.postprocessor.klotz_curve import EDPVR


def test_edpvr():
    v = 120  # mL
    p = 15  # mmHg

    klotz = EDPVR(v, p)

    assert klotz.v0 == pytest.approx(61.2, 0.01)

    assert klotz.get_pressure(100) == pytest.approx(5.0287, 0.001)
    assert np.allclose(klotz.get_pressure(np.array([100, 120])), np.array([5.0287, 15]))

    assert np.allclose(klotz.get_volume(np.array([15])), np.array([120.0]))
    assert np.allclose(klotz.get_volume(np.array([0])), np.array([61.2]))

    fig = klotz.plot_EDPVR()
    assert isinstance(fig, Figure)
