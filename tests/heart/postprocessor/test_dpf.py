# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
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

"""unit test for dpf utils."""
import os

from ansys.heart.postprocessor.dpf_utils import ICVoutReader
import numpy as np
import pytest


def test_icvout():
    path = os.path.dirname(os.path.abspath(__file__))
    fn = os.path.join(path, "binout0000")
    icvout = ICVoutReader(fn)
    assert np.all(icvout._icv_ids == np.array([1, 2]))
    assert np.all(icvout._icvi_ids == np.array([1, 2]))

    assert icvout.get_pressure(1)[-1] == pytest.approx(0.0019017, rel=1e-3)
    assert icvout.get_volume(1)[-1] == pytest.approx(166579.15625, rel=1e-3)
    assert icvout.get_flowrate(1)[-1] == pytest.approx(-105.39063, rel=1e-3)

    with pytest.raises(ValueError):
        icvout.get_flowrate(3)
