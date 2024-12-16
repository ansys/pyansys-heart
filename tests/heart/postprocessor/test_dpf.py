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

os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"

import numpy as np
import pytest

from ansys.dpf import core as dpf
from ansys.heart.postprocessor.dpf_utils import D3plotReader, ICVoutReader


def test_icvout_reader():
    path = os.path.dirname(os.path.abspath(__file__))
    fn = os.path.join(path, "asset", "main", "binout")
    icvout = ICVoutReader(fn)

    assert np.all(icvout._icv_ids == np.array([1]))
    assert np.all(icvout._icvi_ids == np.array([1]))

    assert icvout.get_time()[1] == pytest.approx(5.0, rel=1e-3)
    assert icvout.get_time()[1] == pytest.approx(5.0, rel=1e-3)
    assert icvout.get_pressure(1)[-1] == pytest.approx(0.001700745546258986, rel=1e-3)
    assert icvout.get_volume(1)[-1] == pytest.approx(109582.515625, rel=1e-3)
    assert icvout.get_flowrate(1)[-1] == pytest.approx(-64.24286651611328, rel=1e-1)

    with pytest.raises(ValueError):
        icvout.get_flowrate(3)


def test_d3plot_reader():
    path = os.path.dirname(os.path.abspath(__file__))
    fn = os.path.join(path, "asset", "main", "d3plot")
    d3plot = D3plotReader(fn)

    # just to check all dpf API works
    assert len(d3plot.time) == 2
    assert d3plot.get_material_ids().max() == 8
    assert d3plot.get_initial_coordinates().shape == (8598, 3)
    assert d3plot.get_history_variable([1, 2], at_step=0).shape == (2, 24206)
    assert isinstance(d3plot.get_ep_fields(), dpf.FieldsContainer)


@pytest.xfail(reason="This function fails on Github but not local.")
def test_d3plot_reader2():
    path = os.path.dirname(os.path.abspath(__file__))
    fn = os.path.join(path, "asset", "main", "d3plot")
    d3plot = D3plotReader(fn)

    assert d3plot.get_displacement_at(0.0).shape == (8598, 3)
