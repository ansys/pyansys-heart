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

"""unit test for post-processing."""

import os
from pathlib import Path

import pytest

from ansys.heart.core.models import HeartModel
from ansys.heart.postprocessor.aha17_strain import AhaStrainCalculator
from ansys.heart.postprocessor.auto_process import mech_post, zerop_post
from ansys.heart.postprocessor.system_model_post import SystemModelPost

model: HeartModel
test_dir: str


pytestmark = pytest.mark.local


# TODO: use mock objects to allow testing postprocessing methods.
@pytest.fixture(autouse=True, scope="module")
def get_data():
    global test_dir, model

    test_dir = r"D:\PyAnsys-Heart\test_case\test_lv"
    model = HeartModel.load_model(Path(test_dir) / "model_with_fiber.pickle")
    model.compute_left_ventricle_anatomy_axis()
    model.compute_left_ventricle_aha17()


@pytest.mark.xfail(reason="Test requires local data.")
def test_compute_myocardial_strain():
    d3plot = Path(test_dir) / "main-mechanics" / "d3plot"

    s = AhaStrainCalculator(model, d3plot)
    _, aha_lrc, _ = s._compute_myocardial_strain(1)
    assert aha_lrc[-1, -1] == pytest.approx(0.08878163)


@pytest.mark.xfail(reason="Test requires local data.")
def test_compute_aha_strain():
    d3plot = Path(test_dir) / "main-mechanics" / "d3plot"

    s = AhaStrainCalculator(model, d3plot)
    aha_lrc = s.compute_aha_strain(".")
    assert aha_lrc[1, -1] == pytest.approx(0.08878163)


@pytest.mark.xfail(reason="Test requires local data.")
def test_mech_post():
    mech_post(Path(test_dir) / "main-mechanics", model)
    assert os.path.exists(Path(test_dir) / "main-mechanics" / "post")


@pytest.mark.xfail(reason="Test requires local data.")
def test_zerop_post():
    dct = zerop_post(Path(test_dir) / "zeropressure", model)
    assert dct["True left ventricle volume (mm3)"] == pytest.approx(288876.8)
    assert os.path.exists(Path(test_dir) / "zeropressure" / "post")


class TestSystemModelPost:
    @pytest.fixture
    def system_model(self):
        return SystemModelPost(Path(test_dir) / "main-mechanics")

    @pytest.mark.xfail(reason="Test requires local data.")
    def test_plot_pv_loop(self, system_model):
        ef = system_model.get_ejection_fraction(t_start=2, t_end=3)
        fig = system_model.plot_pv_loop(ef=ef)
        fig.savefig("pv_{0:d}.png".format(0))
        assert os.path.isfile("pv_0.png")
        # ffmpeg -f image2 -i pv_%d.png output.mp4
