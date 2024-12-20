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
import shutil

import pytest

from ansys.heart.core.models import LeftVentricle
from ansys.heart.postprocessor.aha17_strain import AhaStrainCalculator
from ansys.heart.postprocessor.auto_process import zerop_post
import ansys.heart.postprocessor.dpf_utils as dpf_utils
from ansys.heart.postprocessor.system_model_post import SystemModelPost
from tests.heart.conftest import get_assets_folder

# Kills ansyscl. May be necessary to pass tests locally
dpf_utils._KILL_ANSYSCL_ON_DEL = True
# Accept DPF LA
os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"


@pytest.fixture(autouse=True, scope="module")
def get_data():
    test_dir = os.path.join(get_assets_folder(), "post")
    path_to_model = os.path.join(test_dir, "model", "heart_model.vtu")
    model: LeftVentricle = LeftVentricle(working_directory=test_dir)
    model.load_model_from_mesh(path_to_model, path_to_model.replace(".vtu", ".partinfo.json"))

    return test_dir, model


def test_compute_myocardial_strain(get_data):
    test_dir = get_data[0]
    model = get_data[1]
    d3plot = os.path.join(os.path.join(test_dir, "main", "d3plot"))

    s = AhaStrainCalculator(model, d3plot)
    _, aha_lrc, _ = s._compute_myocardial_strain(1)
    assert aha_lrc[-1, -1] == pytest.approx(0.0829107005807934)


def test_compute_aha_strain(get_data):
    test_dir = get_data[0]
    model = get_data[1]
    d3plot = os.path.join(os.path.join(test_dir, "main", "d3plot"))

    s = AhaStrainCalculator(model, d3plot)
    aha_lrc = s.compute_aha_strain()

    assert aha_lrc[1, -1] == pytest.approx(0.0829107005807934)


def test_zerop_post(get_data):
    test_dir = get_data[0]
    model = get_data[1]
    dct = zerop_post(os.path.join(test_dir, "zerop"), model)
    assert dct[0]["True left ventricle volume (mm3)"] == pytest.approx(118078.82768066938)

    # Cleanup
    folder = os.path.join(test_dir, "zerop", "post")
    if os.path.exists(folder):
        shutil.rmtree(folder)


class TestSystemModelPost:
    @pytest.fixture
    def system_model(self, get_data):
        test_dir = get_data[0]
        return SystemModelPost(os.path.join(test_dir, "main"))

    def test_plot_pv_loop(self, system_model):
        ef = system_model.get_ejection_fraction()
        fig = system_model.plot_pv_loop(ef=ef)
        fig.savefig("pv_{0:d}.png".format(0))

        assert os.path.isfile("pv_0.png")
        # clean
        os.remove("pv_0.png")
