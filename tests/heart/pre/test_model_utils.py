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

from ansys.health.heart.landmarks import LandMarks
import ansys.health.heart.models_utils as heart_model_utils
from tests.heart.conftest import get_fourchamber


def test_compute_sa_node():
    fourchamber = get_fourchamber()
    landmarks = LandMarks()

    sa_node = heart_model_utils.define_sino_atrial_node(fourchamber, landmarks)

    assert np.allclose(sa_node.xyz, np.array([-48.80218005, 107.90170883, 423.33688959]))
    assert sa_node.node_id == 105021


def test_compute_av_node():
    fourchamber = get_fourchamber()
    landmarks = LandMarks()

    av_node = heart_model_utils.define_atrio_ventricular_node(fourchamber, landmarks)
    assert np.allclose(av_node.xyz, np.array([-10.16353107, 108.95410155, 371.9505145]))
    assert av_node.node_id == 100501


def test_compute_his_bif_node():
    fourchamber = get_fourchamber()
    landmarks = LandMarks()

    # First define AV node as it is required for the His bifurcation
    heart_model_utils.define_atrio_ventricular_node(fourchamber, landmarks)

    bif_node = heart_model_utils.define_his_bundle_bifurcation_node(
        fourchamber, landmarks, target_coord=np.array([-10.16353107, 108.95410155, 371.9505145])
    )
    assert np.allclose(bif_node.xyz, np.array([1.22510233, 110.31896126, 364.402475]))
    assert bif_node.node_id == 25326


def test_compute_his_end_node():
    fourchamber = get_fourchamber()
    landmarks = LandMarks()

    # Compute prerequisite landmarks
    heart_model_utils.define_atrio_ventricular_node(fourchamber, landmarks)
    heart_model_utils.define_his_bundle_bifurcation_node(fourchamber, landmarks)

    # test for left
    left = heart_model_utils.define_his_bundle_end_node(fourchamber, landmarks, side="left")

    assert np.allclose(left.xyz, np.array([4.15421613, 113.63743565, 369.27104019]))
    assert left.node_id == 49464

    # test for right
    right = heart_model_utils.define_his_bundle_end_node(fourchamber, landmarks, side="right")

    assert np.allclose(right.xyz, np.array([2.93215687, 106.09459183, 365.20590901]))
    assert right.node_id == 43585
