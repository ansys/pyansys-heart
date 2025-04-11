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

"""Test for reading purkinje network as a beam mesh. Uses a mock Mesh object."""

import numpy as np
import pytest
import pyvista as pv

from ansys.health.heart.models_utils import HeartModelUtils
from ansys.health.heart.objects import _ConductionType
from ansys.health.heart.pre.conduction_beam import ConductionSystem
from ansys.health.heart.pre.conduction_beam2 import ConductionBeams, ConductionBeamType
from tests.heart.conftest import get_fourchamber


def test_compute_sa_node2():
    fourchamber = get_fourchamber()

    sa_node = HeartModelUtils.define_sino_atrial_node(fourchamber)

    assert np.allclose(sa_node.xyz, np.array([-48.80218005, 107.90170883, 423.33688959]))
    assert sa_node.node_id == 105021


def test_compute_av_node2():
    fourchamber = get_fourchamber()
    av_node = HeartModelUtils.define_atrio_ventricular_node(fourchamber)
    assert np.allclose(av_node.xyz, np.array([-10.16353107, 108.95410155, 371.9505145]))
    assert av_node.node_id == 100501


def test_compute_his_bif_node():
    fourchamber = get_fourchamber()
    bif_node = HeartModelUtils.define_his_bundle_bifurcation_node(
        fourchamber, target_coord=np.array([-10.16353107, 108.95410155, 371.9505145])
    )
    assert np.allclose(bif_node.xyz, np.array([1.22510233, 110.31896126, 364.402475]))
    assert bif_node.node_id == 25326


def test_compute_his_end_node():
    fourchamber = get_fourchamber()
    # need pre steps
    HeartModelUtils.define_atrio_ventricular_node(fourchamber)
    HeartModelUtils.define_his_bundle_bifurcation_node(fourchamber)
    # test for left
    left = HeartModelUtils.define_his_bundle_end_node(fourchamber, side="left")

    assert np.allclose(left.xyz, np.array([4.15421613, 113.63743565, 369.27104019]))
    assert left.node_id == 49464
    # test for right
    right = HeartModelUtils.define_his_bundle_end_node(fourchamber, side="right")

    assert np.allclose(right.xyz, np.array([2.93215687, 106.09459183, 365.20590901]))
    assert right.node_id == 43585


def test_create_conductionbeams_on_surface():
    model = get_fourchamber()
    sa = HeartModelUtils.define_sino_atrial_node(model)
    av = HeartModelUtils.define_atrio_ventricular_node(model)

    sa_av = ConductionBeams.create_from_keypoints(
        name=ConductionBeamType.SAN_AVN,
        keypoints=[sa.xyz, av.xyz],
        id=2,
        base_mesh=model.right_atrium.endocardium,
        connection="none",
    )

    assert sa_av.mesh.n_lines == 48
    assert np.isclose(sa_av.length, 64.36592438345)


def test_create_conductionbeams_in_solid():
    model = get_fourchamber()
    av = HeartModelUtils.define_atrio_ventricular_node(model)
    bif = HeartModelUtils.define_his_bundle_bifurcation_node(model)
    his_top = ConductionBeams.create_from_keypoints(
        name=ConductionBeamType.HIS_TOP,
        keypoints=[av.xyz, bif.xyz],
        id=1,
        base_mesh=model.mesh,
        connection="none",
    )
    assert np.isclose(his_top.length, 14.276232139149878)
    assert his_top.relying_surface.n_cells == 9


def _mock_purkinje():
    """Create a mock Purkinje mesh."""
    pv.PolyData()


def test_compute_sa_node():
    fourchamber = get_fourchamber()
    cs = ConductionSystem(fourchamber)
    sa_node = cs.compute_sa_node()

    assert np.allclose(sa_node.xyz, np.array([-48.80218005, 107.90170883, 423.33688959]))
    assert sa_node.node_id == 105021


def test_compute_av_node():
    fourchamber = get_fourchamber()
    cs = ConductionSystem(fourchamber)
    av_node = cs.compute_av_node()
    assert np.allclose(av_node.xyz, np.array([-10.16353107, 108.95410155, 371.9505145]))
    assert av_node.node_id == 100501


def test_compute_av_conduction():
    # get a fresh model
    fourchamber = get_fourchamber()
    cs = ConductionSystem(fourchamber)

    cs.compute_sa_node()
    cs.compute_av_node()
    beam = cs.compute_av_conduction()

    # Following assert based on type hint, but will fail:
    # assert isinstance(beam, _BeamsMesh)
    assert beam.n_lines == 48
    assert np.isclose(beam.length, 64.36592438345)


def test_compute_his_conduction():
    # get a fresh model
    fourchamber = get_fourchamber()
    cs = ConductionSystem(fourchamber)

    cs.compute_sa_node()
    cs.compute_av_node()
    beam = cs.compute_av_conduction()
    beam, _, _ = cs.compute_his_conduction()
    his = beam.get_lines_by_name(_ConductionType.HIS.value)

    assert his.n_lines == 20
    assert np.isclose(his.length, 17.961541876776412)


@pytest.mark.xfail(reason="Bachmann bundle not implemented yet.")
def test_compute_bachmann_bundle():
    # get a fresh model
    fourchamber = get_fourchamber()
    cs = ConductionSystem(fourchamber)

    cs.compute_sa_node()
    cs.compute_av_node()
    beam = cs.compute_av_conduction()
    beam, _, _ = cs.compute_his_conduction()

    beam = cs._compute_bachman_bundle(
        start_coord=fourchamber.right_atrium.get_point("SA_node").xyz,
        end_coord=np.array([-34, 163, 413]),
    )

    assert np.all(beam.edges[0] == [108609, 121942])
    assert np.all(beam.edges[-1] == [121994, 94118])


def test_compute_left_right_bundle():
    # get a fresh model
    fourchamber = get_fourchamber()
    cs = ConductionSystem(fourchamber)

    cs.compute_sa_node()
    cs.compute_av_node()
    cs.compute_av_conduction()
    beam, left_end, right_end = cs.compute_his_conduction()
    coord_left = fourchamber.left_ventricle.apex_points[0].xyz
    coord_right = fourchamber.right_ventricle.apex_points[0].xyz
    result = cs.compute_left_right_bundle(
        left_end.xyz, end_coord=coord_left, side=_ConductionType.LEFT_BUNDLE_BRANCH.value
    )
    result = cs.compute_left_right_bundle(
        right_end.xyz, end_coord=coord_right, side=_ConductionType.RIGHT_BUNDLE_BRANCH.value
    )
    left_bundle = result.get_lines_by_name(_ConductionType.LEFT_BUNDLE_BRANCH.value)
    right_bundle = result.get_lines_by_name(_ConductionType.RIGHT_BUNDLE_BRANCH.value)
    assert np.isclose(left_bundle.length, 76.23944721403967)
    assert left_bundle.n_lines == 56
    assert np.isclose(right_bundle.length, 61.4083335120675)
    assert right_bundle.n_lines == 45

    return


def test_connect_to_solid():
    # get a fresh model
    fourchamber = get_fourchamber()
    cs = ConductionSystem(fourchamber)

    cs.compute_sa_node()
    cs.compute_av_node()
    cs.compute_av_conduction()
    _, _, right_end = cs.compute_his_conduction()
    coord_right = fourchamber.right_ventricle.apex_points[0].xyz
    result = cs.compute_left_right_bundle(
        right_end.xyz, end_coord=coord_right, side=_ConductionType.RIGHT_BUNDLE_BRANCH.value
    )
    right_bundle = result.get_lines_by_name(_ConductionType.RIGHT_BUNDLE_BRANCH.value)
    rbb_id = cs.m.conduction_system.get_line_id_from_name(_ConductionType.RIGHT_BUNDLE_BRANCH.value)
    cs._connect_to_solid(rbb_id, [0, -1])
    isconnected = np.zeros(right_bundle.n_points, dtype=int)
    isconnected[[0, -1]] = 1
    assert np.allclose(cs.m.conduction_system.get_lines(rbb_id)["_is-connected"], isconnected)
