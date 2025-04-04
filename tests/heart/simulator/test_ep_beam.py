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

from ansys.heart.core.objects import _ConductionType
from ansys.heart.core.pre.conduction_beam import ConductionSystem
from tests.heart.conftest import get_fourchamber


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
