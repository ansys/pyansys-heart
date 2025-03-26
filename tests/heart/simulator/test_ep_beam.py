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

from ansys.heart.core.objects import _BeamMesh
from ansys.heart.preprocessor.conduction_beam import ConductionSystem
from tests.heart.conftest import get_fourchamber


@pytest.mark.xfail(reason="add_beam_net method removed.")
def test_add_beam_net():
    """Test reading Purkinje from .k files."""
    fourchamber = get_fourchamber()

    nodes = np.array([[0, 0, 0], [10, 10, 10]])
    edges = np.array([[0, 0], [0, 1]])
    mask = np.array([[False, True], [True, True]])  # first node is form solid mesh
    fourchamber.add_beam_net(beam_nodes=nodes, edges=edges.copy(), mask=mask, pid=0, name="test")

    # construct mesh to compare against
    n = fourchamber.mesh.points.shape[0]
    beam_mesh = _BeamMesh(
        nodes=np.vstack((fourchamber.mesh.points, _BeamMesh.all_beam_nodes)),
        edges=np.array([[0, n], [n, n + 1]]),
        beam_nodes_mask=mask,
    )
    beam_mesh.pid = 0
    beam_mesh.name = "test"
    assert fourchamber.beam_network[0] == beam_mesh


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

    # based on type hint: will fail:
    # assert isinstance(beam, _BeamsMesh)
    assert beam.n_lines == 48
    assert np.all(beam.get_cell(0).point_ids == [0, 1])
    assert np.all(beam.get_cell(beam.n_cells - 1).point_ids == [47, 48])


def test_compute_his_conduction():
    # get a fresh model
    fourchamber = get_fourchamber()
    cs = ConductionSystem(fourchamber)

    cs.compute_sa_node()
    cs.compute_av_node()
    beam = cs.compute_av_conduction()
    beam, _, _ = cs.compute_his_conduction()
    his = beam.get_lines_by_name("His")

    assert np.all(his.get_cell(0).point_ids == [0, 1])
    assert np.all(his.get_cell(his.n_cells - 1).point_ids == [19, 20])
    assert np.isclose(his.length, 17.961541876776412)


@pytest.mark.xfail(reason="add_beam_net method was removed")
def test_compute_bachmann_bundle():
    # get a fresh model
    fourchamber = get_fourchamber()
    cs = ConductionSystem(fourchamber)

    cs.compute_sa_node()
    cs.compute_av_node()
    beam = cs.compute_av_conduction()
    beam, _, _ = cs.compute_his_conduction()

    beam = cs.compute_bachman_bundle(
        start_coord=fourchamber.right_atrium.get_point("SA_node").xyz,
        end_coord=np.array([-34, 163, 413]),
    )

    assert np.all(beam.edges[0] == [108609, 121942])
    assert np.all(beam.edges[-1] == [121994, 94118])


@pytest.mark.xfail(reason="Requires Purkinje network to be available.")
def test_compute_left_right_bundle():
    # get a fresh model
    fourchamber = get_fourchamber()
    cs = ConductionSystem(fourchamber)

    cs.compute_sa_node()
    cs.compute_av_node()
    cs.compute_av_conduction()
    beam, left_end, right_end = cs.compute_his_conduction()

    left_bundle = cs.compute_left_right_bundle(left_end.xyz, side="Left")
    right_bundle = cs.compute_left_right_bundle(right_end.xyz, side="Right")

    # fourchamber.plot_purkinje()
    assert np.all(left_bundle.edges[0] == [121937, 121942])
    assert np.all(left_bundle.edges[-1] == [121996, 50025])
    assert np.all(right_bundle.edges[0] == [121941, 121997])
    assert np.all(right_bundle.edges[-1] == [122040, 68079])

    return
