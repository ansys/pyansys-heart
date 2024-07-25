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

"""Test for reading purkinje network as a beam mesh. Uses a mock Mesh object."""

import os
import pathlib

import numpy as np
import pytest

from ansys.heart.preprocessor.mesh.objects import BeamMesh, Point
from ansys.heart.preprocessor.models import FourChamber
from tests.heart.conftest import get_assets_folder

model: FourChamber
model_dir: pathlib.Path

pytestmark = pytest.mark.local


@pytest.fixture(autouse=True, scope="module")
def get_data():
    global model, model_dir

    # TODO: file larger than 100 Mb cannot be added to package

    pickle_file = os.path.join(
        get_assets_folder(),
        "reference_models",
        "strocchi2020",
        "01",
        "FullHeart2",
        "heart_model.pickle",
    )
    model_dir = pathlib.Path(pickle_file)

    model = FourChamber.load_model(model_dir)


@pytest.mark.xfail(reason="Test uses local data.")
def test_add_beam_net():
    """Test reading Purkinje from .k files."""

    node_b = np.array([[0, 0, 0], [10, 10, 10]])

    edges = np.array([[0, 0], [0, 1]])
    mask = np.array([[False, True], [True, True]])  # first node is form solid mesh
    model.add_beam_net(beam_nodes=node_b, edges=edges.copy(), mask=mask, pid=0, name="test")

    # construct mesh to compare against
    beam_mesh = BeamMesh(
        nodes=np.vstack((model.mesh.nodes, BeamMesh.all_beam_nodes)),
        edges=np.array([[0, 69754], [69754, 69755]]),
        beam_nodes_mask=mask,
    )
    beam_mesh.pid = 0
    beam_mesh.name = "test"

    assert model.beam_network[0] == beam_mesh


@pytest.mark.xfail(reason="Test uses local data.")
def test_compute_SA_node():
    p = model.compute_SA_node()
    target = Point(name="SA_node", xyz=[-48.95559814, 108.23159848, 422.91220412], node_id=22056)
    assert p.node_id == target.node_id


@pytest.mark.xfail(reason="Test uses local data.")
def test_compute_AV_node():
    p = model.compute_AV_node()
    target = Point(
        name="AV_node", xyz=np.array([-8.20742556, 106.99512699, 373.32172823]), node_id=28424
    )
    assert p.node_id == target.node_id


@pytest.mark.xfail(reason="Test uses local data.")
def test_compute_av_conduction():
    beam = model.compute_av_conduction()

    assert np.all(beam.edges[0] == np.array([22056, 69754]))
    assert np.all(beam.edges[-1] == np.array([69787, 69788]))


@pytest.mark.xfail(reason="Test uses local data.")
def test_compute_His_conduction():
    model.compute_AV_node()
    model.compute_av_conduction()
    beam, _ = model.compute_His_conduction()

    assert np.all(beam.edges[0] == [69788, 69789])
    assert np.all(beam.edges[-1] == [69793, 69795])


@pytest.mark.xfail(reason="Test uses local data.")
def test__define_hisbundle_start_end_point():
    model.compute_AV_node()
    start_p, bifurcation_p = model._define_hisbundle_start_bifurcation()

    assert np.allclose(start_p, [2.96411229, 107.68694446, 367.01330368])
    assert np.allclose(bifurcation_p, [7.24693771, 106.92654286, 365.52226939])
