"""Test for reading purkinje network as a beam mesh. Uses a mock Mesh object."""
import os

from ansys.heart.preprocessor.mesh.objects import BeamMesh, Mesh, Point
from ansys.heart.preprocessor.models import FourChamber
import numpy as np
import pytest
from pyvista import examples

from .conftest import get_workdir

model: FourChamber


@pytest.fixture(autouse=True, scope="module")
def get_data():
    global model

    model_dir = r"D:\pyheart-lib\test_case\test_4c\heart_model.pickle"
    model = FourChamber.load_model(model_dir)


def test_read_purkinje_from_kfile_001():
    """Test reading Purkinje from .k files."""

    # prepare example (cast directly to Mesh object)
    grid: Mesh = Mesh(examples.load_tetbeam())

    # assume beams are using already existed nodes except the last one
    offset = grid.nodes.shape[0] + 1
    beams = np.array(
        [[3, 12], [12, 24], [24, 33], [12, 18], [12, 21], [21, 30], [3, offset]], dtype=int
    )

    # add additional point to mesh
    extra_nodes = [[0.5, 0.5, -1], [0.5, 0.5, -0.5]]
    # grid.nodes = np.vstack([grid.points, extra_nodes])

    # mock .k file
    keyword = "*KEYWORD\n"
    keyword += "*NODE +\n"
    keyword += "{:>20d}{:>20f}{:>20f}{:>20f}\n".format(offset, *extra_nodes[0])
    keyword += "{:>20d}{:>20f}{:>20f}{:>20f}\n".format(offset + 1, *extra_nodes[1])
    keyword += "*ELEMENT_BEAM\n"
    # keyword += "$#   eid     pid      n1      n2\n"
    elid = 1
    for beam in beams:
        keyword = keyword + "{:>8d}{:>8d}{:>8d}{:>8d}\n".format(elid, 4, *beam + 1)
        elid += 1
    keyword += "*NODE\n"
    keyword += "*END\n"

    # write to disc
    temp_k_file = os.path.join(get_workdir(), "test_beams.k")
    with open(temp_k_file, "w") as f:
        f.write(keyword)

    grid.add_purkinje_from_kfile(temp_k_file)

    # construct mesh to compare against
    beam_mesh = BeamMesh()
    beam_mesh.pid = 4
    beam_mesh.id = 1
    beam_mesh.nodes = grid.nodes
    beam_mesh.edges = beams
    assert grid.beam_network[0] == beam_mesh

    # make sure IDs are right for multi calls
    grid.add_purkinje_from_kfile(temp_k_file)

    beam_mesh2 = BeamMesh()
    beam_mesh2.pid = 5
    beam_mesh2.id = 2
    beam_mesh2.nodes = grid.nodes
    beams[-1, -1] += 2
    beam_mesh2.edges = beams

    assert grid.beam_network[1] == beam_mesh2


def test_compute_SA_node():
    p = model.compute_SA_node()
    target = Point(name="SA_node", xyz=[-49.53661854, 108.12227932, 422.97088272], node_id=19705)
    assert p.node_id == target.node_id


def test_compute_AV_node():
    p = model.compute_AV_node()
    target = Point(name="AV_node", xyz=[-8.22263189, 106.95353898, 373.34239855], node_id=26409)
    assert p.node_id == target.node_id
