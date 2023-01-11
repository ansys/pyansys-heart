"""Test for reading purkinje network as a beam mesh. Uses a mock Mesh object."""
import pyvista as pv
from pyvista import examples
import numpy as np
from ansys.heart.preprocessor.mesh.objects import Mesh, BeamMesh


def test_read_purkinje_from_kfile_001():
    """Test reading Purkinje from .k files."""

    # prepare example (cast directly to Mesh object)
    grid: Mesh = Mesh(examples.load_tetbeam())

    offset = grid.nodes.shape[0] + 1
    beams = np.array(
        [[3, 12], [12, 24], [24, 33], [12, 18], [12, 21], [21, 30], [3, offset]], dtype=int
    )

    # add additional point to mesh
    extra_nodes = [[0.5, 0.5, -1], [0.5, 0.5, -0.5]]
    grid.nodes = np.vstack([grid.points, extra_nodes])

    # construct mesh to compare against
    beam_mesh = BeamMesh()
    beam_mesh.pid = 4
    beam_mesh.id = 1
    beam_mesh.nodes = grid.nodes
    beam_mesh.edges = beams

    # mock .k file
    keyword = "*KEYWORD\n"
    keyword += "*NODE +\n"
    keyword += "{:>20d}{:>20f}{:>20f}{:>20f}\n".format(offset, *extra_nodes[0])
    keyword += "{:>20d}{:>20f}{:>20f}{:>20f}\n".format(offset + 1, *extra_nodes[1])
    keyword += "*ELEMENT_BEAM\n"
    keyword += "$#   eid     pid      n1      n2\n"
    elid = 1
    for beam in beams:
        keyword = keyword + "{:>8d}{:>8d}{:>8d}{:>8d}\n".format(elid, 4, *beam + 1)
        elid += 1
    keyword += "*NODE\n"
    keyword += "*END\n"

    # write to disc
    with open("test_beams.k", "w") as f:
        f.write(keyword)

    # start test:
    grid.add_purkinje_from_kfile("test_beams.k")

    # assert whether result is as expected
    assert grid.beam_network[0] == beam_mesh
