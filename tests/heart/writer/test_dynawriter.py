import unittest.mock as mock

import numpy as np
import pytest
import pyvista as pv
import pyvista.examples as examples

from ansys.heart.core.models import FullHeart
from ansys.heart.core.objects import Mesh, Point
from ansys.heart.simulator.settings.settings import SimulationSettings, Stimulation
import ansys.heart.writer.dynawriter as writers


@pytest.fixture()
def _mock_model():
    """Create a mock EP Writer."""
    model: FullHeart = mock.MagicMock(FullHeart)
    model.mesh = Mesh(examples.load_tetbeam())

    p1 = Point(name="electrode-1", xyz=np.array([0.0, 0.0, 0.0]), node_id=1)
    p2 = Point(name="electrode-1", xyz=np.array([1.0, 1.0, 1.0]), node_id=1)

    model.electrodes = [p1, p2]

    yield model


def test_update_ECG_coordinates(_mock_model):  # noqa: N802
    """Test updating ECG Coordinates."""
    model = _mock_model

    writer = writers.ElectroMechanicsDynaWriter(model)

    writer._update_ECG_coordinates()

    strings_kw = " ".join(writer.kw_database.ep_settings.string_keywords)

    assert len(writer.kw_database.ep_settings.string_keywords) == 4
    assert "*EM_POINT_SET" in strings_kw
    assert "*EM_EP_EKG" in strings_kw


def test_add_segment_from_surface(_mock_model):
    """Test adding a segment set from a Mesh surface."""
    model = _mock_model
    writer = writers.ElectroMechanicsDynaWriter(model)
    model.mesh.add_surface(pv.Sphere(), name="test", id=1)
    writer._add_segment_from_surface(name="test")

    assert len(writer.kw_database.segment_sets.keywords) == 1
    assert writer.kw_database.segment_sets.keywords[0].get_title() == "*SET_SEGMENT_TITLE"

    pass


@pytest.mark.parametrize(
    "solvertype,expected_kw",
    [
        ("Monodomain", "*EM_EP_TENTUSSCHER_STIMULUS"),
        ("Eikonal", "*EM_EP_EIKONAL"),
        ("ReactionEikonal", "*EM_EP_EIKONAL"),
    ],
)
def test_add_stimulation_keyword(_mock_model, solvertype, expected_kw):
    """Test adding a stimulation keyword."""
    model = _mock_model

    settings = SimulationSettings()
    settings.load_defaults()
    settings.electrophysiology.analysis.solvertype = solvertype
    # set up stimulation
    stimulation = Stimulation([1, 2])

    writer = writers.ElectroMechanicsDynaWriter(model, settings)

    nodeset_kw, stim_kw = writer._add_stimulation_keyword(stimulation)

    if solvertype == "Monodomain":
        assert stim_kw.get_title() == expected_kw
    else:
        assert "*EM_EP_EIKONAL" in stim_kw
