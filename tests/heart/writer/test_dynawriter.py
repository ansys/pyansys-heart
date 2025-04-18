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
import os
import tempfile
import unittest.mock as mock

import numpy as np
import pytest
import pyvista as pv
import pyvista.examples as examples

from ansys.health.heart.models import FullHeart
from ansys.health.heart.objects import (
    Mesh,
    Part,
    PartType,
    Point,
)
from ansys.health.heart.pre.conduction_beams import ConductionPath, ConductionPathType
from ansys.health.heart.settings.settings import Mechanics, SimulationSettings, Stimulation
import ansys.health.heart.writer.dynawriter as writers


def _get_mock_conduction_system() -> Mesh:
    """Get a mock conduction system."""
    edges = examples.load_tetbeam().extract_feature_edges()
    conduction_system = Mesh(edges)
    conduction_system.point_data["_is-connected"] = 0
    conduction_system.cell_data["_line-id"] = 10
    conduction_system.point_data["_is-connected"][0:10] = 1
    conduction_system.point_data["_shifted_id"] = 1

    return conduction_system


@pytest.fixture()
def _mock_model():
    """Create a mock EP Writer."""
    model: FullHeart = mock.MagicMock(FullHeart)
    model.mesh = Mesh(examples.load_tetbeam())

    p1 = Point(name="electrode-1", xyz=np.array([0.0, 0.0, 0.0]), node_id=1)
    p2 = Point(name="electrode-1", xyz=np.array([1.0, 1.0, 1.0]), node_id=1)

    model.electrodes = [p1, p2]

    lines = pv.line_segments_from_points([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    model.conduction_paths = [
        ConductionPath(ConductionPathType.LEFT_PURKINJE, lines, 1, [False, False], pv.Sphere())
    ]

    model.conduction_mesh = _get_mock_conduction_system()

    yield model


def _add_parts(model: FullHeart):
    """Add parts to model."""
    model.parts = [
        Part(name="left_ventricle", part_type=PartType.VENTRICLE),
        Part(name="Right ventricle", part_type=PartType.VENTRICLE),
        Part(name="Septum", part_type=PartType.SEPTUM),
    ]
    for ii, part in enumerate(model.parts):
        part.pid = ii
    return model


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
    writer._add_segment_from_surface(pv.Sphere(), name="test")

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


@pytest.mark.parametrize(
    "solvertype,expected_num_keywords",
    [
        ("Monodomain", 6),
        ("Eikonal", 6),
        ("ReactionEikonal", 6),
    ],
)
def test_update_ep_settings(_mock_model, solvertype, expected_num_keywords):
    """Test updating EP settings."""
    model = _mock_model

    model = _add_parts(model)

    settings = SimulationSettings()
    settings.load_defaults()
    settings.electrophysiology.analysis.solvertype = solvertype

    writer = writers.ElectroMechanicsDynaWriter(model, settings)
    writer._update_ep_settings(beam_pid=[1, 2, 3])

    assert len(writer.kw_database.ep_settings.keywords) == expected_num_keywords

    del writer


def test_add_myocardial_nodeset_layer(_mock_model):
    """Test adding myocardial nodeset layer."""
    model = _mock_model
    # add a dummy transmural direction.
    model.mesh.point_data["transmural"] = model.mesh.points[:, -1] / np.max(
        model.mesh.points[:, -1]
    )

    settings = SimulationSettings()
    settings.load_defaults()

    writer = writers.ElectroMechanicsDynaWriter(model, settings)
    # assert node-set ids (no other nodesets present, so expecting 1,2,3)
    assert writer._create_myocardial_nodeset_layers() == (1, 2, 3)


@pytest.mark.xfail(reason="_update_create_fibers un-testable, needs refactoring.")
def test_update_create_fibers(_mock_model):
    """Test updating the create fibers method."""
    model = _mock_model
    settings = SimulationSettings()
    settings.load_defaults()

    writer = writers.FiberGenerationDynaWriter(model, settings)

    writer._update_create_fibers()


def test_update_use_purkinje(_mock_model: FullHeart):
    """Test update use purkinje."""
    model = _mock_model
    model = _add_parts(model)

    settings = SimulationSettings()
    settings.load_defaults()

    writer = writers.ElectroMechanicsDynaWriter(model, settings)

    writer._update_use_Purkinje()

    kw_titles = [kw.get_title() for kw in writer.kw_database.beam_networks.keywords]
    expected_kw_titles = [
        "*SECTION_BEAM",
        "*NODE",
        "*EM_EP_PURKINJE_NETWORK2",
        "*PART",
        "*MAT_NULL",
        "*EM_MAT_001",
        "*ELEMENT_BEAM",
    ]
    for expected_kw in expected_kw_titles:
        assert expected_kw in kw_titles, f"Did not find {expected_kw} in keywords"


def test_export(_mock_model):
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        setting = mock.Mock(spec=Mechanics).return_value
        setting.mechanics.system.name = "ConstantPreloadWindkesselAfterload"
        w = writers.MechanicsDynaWriter(_mock_model, setting)
        w.kw_database.main.append("*END")
        w.export(tempdir)
        # test export
        assert os.listdir(tempdir) == ["main.k"]

        user_file = os.path.join(tempdir, "..", "user.k")
        with open(user_file, "w") as tmpfile:
            tmpfile.write("*END")

        # test raise with not found file
        with pytest.raises(FileNotFoundError):
            w.export(tempdir, user_k=[user_file + "0"])
        # test write with user file
        w.export(tempdir, user_k=[user_file])
        assert set(os.listdir(tempdir)) == {"main.k", "user.k"}
