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

import importlib.util as importlib_util
import sys

import pytest

from ansys.health.heart.settings.material import ep_material_factory as factory
from ansys.health.heart.settings.material.ep_material import EPSolverType


def _make_dummy_defaults():
    """Create a dummy defaults module with eikonal and monodomain entries."""
    module_name = "ansys.health.heart.settings.defaults.electrophysiology"
    spec = importlib_util.spec_from_loader(module_name, loader=None)
    dummy = importlib_util.module_from_spec(spec)

    dummy.default_myocardium_material_eikonal = {
        "sigma_fiber": 0.7,
        "sigma_sheet": 0.2,
    }
    dummy.default_myocardium_material_monodomain = {
        "sigma_fiber": 0.5,
        "sigma_sheet": 0.1,
    }
    dummy.default_beam_material_eikonal = {"sigma_fiber": 1.0}
    dummy.default_beam_material_monodomain = {"sigma_fiber": 1.0}

    return module_name, dummy


# Mocks ActiveBeam and Active EP Material classes.
class MockActive:
    def __init__(self, **kwargs):
        self.kwargs = kwargs


class MockActiveBeam:
    def __init__(self, **kwargs):
        self.kwargs = kwargs


# Mock classes for testing assign_default_ep_materials
class MockPart:
    """Mock Part class for testing material assignment."""

    def __init__(self, name: str, active: bool = True, ep_material=None):
        self.name = name
        self.active = active
        self.ep_material = ep_material


class MockArtery:
    """Mock Artery class for testing insulator assignment."""

    def __init__(self, name: str, ep_material=None):
        self.name = name
        self.active = False
        self.ep_material = ep_material


class MockConductionPath:
    """Mock ConductionPath class for testing."""

    def __init__(self, name: str, ep_material=None):
        self.name = name
        self.ep_material = ep_material


class MockHeartModel:
    """Mock HeartModel class for testing."""

    def __init__(self, parts=None, conduction_paths=None):
        self.parts = parts or []
        self.conduction_paths = conduction_paths or []


class MockEPMaterial:
    """Mock EP material base class."""

    def __init__(self, **kwargs):
        self.kwargs = kwargs


class MockInsulator:
    """Mock Insulator material class."""

    def __init__(self):
        self.material_type = "Insulator"


def test_get_default_myocardium_material_selects_correct_defaults(monkeypatch):
    module_name, dummy = _make_dummy_defaults()
    monkeypatch.setitem(sys.modules, module_name, dummy)

    # Capture calls to Active
    monkeypatch.setattr(factory, "Active", MockActive)

    # Eikonal via enum
    res = factory.get_default_myocardium_material(EPSolverType.EIKONAL)
    assert isinstance(res, MockActive)
    assert res.kwargs == dummy.default_myocardium_material_eikonal

    # Reaction-Eikonal should behave like Eikonal
    res = factory.get_default_myocardium_material(EPSolverType.REACTION_EIKONAL)
    assert isinstance(res, MockActive)
    assert res.kwargs == dummy.default_myocardium_material_eikonal

    # Monodomain via string input
    res = factory.get_default_myocardium_material("Monodomain")
    assert isinstance(res, MockActive)
    assert res.kwargs == dummy.default_myocardium_material_monodomain


def test_get_default_conduction_system_material_selects_correct_defaults(monkeypatch):
    module_name, dummy = _make_dummy_defaults()
    monkeypatch.setitem(sys.modules, module_name, dummy)

    monkeypatch.setattr(factory, "ActiveBeam", MockActiveBeam)

    # Eikonal via enum
    res = factory.get_default_conduction_system_material(EPSolverType.EIKONAL)
    assert isinstance(res, MockActiveBeam)
    assert res.kwargs == dummy.default_beam_material_eikonal

    # Monodomain via string
    res = factory.get_default_conduction_system_material("Monodomain")
    assert isinstance(res, MockActiveBeam)
    assert res.kwargs == dummy.default_beam_material_monodomain


def test_assign_default_ep_materials_with_enum_solver_type(monkeypatch):
    """Test assign_default_ep_materials with EPSolverType enum input.

    This test verifies that the function correctly handles enum solver types
    and assigns appropriate materials to different part types.
    """
    # Mock the factory functions
    mock_active_material = MockEPMaterial(type="active")
    mock_passive_material = MockEPMaterial(type="passive")
    mock_beam_material = MockEPMaterial(type="beam")
    mock_insulator = MockInsulator()

    monkeypatch.setattr(factory, "get_default_myocardium_material", lambda x: mock_active_material)
    monkeypatch.setattr(factory, "get_default_passive_material", lambda x: mock_passive_material)
    monkeypatch.setattr(
        factory, "get_default_conduction_system_material", lambda x: mock_beam_material
    )
    monkeypatch.setattr(factory, "Insulator", lambda: mock_insulator)
    monkeypatch.setattr(factory, "Artery", MockArtery)

    # Create mock model with different part types
    active_part = MockPart("left_ventricle", active=True)
    passive_part = MockPart("scar_tissue", active=False)
    artery_part = MockArtery("coronary_artery")
    conduction_path = MockConductionPath("purkinje_network")

    model = MockHeartModel(
        parts=[active_part, passive_part, artery_part],
        conduction_paths=[conduction_path],
    )

    # Test with enum solver type
    factory.assign_default_ep_materials(model, EPSolverType.MONODOMAIN)

    # Verify materials were assigned correctly
    assert active_part.ep_material is mock_active_material
    assert passive_part.ep_material is mock_passive_material
    assert artery_part.ep_material is mock_insulator
    assert conduction_path.ep_material is mock_beam_material


def test_assign_default_ep_materials_with_string_solver_type(monkeypatch):
    """Test assign_default_ep_materials with string solver type input.

    This test verifies that the function correctly converts string solver types
    to enums before processing.
    """
    # Mock the factory functions
    mock_material = MockEPMaterial(type="active")

    monkeypatch.setattr(factory, "get_default_myocardium_material", lambda x: mock_material)
    monkeypatch.setattr(factory, "get_default_passive_material", lambda x: mock_material)
    monkeypatch.setattr(factory, "get_default_conduction_system_material", lambda x: mock_material)

    # Create mock model
    part = MockPart("test_part", active=True)
    model = MockHeartModel(parts=[part])

    # Test with string solver type
    factory.assign_default_ep_materials(model, "Eikonal")

    # Verify material was assigned
    assert part.ep_material is mock_material


def test_assign_default_ep_materials_with_invalid_solver_type():
    """Test assign_default_ep_materials with invalid solver type.

    This test verifies that the function raises appropriate errors for
    invalid solver type inputs.
    """
    # Create mock model
    part = MockPart("test_part", active=True)
    model = MockHeartModel(parts=[part])

    # Test with invalid string solver type
    with pytest.raises(ValueError):
        factory.assign_default_ep_materials(model, "InvalidSolver")

    # Test with invalid type
    with pytest.raises(ValueError):
        factory.assign_default_ep_materials(model, 123)


def test_assign_default_ep_materials_skips_existing_materials(monkeypatch):
    """Test that assign_default_ep_materials skips parts with existing valid materials.

    This test verifies that the function doesn't overwrite existing EP materials
    when they are already properly assigned.
    """
    # Mock existing material
    existing_material = MockEPMaterial(type="existing")

    # Mock the factory functions (should not be called)
    def should_not_be_called(x):
        pytest.fail("Factory function should not be called for parts with existing materials")

    monkeypatch.setattr(factory, "get_default_myocardium_material", should_not_be_called)
    monkeypatch.setattr(factory, "EPMaterialModel", MockEPMaterial)

    # Create mock part with existing material
    part = MockPart("test_part", active=True, ep_material=existing_material)
    model = MockHeartModel(parts=[part])

    # Function should not modify the existing material
    factory.assign_default_ep_materials(model, EPSolverType.MONODOMAIN)

    # Verify existing material is preserved
    assert part.ep_material is existing_material


def test_assign_default_ep_materials_handles_empty_model(monkeypatch):
    """Test assign_default_ep_materials with empty model.

    This test verifies that the function handles models with no parts
    or conduction paths gracefully.
    """

    # Mock the factory functions (should not be called)
    def should_not_be_called(x):
        pytest.fail("Factory function should not be called for empty model")

    monkeypatch.setattr(factory, "get_default_myocardium_material", should_not_be_called)

    # Create empty model
    model = MockHeartModel(parts=[], conduction_paths=[])

    # Should handle empty model without errors
    factory.assign_default_ep_materials(model, EPSolverType.MONODOMAIN)


def test_assign_default_ep_materials_handles_mixed_part_states(monkeypatch):
    """Test assign_default_ep_materials with mixed part material states.

    This test verifies correct handling when some parts have materials
    and others don't.
    """
    # Mock materials
    mock_active_material = MockEPMaterial(type="active")
    existing_material = MockEPMaterial(type="existing")

    monkeypatch.setattr(factory, "get_default_myocardium_material", lambda x: mock_active_material)
    monkeypatch.setattr(factory, "EPMaterialModel", MockEPMaterial)

    # Create parts with mixed states
    part_no_material = MockPart("part1", active=True, ep_material=None)
    part_with_material = MockPart("part2", active=True, ep_material=existing_material)

    model = MockHeartModel(parts=[part_no_material, part_with_material])

    # Assign materials
    factory.assign_default_ep_materials(model, EPSolverType.MONODOMAIN)

    # Verify correct assignment
    assert part_no_material.ep_material is mock_active_material
    assert part_with_material.ep_material is existing_material


@pytest.mark.parametrize(
    "solver_type",
    [
        EPSolverType.MONODOMAIN,
        EPSolverType.EIKONAL,
        EPSolverType.REACTION_EIKONAL,
        "Monodomain",
        "Eikonal",
        "ReactionEikonal",
    ],
)
def test_assign_default_ep_materials_all_solver_types(monkeypatch, solver_type):
    """Test assign_default_ep_materials with all valid solver types.

    This parametrized test verifies that the function works correctly
    with all supported solver types (both enum and string forms).
    """
    # Mock materials
    mock_material = MockEPMaterial(type="test")

    monkeypatch.setattr(factory, "get_default_myocardium_material", lambda x: mock_material)
    monkeypatch.setattr(factory, "get_default_passive_material", lambda x: mock_material)
    monkeypatch.setattr(factory, "get_default_conduction_system_material", lambda x: mock_material)

    # Create test model
    active_part = MockPart("active_part", active=True)
    passive_part = MockPart("passive_part", active=False)
    conduction_path = MockConductionPath("test_path")

    model = MockHeartModel(
        parts=[active_part, passive_part],
        conduction_paths=[conduction_path],
    )

    # Test assignment
    factory.assign_default_ep_materials(model, solver_type)

    # Verify all materials were assigned
    assert active_part.ep_material is mock_material
    assert passive_part.ep_material is mock_material
    assert conduction_path.ep_material is mock_material
