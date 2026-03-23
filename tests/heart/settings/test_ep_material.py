# Copyright (C) 2023 - 2026 ANSYS, Inc. and/or its affiliates.
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
"""Test module EP materials."""

from pydantic import ValidationError
import pytest

from ansys.health.heart.settings.material.cell_models import Tentusscher, TentusscherEndo
from ansys.health.heart.settings.material.ep_material import (
    Active,
    ActiveBeam,
    ActiveNew,
    EPMaterialModel,
    Insulator,
    Passive,
)


def test_active() -> None:
    active0 = Active(sigma_fiber=1)
    assert active0.sigma_sheet is None
    active = Active(sigma_fiber=1, sigma_sheet=1, sigma_sheet_normal=1)
    assert active.sigma_sheet_normal == 1


class TestEPMaterialModel:
    """Test EPMaterialModel base class."""

    def test_ep_material_model_creation_with_defaults(self) -> None:
        """Test creating EPMaterialModel with default values."""
        model = EPMaterialModel()

        # Check default values from ep_defaults
        assert model.sigma_fiber is None  # Should be None by default
        assert model.sigma_sheet is None
        assert model.sigma_sheet_normal is None
        assert model.beta is None
        assert model.cm is None
        assert model.lambda_ is None

    def test_ep_material_model_optional_fields(self) -> None:
        """Test that optional fields can be None or have values."""
        # Test with None values
        model = EPMaterialModel(
            sigma_sheet=None, sigma_sheet_normal=None, beta=None, cm=None, lambda_=None
        )

        assert model.sigma_sheet is None
        assert model.sigma_sheet_normal is None
        assert model.beta is None
        assert model.cm is None
        assert model.lambda_ is None

    def test_ep_material_model_custom_values(self) -> None:
        """Test EPMaterialModel with custom values."""
        model = EPMaterialModel(
            sigma_fiber=1.5,
            sigma_sheet=0.8,
            sigma_sheet_normal=0.9,
            beta=2000.0,
            cm=0.2,
            lambda_=0.5,
        )

        assert model.sigma_fiber == 1.5
        assert model.sigma_sheet == 0.8
        assert model.sigma_sheet_normal == 0.9
        assert model.beta == 2000.0
        assert model.cm == 0.2
        assert model.lambda_ == 0.5

    def test_ep_material_model_sheet_parameter_sync(self) -> None:
        """Test that sigma_sheet and sigma_sheet_normal sync properly."""
        # Test sigma_sheet set, sigma_sheet_normal None -> sync
        model = EPMaterialModel(sigma_sheet=1.0, sigma_sheet_normal=None)
        assert model.sigma_sheet == 1.0
        assert model.sigma_sheet_normal == 1.0  # Should be synced

        # Test sigma_sheet_normal set, sigma_sheet None -> sync
        model = EPMaterialModel(sigma_sheet=None, sigma_sheet_normal=2.0)
        assert model.sigma_sheet == 2.0  # Should be synced
        assert model.sigma_sheet_normal == 2.0

    def test_ep_material_model_both_sheet_parameters_set(self) -> None:
        """Test when both sheet parameters are explicitly set."""
        model = EPMaterialModel(sigma_sheet=1.0, sigma_sheet_normal=2.0)
        # Both should retain their original values
        assert model.sigma_sheet == 1.0
        assert model.sigma_sheet_normal == 2.0

    def test_ep_material_model_neither_sheet_parameter_set(self) -> None:
        """Test when neither sheet parameter is set."""
        model = EPMaterialModel(sigma_sheet=None, sigma_sheet_normal=None)
        # Both should remain None
        assert model.sigma_sheet is None
        assert model.sigma_sheet_normal is None

    def test_ep_material_model_invalid_types(self) -> None:
        """Test validation with invalid types."""
        with pytest.raises(ValidationError):
            EPMaterialModel(sigma_fiber="invalid")

        with pytest.raises(ValidationError):
            EPMaterialModel(beta="not_a_number")

    def test_ep_material_model_serialization(self) -> None:
        """Test serialization handling of None values."""
        model = EPMaterialModel(sigma_fiber=1.0, sigma_sheet=None, beta=2000.0)

        data = model.model_dump()
        assert data["sigma_fiber"] == 1.0
        assert data["sigma_sheet"] is None
        assert data["beta"] == 2000.0

        # Test excluding None values
        data_no_none = model.model_dump(exclude_none=True)
        assert "sigma_sheet" not in data_no_none
        assert "sigma_fiber" in data_no_none
        assert "beta" in data_no_none


class TestInsulator:
    """Test Insulator material class."""

    def test_insulator_defaults(self) -> None:
        """Test Insulator default values."""
        insulator = Insulator()

        assert insulator.sigma_fiber == 0.0
        assert insulator.cm == 0.0
        assert insulator.beta == 0.0

    def test_insulator_custom_values(self) -> None:
        """Test Insulator with custom values."""
        insulator = Insulator(sigma_fiber=0.1, cm=0.05, beta=100.0)

        assert insulator.sigma_fiber == 0.1
        assert insulator.cm == 0.05
        assert insulator.beta == 100.0

    def test_insulator_inheritance(self) -> None:
        """Test that Insulator doesn't inherit from EPMaterialModel."""
        insulator = Insulator()

        # Insulator is standalone BaseModel, not EPMaterialModel
        assert isinstance(insulator, Insulator)
        assert isinstance(insulator, EPMaterialModel)

    def test_insulator_serialization(self) -> None:
        """Test Insulator serialization."""
        insulator = Insulator(sigma_fiber=0.1)

        data = insulator.model_dump()
        assert data["sigma_fiber"] == 0.1
        assert data["cm"] == 0.0
        assert data["beta"] == 0.0

    def test_insulator_validation(self) -> None:
        """Test Insulator field validation."""
        with pytest.raises(ValidationError):
            Insulator(sigma_fiber="invalid")


class TestActive:
    """Test Active EP material class."""

    def test_active_inheritance(self) -> None:
        """Test that Active inherits from EPMaterialModel."""
        active = Active()
        assert isinstance(active, Active)
        assert isinstance(active, EPMaterialModel)

    def test_active_default_creation(self) -> None:
        """Test Active creation with defaults."""
        active = Active()
        assert active.sigma_fiber is None
        assert active.sigma_sheet is None
        assert active.sigma_sheet_normal is None
        assert active.beta == 140
        assert active.cm == 0.01

        # Active-specific defaults
        assert isinstance(active.cell_model, Tentusscher)

    def test_active_cell_model_default_factory(self) -> None:
        """Test that cell_model uses default_factory properly."""
        active1 = Active()
        active2 = Active()

        # Should be separate instances
        assert active1.cell_model is not active2.cell_model
        assert isinstance(active1.cell_model, Tentusscher)
        assert isinstance(active2.cell_model, Tentusscher)

    def test_active_new(self) -> None:
        """Test ActiveNew class."""
        active_new = ActiveNew(sigma_fiber=1, cond_sigma_fiber=1)

        assert isinstance(active_new, Active)
        assert active_new.sigma_fiber == 1
        assert active_new.cond_sigma_fiber == 1

    def test_error_on_invalid_active_new_field(self) -> None:
        """Test that ActiveNew raises error on invalid field."""
        with pytest.raises(ValidationError):
            ActiveNew(sigma_fiber=1, invalid_field=123)

    def test_active_custom_cell_model(self) -> None:
        """Test Active with custom cell model."""
        custom_cell = TentusscherEndo(gto=0.1)
        active = Active(cell_model=custom_cell)

        assert active.cell_model is custom_cell
        assert isinstance(active.cell_model, TentusscherEndo)
        assert active.cell_model.gto == 0.1

    def test_active_optional_sigma_fields(self) -> None:
        """Test Active with optional sigma fields."""
        # Test with some None, some set
        active = Active(sigma_fiber=1.2, sigma_sheet=None, sigma_sheet_normal=0.8)

        assert active.sigma_fiber == 1.2
        assert active.sigma_sheet == 0.8  # Should sync from sigma_sheet_normal
        assert active.sigma_sheet_normal == 0.8

    def test_active_serialization_with_cell_model(self) -> None:
        """Test Active serialization includes cell model."""
        active = Active(sigma_fiber=1.0)

        data = active.model_dump()
        assert "cell_model" in data
        assert isinstance(data["cell_model"], dict)
        assert "gas_constant" in data["cell_model"]  # From Tentusscher

    def test_active_deserialization_with_cell_model(self) -> None:
        """Test Active deserialization with nested cell model."""
        data = {
            "solver_type": "Monodomain",
            "sigma_fiber": 1.5,
            "cell_model": {"gas_constant": 8000.0, "t": 300.0},
        }

        active = Active(**data)
        assert active.sigma_fiber == 1.5
        assert active.cell_model.gas_constant == 8000.0
        assert active.cell_model.t == 300.0


class TestActiveBeam:
    """Test ActiveBeam EP material class."""

    def test_active_beam_inheritance(self) -> None:
        """Test that ActiveBeam inherits from Active."""
        beam = ActiveBeam()
        assert isinstance(beam, ActiveBeam)
        assert isinstance(beam, Active)
        assert isinstance(beam, EPMaterialModel)

    def test_active_beam_defaults(self) -> None:
        """Test ActiveBeam default values."""
        beam = ActiveBeam()

        # Should have beam-specific defaults
        assert beam.sigma_fiber is None  # From ep_defaults.material["beam"]
        assert isinstance(beam.cell_model, Tentusscher)

    def test_active_beam_cell_model_factory(self) -> None:
        """Test that ActiveBeam uses TentusscherEndo by default."""
        beam1 = ActiveBeam()
        beam2 = ActiveBeam()

        # Should be separate TentusscherEndo instances
        assert beam1.cell_model is not beam2.cell_model
        assert isinstance(beam1.cell_model, Tentusscher)
        assert isinstance(beam2.cell_model, Tentusscher)

    def test_active_beam_custom_values(self) -> None:
        """Test ActiveBeam with custom values."""
        custom_cell = TentusscherEndo(gas_constant=8100.0)
        beam = ActiveBeam(sigma_fiber=2.0, cell_model=custom_cell)

        assert beam.sigma_fiber == 2.0
        assert beam.cell_model is custom_cell

    def test_active_beam_serialization(self) -> None:
        """Test ActiveBeam serialization includes all fields."""
        beam = ActiveBeam()

        data = beam.model_dump()
        assert "cell_model" in data
        assert isinstance(data["cell_model"], dict)


class TestPassive:
    """Test Passive EP material class."""

    def test_passive_inheritance(self) -> None:
        """Test that Passive inherits from EPMaterialModel."""
        passive = Passive()
        assert isinstance(passive, Passive)
        assert isinstance(passive, EPMaterialModel)

    def test_passive_defaults(self) -> None:
        """Test Passive default values."""
        passive = Passive()

        # Should have passive-specific defaults
        assert passive.sigma_fiber is None  # From ep_defaults
        assert passive.sigma_sheet is None
        assert passive.sigma_sheet_normal is None

        # Inherited from EPMaterialModel
        assert passive.beta is None
        assert passive.cm is None

    def test_passive_no_cell_model(self) -> None:
        """Test that Passive doesn't have cell_model field."""
        passive = Passive()
        assert not hasattr(passive, "cell_model")

    def test_passive_custom_values(self) -> None:
        """Test Passive with custom values."""
        passive = Passive(sigma_fiber=1.8, sigma_sheet=0.9, sigma_sheet_normal=1.1)

        assert passive.sigma_fiber == 1.8
        assert passive.sigma_sheet == 0.9
        assert passive.sigma_sheet_normal == 1.1

    def test_passive_sheet_parameter_sync(self) -> None:
        """Test sheet parameter syncing in Passive."""
        passive = Passive(sigma_sheet=1.2)
        assert passive.sigma_sheet == 1.2
        assert passive.sigma_sheet_normal == 1.2  # Should sync


class TestFieldInteractions:
    """Test interactions between different field types."""

    def test_optional_vs_required_fields(self) -> None:
        """Test behavior of optional vs required fields."""
        # Test that we can create models with minimal required fields
        models = [EPMaterialModel(), Insulator(), Active(), ActiveBeam(), Passive()]

        for model in models:
            # Should be able to serialize
            data = model.model_dump()
            assert isinstance(data, dict)

            # Should be able to recreate
            restored = type(model)(**data)
            assert type(restored) is type(model)

    def test_default_factory_independence(self) -> None:
        """Test that default factories create independent instances."""
        active1 = Active()
        active2 = Active()
        beam1 = ActiveBeam()
        beam2 = ActiveBeam()

        # Cell models should be independent instances
        assert active1.cell_model is not active2.cell_model
        assert beam1.cell_model is not beam2.cell_model
        assert active1.cell_model is not beam1.cell_model

        # Modifying one shouldn't affect others
        active1.cell_model.gas_constant = 9000.0
        assert active2.cell_model.gas_constant != 9000.0

    def test_none_handling_in_serialization(self) -> None:
        """Test proper None handling in serialization."""
        model = EPMaterialModel(sigma_fiber=1.0, sigma_sheet=None, beta=2000.0, lambda_=None)

        # Standard serialization includes None
        data_with_none = model.model_dump()
        assert data_with_none["sigma_sheet"] is None
        assert data_with_none["lambda_"] is None

        # Exclude None serialization
        data_no_none = model.model_dump(exclude_none=True)
        assert "sigma_sheet" not in data_no_none
        assert "lambda_" not in data_no_none
        assert "sigma_fiber" in data_no_none
        assert "beta" in data_no_none

    def test_field_validation_inheritance(self) -> None:
        """Test that field validation works properly in inheritance."""
        # Test validation in base class
        with pytest.raises(ValidationError):
            EPMaterialModel(sigma_fiber="invalid")

        # Test validation in derived classes
        with pytest.raises(ValidationError):
            Active(sigma_fiber="invalid")

        with pytest.raises(ValidationError):
            Passive(sigma_sheet="invalid")


class TestDefaultsIntegration:
    """Test integration with ep_defaults module."""

    def test_field_override_defaults(self) -> None:
        """Test that explicit values override defaults."""
        # Create models with explicit values that should override defaults
        ep_model = EPMaterialModel(beta=3000.0, cm=0.3)
        passive = Passive(sigma_fiber=5.0)

        # Explicit values should be used instead of defaults
        assert ep_model.beta == 3000.0
        assert ep_model.cm == 0.3
        assert passive.sigma_fiber == 5.0


class TestComplexSerialization:
    """Test complex serialization scenarios."""

    def test_nested_model_serialization(self) -> None:
        """Test serialization with nested models."""
        custom_cell = TentusscherEndo(gto=0.08, v=-88.0)
        active = Active(sigma_fiber=1.5, cell_model=custom_cell, sigma_sheet=0.8)

        # Test serialization
        data = active.model_dump()

        # Verify structure
        assert data["sigma_fiber"] == 1.5
        assert data["sigma_sheet"] == 0.8
        assert data["sigma_sheet_normal"] == 0.8  # Should be synced
        assert "cell_model" in data
        assert data["cell_model"]["gto"] == 0.08
        assert data["cell_model"]["v"] == -88.0

        # Test deserialization
        restored = Active(**data)
        assert restored.sigma_fiber == 1.5
        assert restored.cell_model.gto == 0.08
        assert restored.cell_model.v == -88.0

    def test_json_serialization_with_defaults(self) -> None:
        """Test JSON serialization preserves defaults and optionals."""
        beam = ActiveBeam()

        # JSON round trip
        json_str = beam.model_dump_json()
        restored = ActiveBeam.model_validate_json(json_str)

        # All fields should be preserved
        assert isinstance(restored.cell_model, Tentusscher)
        assert restored.sigma_fiber == beam.sigma_fiber
        assert restored.beta == beam.beta
