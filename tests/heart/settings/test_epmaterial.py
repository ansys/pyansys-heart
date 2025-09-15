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
"""Test module EP materials."""

from unittest.mock import MagicMock, patch

from pydantic import ValidationError
import pytest

import ansys.health.heart.settings.material.cell_models as cell_models
from ansys.health.heart.settings.material.cell_models import Tentusscher, TentusscherEndo
import ansys.health.heart.settings.material.ep_material as ep_materials
from ansys.health.heart.settings.material.ep_material import (
    Active,
    ActiveBeam,
    EPMaterialModel,
    Insulator,
    Passive,
)

midcelldata = dict(
    [
        ("gas_constant", 8314.472),
        ("t", 310),
        ("faraday_constant", 96485.3415),
        ("cm", 0.185),
        ("vc", 0.016404),
        ("vsr", 0.001094),
        ("vss", 5.468e-05),
        ("pkna", 0.03),
        ("ko", 5.4),
        ("nao", 140.0),
        ("cao", 2.0),
        ("gk1", 5.405),
        ("gkr", 0.153),
        ("gna", 14.838),
        ("gbna", 0.0002),
        ("gcal", 3.98e-05),
        ("gbca", 0.000592),
        ("gpca", 0.1238),
        ("gpk", 0.0146),
        ("pnak", 2.724),
        ("km", 1.0),
        ("kmna", 40.0),
        ("knaca", 1000.0),
        ("ksat", 0.1),
        ("alpha", 2.5),
        ("gamma", 0.35),
        ("kmca", 1.38),
        ("kmnai", 87.5),
        ("kpca", 0.0005),
        ("k1", 0.15),
        ("k2", 0.045),
        ("k3", 0.06),
        ("k4", 0.005),
        ("ec", 1.5),
        ("maxsr", 2.5),
        ("minsr", 1.0),
        ("vrel", 0.102),
        ("vleak", 0.00036),
        ("vxfer", 0.0038),
        ("vmaxup", 0.006375),
        ("kup", 0.00025),
        ("bufc", 0.2),
        ("kbufc", 0.001),
        ("bufsr", 10.0),
        ("kbufsf", 0.3),
        ("bufss", 0.4),
        ("kbufss", 0.00025),
        ("gks", 0.098),
        ("gto", 0.294),
        ("v", -85.423),
        ("ki", 138.52),
        ("nai", 10.132),
        ("cai", 0.000153),
        ("cass", 0.00042),
        ("casr", 4.272),
        ("rpri", 0.8978),
        ("xr1", 0.0165),
        ("xr2", 0.473),
        ("xs", 0.0174),
        ("m", 0.00165),
        ("h", 0.749),
        ("j", 0.6788),
        ("d", 3.288e-05),
        ("f", 0.7026),
        ("f2", 0.9526),
        ("fcass", 0.9942),
        ("s", 0.999998),
        ("r", 2.347e-08),
    ]
)

endocelldata = dict(
    [
        ("gas_constant", 8314.472),
        ("t", 310),
        ("faraday_constant", 96485.3415),
        ("cm", 0.185),
        ("vc", 0.016404),
        ("vsr", 0.001094),
        ("vss", 5.468e-05),
        ("pkna", 0.03),
        ("ko", 5.4),
        ("nao", 140.0),
        ("cao", 2.0),
        ("gk1", 5.405),
        ("gkr", 0.153),
        ("gna", 14.838),
        ("gbna", 0.0002),
        ("gcal", 3.98e-05),
        ("gbca", 0.000592),
        ("gpca", 0.1238),
        ("gpk", 0.0146),
        ("pnak", 2.724),
        ("km", 1.0),
        ("kmna", 40.0),
        ("knaca", 1000.0),
        ("ksat", 0.1),
        ("alpha", 2.5),
        ("gamma", 0.35),
        ("kmca", 1.38),
        ("kmnai", 87.5),
        ("kpca", 0.0005),
        ("k1", 0.15),
        ("k2", 0.045),
        ("k3", 0.06),
        ("k4", 0.005),
        ("ec", 1.5),
        ("maxsr", 2.5),
        ("minsr", 1.0),
        ("vrel", 0.102),
        ("vleak", 0.00036),
        ("vxfer", 0.0038),
        ("vmaxup", 0.006375),
        ("kup", 0.00025),
        ("bufc", 0.2),
        ("kbufc", 0.001),
        ("bufsr", 10.0),
        ("kbufsf", 0.3),
        ("bufss", 0.4),
        ("kbufss", 0.00025),
        ("gks", 0.392),
        ("gto", 0.073),
        ("v", -86.709),
        ("ki", 138.4),
        ("nai", 10.355),
        ("cai", 0.00013),
        ("cass", 0.00036),
        ("casr", 3.715),
        ("rpri", 0.9068),
        ("xr1", 0.00448),
        ("xr2", 0.476),
        ("xs", 0.0087),
        ("m", 0.00155),
        ("h", 0.7573),
        ("j", 0.7225),
        ("d", 3.164e-05),
        ("f", 0.8009),
        ("f2", 0.9778),
        ("fcass", 0.9953),
        ("s", 0.3212),
        ("r", 2.235e-08),
    ]
)

epicelldata = dict(
    [
        ("gas_constant", 8314.472),
        ("t", 310),
        ("faraday_constant", 96485.3415),
        ("cm", 0.185),
        ("vc", 0.016404),
        ("vsr", 0.001094),
        ("vss", 5.468e-05),
        ("pkna", 0.03),
        ("ko", 5.4),
        ("nao", 140.0),
        ("cao", 2.0),
        ("gk1", 5.405),
        ("gkr", 0.153),
        ("gna", 14.838),
        ("gbna", 0.0002),
        ("gcal", 3.98e-05),
        ("gbca", 0.000592),
        ("gpca", 0.1238),
        ("gpk", 0.0146),
        ("pnak", 2.724),
        ("km", 1.0),
        ("kmna", 40.0),
        ("knaca", 1000.0),
        ("ksat", 0.1),
        ("alpha", 2.5),
        ("gamma", 0.35),
        ("kmca", 1.38),
        ("kmnai", 87.5),
        ("kpca", 0.0005),
        ("k1", 0.15),
        ("k2", 0.045),
        ("k3", 0.06),
        ("k4", 0.005),
        ("ec", 1.5),
        ("maxsr", 2.5),
        ("minsr", 1.0),
        ("vrel", 0.102),
        ("vleak", 0.00036),
        ("vxfer", 0.0038),
        ("vmaxup", 0.006375),
        ("kup", 0.00025),
        ("bufc", 0.2),
        ("kbufc", 0.001),
        ("bufsr", 10.0),
        ("kbufsf", 0.3),
        ("bufss", 0.4),
        ("kbufss", 0.00025),
        ("gks", 0.392),
        ("gto", 0.294),
        ("v", -85.23),
        ("ki", 136.89),
        ("nai", 8.604),
        ("cai", 0.000126),
        ("cass", 0.00036),
        ("casr", 3.64),
        ("rpri", 0.9073),
        ("xr1", 0.00621),
        ("xr2", 0.4712),
        ("xs", 0.0095),
        ("m", 0.00172),
        ("h", 0.7444),
        ("j", 0.7045),
        ("d", 3.373e-05),
        ("f", 0.7888),
        ("f2", 0.9755),
        ("fcass", 0.9953),
        ("s", 0.999998),
        ("r", 2.42e-08),
    ]
)


def test_cellmodel():
    tentusendo = cell_models.TentusscherEndo()
    tentusmid = cell_models.TentusscherMid()
    tentusepi = cell_models.TentusscherEpi()

    assert tentusendo.model_dump() == endocelldata
    assert tentusmid.model_dump() == midcelldata
    assert tentusepi.model_dump() == epicelldata


def test_active():
    active0 = ep_materials.Active(sigma_fiber=1)
    assert active0.sigma_sheet is not None
    active = ep_materials.Active(sigma_fiber=1, sigma_sheet=1, sigma_sheet_normal=1)
    assert active.sigma_sheet_normal == 1
    active_beam = ep_materials.ActiveBeam(sigma_fiber=1)
    assert active_beam.pmjres is not None


def test_passive():
    passive = ep_materials.Passive(sigma_fiber=1, sigma_sheet_normal=4)
    assert not hasattr(passive, "cell_model")
    assert passive.sigma_sheet_normal == 4


def test_insulator():
    insulator = ep_materials.Insulator()
    assert insulator.sigma_fiber == 0
    assert insulator.beta == 0
    assert insulator.cm == 0


class TestEPMaterialModel:
    """Test EPMaterialModel base class."""

    def test_ep_material_model_creation_with_defaults(self):
        """Test creating EPMaterialModel with default values."""
        model = EPMaterialModel()

        # Check default values from ep_defaults
        assert model.sigma_fiber is None  # Should be None by default
        assert model.sigma_sheet is None
        assert model.sigma_sheet_normal is None
        assert model.beta is not None  # Should get value from ep_defaults
        assert model.cm is not None  # Should get value from ep_defaults
        assert model.lambda_ is None

    def test_ep_material_model_optional_fields(self):
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

    def test_ep_material_model_custom_values(self):
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

    def test_ep_material_model_sheet_parameter_sync(self):
        """Test that sigma_sheet and sigma_sheet_normal sync properly."""
        # Test sigma_sheet set, sigma_sheet_normal None -> sync
        model = EPMaterialModel(sigma_sheet=1.0, sigma_sheet_normal=None)
        assert model.sigma_sheet == 1.0
        assert model.sigma_sheet_normal == 1.0  # Should be synced

        # Test sigma_sheet_normal set, sigma_sheet None -> sync
        model = EPMaterialModel(sigma_sheet=None, sigma_sheet_normal=2.0)
        assert model.sigma_sheet == 2.0  # Should be synced
        assert model.sigma_sheet_normal == 2.0

    def test_ep_material_model_both_sheet_parameters_set(self):
        """Test when both sheet parameters are explicitly set."""
        model = EPMaterialModel(sigma_sheet=1.0, sigma_sheet_normal=2.0)
        # Both should retain their original values
        assert model.sigma_sheet == 1.0
        assert model.sigma_sheet_normal == 2.0

    def test_ep_material_model_neither_sheet_parameter_set(self):
        """Test when neither sheet parameter is set."""
        model = EPMaterialModel(sigma_sheet=None, sigma_sheet_normal=None)
        # Both should remain None
        assert model.sigma_sheet is None
        assert model.sigma_sheet_normal is None

    def test_ep_material_model_invalid_types(self):
        """Test validation with invalid types."""
        with pytest.raises(ValidationError):
            EPMaterialModel(sigma_fiber="invalid")

        with pytest.raises(ValidationError):
            EPMaterialModel(beta="not_a_number")

    def test_ep_material_model_serialization(self):
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

    def test_insulator_defaults(self):
        """Test Insulator default values."""
        insulator = Insulator()

        assert insulator.sigma_fiber == 0.0
        assert insulator.cm == 0.0
        assert insulator.beta == 0.0

    def test_insulator_custom_values(self):
        """Test Insulator with custom values."""
        insulator = Insulator(sigma_fiber=0.1, cm=0.05, beta=100.0)

        assert insulator.sigma_fiber == 0.1
        assert insulator.cm == 0.05
        assert insulator.beta == 100.0

    def test_insulator_inheritance(self):
        """Test that Insulator doesn't inherit from EPMaterialModel."""
        insulator = Insulator()

        # Insulator is standalone BaseModel, not EPMaterialModel
        assert isinstance(insulator, Insulator)
        assert not isinstance(insulator, EPMaterialModel)

    def test_insulator_serialization(self):
        """Test Insulator serialization."""
        insulator = Insulator(sigma_fiber=0.1)

        data = insulator.model_dump()
        assert data["sigma_fiber"] == 0.1
        assert data["cm"] == 0.0
        assert data["beta"] == 0.0

    def test_insulator_validation(self):
        """Test Insulator field validation."""
        with pytest.raises(ValidationError):
            Insulator(sigma_fiber="invalid")


class TestActive:
    """Test Active EP material class."""

    def test_active_inheritance(self):
        """Test that Active inherits from EPMaterialModel."""
        active = Active()
        assert isinstance(active, Active)
        assert isinstance(active, EPMaterialModel)

    def test_active_default_creation(self):
        """Test Active creation with defaults."""
        active = Active()
        assert active.solver_type == "Monodomain"
        assert active.sigma_fiber == 0.5
        assert active.sigma_sheet == 0.1
        assert active.sigma_sheet_normal == 0.1
        assert active.beta == 140  # From ep_defaults
        assert active.cm == 0.01

        # Active-specific defaults
        assert isinstance(active.cell_model, Tentusscher)

    def test_active_cell_model_default_factory(self):
        """Test that cell_model uses default_factory properly."""
        active1 = Active()
        active2 = Active()

        # Should be separate instances
        assert active1.cell_model is not active2.cell_model
        assert isinstance(active1.cell_model, Tentusscher)
        assert isinstance(active2.cell_model, Tentusscher)

    def test_active_custom_cell_model(self):
        """Test Active with custom cell model."""
        custom_cell = TentusscherEndo(gto=0.1)
        active = Active(cell_model=custom_cell)

        assert active.cell_model is custom_cell
        assert isinstance(active.cell_model, TentusscherEndo)
        assert active.cell_model.gto == 0.1

    def test_active_solver_type_override(self):
        """Test Active with custom solver type."""
        active = Active(solver_type="Eikonal")
        assert active.solver_type == "Eikonal"

    def test_active_optional_sigma_fields(self):
        """Test Active with optional sigma fields."""
        # Test with some None, some set
        active = Active(sigma_fiber=1.2, sigma_sheet=None, sigma_sheet_normal=0.8)

        assert active.sigma_fiber == 1.2
        assert active.sigma_sheet == 0.8  # Should sync from sigma_sheet_normal
        assert active.sigma_sheet_normal == 0.8

    def test_active_serialization_with_cell_model(self):
        """Test Active serialization includes cell model."""
        active = Active(sigma_fiber=1.0)

        data = active.model_dump()
        assert "cell_model" in data
        assert isinstance(data["cell_model"], dict)
        assert "gas_constant" in data["cell_model"]  # From Tentusscher

    def test_active_deserialization_with_cell_model(self):
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

    def test_active_beam_inheritance(self):
        """Test that ActiveBeam inherits from Active."""
        beam = ActiveBeam()
        assert isinstance(beam, ActiveBeam)
        assert isinstance(beam, Active)
        assert isinstance(beam, EPMaterialModel)

    def test_active_beam_defaults(self):
        """Test ActiveBeam default values."""
        beam = ActiveBeam()

        # Should have beam-specific defaults
        assert beam.sigma_fiber is not None  # From ep_defaults.material["beam"]
        assert isinstance(beam.cell_model, Tentusscher)
        assert beam.pmjres is not None  # Beam-specific field

    def test_active_beam_cell_model_factory(self):
        """Test that ActiveBeam uses TentusscherEndo by default."""
        beam1 = ActiveBeam()
        beam2 = ActiveBeam()

        # Should be separate TentusscherEndo instances
        assert beam1.cell_model is not beam2.cell_model
        assert isinstance(beam1.cell_model, Tentusscher)
        assert isinstance(beam2.cell_model, Tentusscher)

    def test_active_beam_custom_values(self):
        """Test ActiveBeam with custom values."""
        custom_cell = TentusscherEndo(gas_constant=8100.0)
        beam = ActiveBeam(sigma_fiber=2.0, cell_model=custom_cell, pmjres=50.0)

        assert beam.sigma_fiber == 2.0
        assert beam.cell_model is custom_cell
        assert beam.pmjres == 50.0

    def test_active_beam_pmjres_field(self):
        """Test pmjres field specific to ActiveBeam."""
        beam = ActiveBeam(pmjres=100.0)
        assert beam.pmjres == 100.0

        # Regular Active shouldn't have this field
        active = Active()
        assert not hasattr(active, "pmjres")

    def test_active_beam_serialization(self):
        """Test ActiveBeam serialization includes all fields."""
        beam = ActiveBeam(pmjres=80.0)

        data = beam.model_dump()
        assert "pmjres" in data
        assert data["pmjres"] == 80.0
        assert "cell_model" in data
        assert isinstance(data["cell_model"], dict)


class TestPassive:
    """Test Passive EP material class."""

    def test_passive_inheritance(self):
        """Test that Passive inherits from EPMaterialModel."""
        passive = Passive()
        assert isinstance(passive, Passive)
        assert isinstance(passive, EPMaterialModel)

    def test_passive_defaults(self):
        """Test Passive default values."""
        passive = Passive()

        # Should have passive-specific defaults
        assert passive.sigma_fiber is not None  # From ep_defaults
        assert passive.sigma_sheet is None
        assert passive.sigma_sheet_normal is None

        # Inherited from EPMaterialModel
        assert passive.beta is not None
        assert passive.cm is not None

    def test_passive_no_cell_model(self):
        """Test that Passive doesn't have cell_model field."""
        passive = Passive()
        assert not hasattr(passive, "cell_model")

    def test_passive_custom_values(self):
        """Test Passive with custom values."""
        passive = Passive(sigma_fiber=1.8, sigma_sheet=0.9, sigma_sheet_normal=1.1)

        assert passive.sigma_fiber == 1.8
        assert passive.sigma_sheet == 0.9
        assert passive.sigma_sheet_normal == 1.1

    def test_passive_sheet_parameter_sync(self):
        """Test sheet parameter syncing in Passive."""
        passive = Passive(sigma_sheet=1.2)
        assert passive.sigma_sheet == 1.2
        assert passive.sigma_sheet_normal == 1.2  # Should sync


class TestFieldInteractions:
    """Test interactions between different field types."""

    def test_optional_vs_required_fields(self):
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

    def test_default_factory_independence(self):
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

    def test_none_handling_in_serialization(self):
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

    def test_field_validation_inheritance(self):
        """Test that field validation works properly in inheritance."""
        # Test validation in base class
        with pytest.raises(ValidationError):
            EPMaterialModel(sigma_fiber="invalid")

        # Test validation in derived classes
        with pytest.raises(ValidationError):
            Active(sigma_fiber="invalid")

        with pytest.raises(ValidationError):
            ActiveBeam(pmjres="invalid")

        with pytest.raises(ValidationError):
            Passive(sigma_sheet="invalid")


class TestDefaultsIntegration:
    """Test integration with ep_defaults module."""

    @pytest.mark.xfail(reason="ep_defaults not properly mocked.1")
    @patch("ansys.health.heart.settings.material._pyd_ep_material.ep_defaults")
    def test_defaults_integration(self, mock_defaults):
        """Test proper integration with ep_defaults."""
        # Set up mock defaults
        mock_defaults.analysis = {"solvertype": "TestSolver"}
        mock_defaults.material = {
            "myocardium": {
                "sigma_fiber": MagicMock(m=1.5),
                "beta": MagicMock(m=1500.0),
                "cm": MagicMock(m=0.15),
            },
            "beam": {"sigma": MagicMock(m=2.0), "pmjres": MagicMock(m=50.0)},
        }

        # Test various models use defaults correctly
        ep_model = EPMaterialModel()
        assert ep_model.beta == 1500.0
        assert ep_model.cm == 0.15

        active = Active()
        assert active.solver_type == "TestSolver"
        assert active.beta == 1500.0

        beam = ActiveBeam()
        assert beam.sigma_fiber == 2.0
        assert beam.pmjres == 50.0

        passive = Passive()
        assert passive.sigma_fiber == 1.5

    def test_field_override_defaults(self):
        """Test that explicit values override defaults."""
        # Create models with explicit values that should override defaults
        ep_model = EPMaterialModel(beta=3000.0, cm=0.3)
        active = Active(solver_type="Eikonal")
        beam = ActiveBeam(pmjres=100.0)
        passive = Passive(sigma_fiber=5.0)

        # Explicit values should be used instead of defaults
        assert ep_model.beta == 3000.0
        assert ep_model.cm == 0.3
        assert active.solver_type == "Eikonal"
        assert beam.pmjres == 100.0
        assert passive.sigma_fiber == 5.0


class TestComplexSerialization:
    """Test complex serialization scenarios."""

    def test_nested_model_serialization(self):
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

    def test_json_serialization_with_defaults(self):
        """Test JSON serialization preserves defaults and optionals."""
        beam = ActiveBeam()

        # JSON round trip
        json_str = beam.model_dump_json()
        restored = ActiveBeam.model_validate_json(json_str)

        # All fields should be preserved
        assert isinstance(restored.cell_model, Tentusscher)
        assert restored.pmjres == beam.pmjres
        assert restored.sigma_fiber == beam.sigma_fiber
        assert restored.beta == beam.beta
