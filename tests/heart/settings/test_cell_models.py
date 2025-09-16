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

"""Test module for _pyd_cell_models.py."""

from pydantic import ValidationError
import pytest

from ansys.health.heart.settings.material.cell_models import (
    Tentusscher,
    TentusscherEndo,
    TentusscherEpi,
    TentusscherMid,
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
    tentusendo = TentusscherEndo()
    tentusmid = TentusscherMid()
    tentusepi = TentusscherEpi()

    assert tentusendo.model_dump() == endocelldata
    assert tentusmid.model_dump() == midcelldata
    assert tentusepi.model_dump() == epicelldata


class TestTentusscher:
    """Test base Tentusscher cell model."""

    def test_tentusscher_default_creation(self):
        """Test creating Tentusscher with default values."""
        model = Tentusscher()

        # Test some key default values
        assert model.gas_constant == 8314.472
        assert model.t == 310
        assert model.faraday_constant == 96485.3415
        assert model.cm == 0.185
        assert model.ko == 5.4
        assert model.nao == 140.0
        assert model.cao == 2.0
        assert model.gks == 0.392
        assert model.gto == 0.294

    def test_tentusscher_custom_values(self):
        """Test creating Tentusscher with custom values."""
        model = Tentusscher(gas_constant=8000.0, t=300.0, ko=5.0, nao=130.0)

        assert model.gas_constant == 8000.0
        assert model.t == 300.0
        assert model.ko == 5.0
        assert model.nao == 130.0
        # Default values should remain unchanged
        assert model.faraday_constant == 96485.3415
        assert model.cm == 0.185

    def test_tentusscher_all_required_fields(self):
        """Test that all expected fields are present."""
        model = Tentusscher()

        # Check presence of key physiological parameters
        required_fields = [
            "gas_constant",
            "t",
            "faraday_constant",
            "cm",
            "vc",
            "vsr",
            "vss",
            "pkna",
            "ko",
            "nao",
            "cao",
            "gk1",
            "gkr",
            "gna",
            "gbna",
            "gcal",
            "gbca",
            "gpca",
            "gpk",
            "pnak",
            "km",
            "kmna",
            "knaca",
            "ksat",
            "alpha",
            "gamma",
            "kmca",
            "kmnai",
            "kpca",
            "k1",
            "k2",
            "k3",
            "k4",
            "ec",
            "maxsr",
            "minsr",
            "vrel",
            "vleak",
            "vxfer",
            "vmaxup",
            "kup",
            "bufc",
            "kbufc",
            "bufsr",
            "kbufsf",
            "bufss",
            "kbufss",
            "gks",
            "gto",
            "v",
            "ki",
            "nai",
            "cai",
            "cass",
            "casr",
            "rpri",
            "xr1",
            "xr2",
            "xs",
            "m",
            "h",
            "j",
            "d",
            "f",
            "f2",
            "fcass",
            "s",
            "r",
        ]

        for field in required_fields:
            assert hasattr(model, field), f"Missing field: {field}"
            assert getattr(model, field) is not None, f"Field {field} is None"

    def test_tentusscher_field_types(self):
        """Test that all fields have correct types."""
        model = Tentusscher()

        # All fields should be float
        for field_name, field_value in model.model_dump().items():
            assert isinstance(field_value, (float, int)), (
                f"Field {field_name} is not numeric: {type(field_value)}"
            )

    def test_tentusscher_validation_negative_values(self):
        """Test validation with negative values where appropriate."""
        # Some parameters like voltage can be negative
        model = Tentusscher(v=-90.0, h=0.5)
        assert model.v == -90.0
        assert model.h == 0.5

    def test_tentusscher_validation_zero_values(self):
        """Test validation with zero values."""
        model = Tentusscher(gbna=0.0, vleak=0.0)
        assert model.gbna == 0.0
        assert model.vleak == 0.0

    def test_tentusscher_scientific_notation(self):
        """Test handling of scientific notation values."""
        model = Tentusscher(d=3.373e-5, r=2.42e-8)
        assert model.d == 3.373e-5
        assert model.r == 2.42e-8

    def test_tentusscher_serialization(self):
        """Test serialization to dictionary."""
        model = Tentusscher(gas_constant=8000.0, t=305.0)
        data = model.model_dump()

        assert isinstance(data, dict)
        assert data["gas_constant"] == 8000.0
        assert data["t"] == 305.0
        assert "faraday_constant" in data

    def test_tentusscher_deserialization(self):
        """Test deserialization from dictionary."""
        data = {"gas_constant": 8000.0, "t": 305.0, "ko": 6.0}

        model = Tentusscher(**data)
        assert model.gas_constant == 8000.0
        assert model.t == 305.0
        assert model.ko == 6.0
        # Default values should be preserved
        assert model.faraday_constant == 96485.3415

    def test_tentusscher_json_serialization(self):
        """Test JSON serialization."""
        model = Tentusscher(gas_constant=8100.0)

        # Test JSON serialization
        json_str = model.model_dump_json()
        assert isinstance(json_str, str)

        # Test JSON deserialization
        model_restored = Tentusscher.model_validate_json(json_str)
        assert model_restored.gas_constant == 8100.0
        assert model_restored.faraday_constant == model.faraday_constant


class TestTentusscherEndo:
    """Test TentusscherEndo cell model."""

    def test_tentusscher_endo_inheritance(self):
        """Test that TentusscherEndo inherits from Tentusscher."""
        model = TentusscherEndo()
        assert isinstance(model, Tentusscher)
        assert isinstance(model, TentusscherEndo)

    def test_tentusscher_endo_overrides(self):
        """Test that TentusscherEndo overrides specific values."""
        base_model = Tentusscher()
        endo_model = TentusscherEndo()

        # Test specific overrides for endocardium
        assert endo_model.gks == 0.392  # Same as base
        assert endo_model.gto == 0.073  # Different from base (0.294)
        assert endo_model.v == -86.709  # Different from base (-85.23)
        assert endo_model.ki == 138.4  # Different from base (136.89)

        # Inherited values should be the same
        assert endo_model.gas_constant == base_model.gas_constant
        assert endo_model.faraday_constant == base_model.faraday_constant

    def test_tentusscher_endo_unique_parameters(self):
        """Test endocardium-specific parameter values."""
        model = TentusscherEndo()

        # Test key endocardium-specific values
        assert model.gto == 0.073
        assert model.v == -86.709
        assert model.ki == 138.4
        assert model.nai == 10.355
        assert model.s == 0.3212  # Very different from base (0.999998)

    def test_tentusscher_endo_custom_values(self):
        """Test TentusscherEndo with custom values."""
        model = TentusscherEndo(gto=0.1, v=-90.0)

        assert model.gto == 0.1
        assert model.v == -90.0
        # Other endo-specific defaults should remain
        assert model.ki == 138.4
        assert model.nai == 10.355


class TestTentusscherEpi:
    """Test TentusscherEpi cell model."""

    def test_tentusscher_epi_inheritance(self):
        """Test that TentusscherEpi inherits from Tentusscher."""
        model = TentusscherEpi()
        assert isinstance(model, Tentusscher)
        assert isinstance(model, TentusscherEpi)

    def test_tentusscher_epi_same_as_base(self):
        """Test that TentusscherEpi has same values as base Tentusscher."""
        base_model = Tentusscher()
        epi_model = TentusscherEpi()

        # TentusscherEpi should have identical values to base
        assert epi_model.gks == base_model.gks
        assert epi_model.gto == base_model.gto
        assert epi_model.v == base_model.v
        assert epi_model.ki == base_model.ki
        assert epi_model.nai == base_model.nai
        assert epi_model.s == base_model.s

    def test_tentusscher_epi_explicit_values(self):
        """Test TentusscherEpi explicit parameter values."""
        model = TentusscherEpi()

        # Test key epicardium values (same as base Tentusscher)
        assert model.gks == 0.392
        assert model.gto == 0.294
        assert model.v == -85.23
        assert model.ki == 136.89
        assert model.nai == 8.604
        assert model.s == 0.999998


class TestTentusscherMid:
    """Test TentusscherMid cell model."""

    def test_tentusscher_mid_inheritance(self):
        """Test that TentusscherMid inherits from Tentusscher."""
        model = TentusscherMid()
        assert isinstance(model, Tentusscher)
        assert isinstance(model, TentusscherMid)

    def test_tentusscher_mid_overrides(self):
        """Test that TentusscherMid overrides specific values."""
        mid_model = TentusscherMid()

        # Test specific overrides for mid-myocardium
        assert mid_model.gks == 0.098  # Very different from base (0.392)
        assert mid_model.gto == 0.294  # Same as base
        assert mid_model.v == -85.423  # Slightly different from base (-85.23)
        assert mid_model.ki == 138.52  # Different from base (136.89)

    def test_tentusscher_mid_unique_parameters(self):
        """Test mid-myocardium-specific parameter values."""
        model = TentusscherMid()

        # Test key mid-myocardium-specific values
        assert model.gks == 0.098
        assert model.v == -85.423
        assert model.ki == 138.52
        assert model.nai == 10.132
        assert model.cai == 0.000153
        assert model.xr1 == 0.0165


class TestModelComparisons:
    """Test comparisons between different Tentusscher variants."""

    def test_model_differences(self):
        """Test key differences between model variants."""
        endo = TentusscherEndo()
        epi = TentusscherEpi()
        mid = TentusscherMid()

        # gks values should be different
        assert endo.gks == 0.392
        assert epi.gks == 0.392
        assert mid.gks == 0.098  # Significantly lower

        # gto values
        assert endo.gto == 0.073  # Lower
        assert epi.gto == 0.294
        assert mid.gto == 0.294

        # s values (repolarization parameter)
        assert endo.s == 0.3212  # Much lower
        assert epi.s == 0.999998  # High
        assert mid.s == 0.999998  # High

    def test_model_serialization_differences(self):
        """Test that serialized models show expected differences."""
        endo = TentusscherEndo()
        epi = TentusscherEpi()
        mid = TentusscherMid()

        endo_data = endo.model_dump()
        epi_data = epi.model_dump()
        mid_data = mid.model_dump()

        # Key differentiating parameters
        assert endo_data["gto"] != epi_data["gto"]
        assert mid_data["gks"] != epi_data["gks"]
        assert endo_data["s"] != epi_data["s"]


class TestValidationAndErrors:
    """Test validation and error handling."""

    def test_invalid_field_types(self):
        """Test validation with invalid field types."""
        with pytest.raises(ValidationError):
            Tentusscher(gas_constant="invalid")

        with pytest.raises(ValidationError):
            TentusscherEndo(gks="not_a_number")

    def test_model_validation(self):
        """Test model validation with edge cases."""
        # Very large values should be accepted
        model = Tentusscher(gas_constant=1e10, faraday_constant=1e10)
        assert model.gas_constant == 1e10

        # Very small positive values
        model = Tentusscher(d=1e-20, r=1e-20)
        assert model.d == 1e-20
        assert model.r == 1e-20

    def test_field_assignment_after_creation(self):
        """Test field assignment after model creation."""
        # Pydantic models are immutable by default, but we can test validation
        with pytest.raises(ValidationError):
            # This should fail if we try to create with invalid type
            Tentusscher.model_validate({"gks": "invalid"})


class TestSerialization:
    """Test comprehensive serialization scenarios."""

    def test_round_trip_serialization(self):
        """Test complete serialization round trip for all models."""
        models = [
            Tentusscher(gas_constant=8200.0),
            TentusscherEndo(gto=0.08),
            TentusscherEpi(v=-84.0),
            TentusscherMid(gks=0.1),
        ]

        for model in models:
            # Test dict serialization
            data = model.model_dump()
            restored = type(model)(**data)

            assert type(restored) is type(model)
            assert restored.model_dump() == data

            # Test JSON serialization
            json_str = model.model_dump_json()
            restored_json = type(model).model_validate_json(json_str)

            assert type(restored_json) is type(model)
            assert restored_json.model_dump() == model.model_dump()

    def test_partial_serialization(self):
        """Test serialization with only some fields specified."""
        # Create model with only some custom values
        model = TentusscherEndo(gto=0.1, v=-88.0)

        data = model.model_dump()
        assert data["gto"] == 0.1
        assert data["v"] == -88.0
        assert data["gas_constant"] == 8314.472  # Default value preserved

    def test_json_schema_generation(self):
        """Test JSON schema generation."""
        schema = Tentusscher.model_json_schema()

        assert isinstance(schema, dict)
        assert "properties" in schema
        assert "gas_constant" in schema["properties"]
        assert "t" in schema["properties"]

        # Check that field types are correctly specified
        assert schema["properties"]["gas_constant"]["type"] == "number"
        assert schema["properties"]["t"]["type"] == "number"
