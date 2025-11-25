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

"""Test Quantity serialization and validation with Pydantic v2."""

import json
from pathlib import Path
import tempfile

from pint import Quantity
import pytest
import yaml

from ansys.health.heart.settings.settings import (
    Analysis,
    Mechanics,
    SimulationSettings,
    Stimulation,
)


class TestQuantitySerializationPydantic:
    """Test Pydantic v2 based Quantity serialization and validation."""

    def test_quantity_field_serializer_json(self) -> None:
        """Test that Quantity fields are properly serialized to JSON strings."""
        analysis = Analysis(
            end_time=Quantity(1000, "ms"),
            dtmin=Quantity(0.1, "ms"),
            dtmax=Quantity(10, "ms"),
            global_damping=Quantity(0.5, "1/s"),
        )

        # Test JSON serialization
        json_data = analysis.model_dump_json()
        parsed_data = json.loads(json_data)

        # Verify Quantity objects are serialized as strings with units
        assert parsed_data["end_time"] == "1000 millisecond"
        assert parsed_data["dtmin"] == "0.1 millisecond"
        assert parsed_data["dtmax"] == "10 millisecond"
        assert parsed_data["global_damping"] == "0.5 / second"

    def test_quantity_field_serializer_dict(self) -> None:
        """Test that Quantity fields are properly serialized to dict format."""
        analysis = Analysis(
            end_time=Quantity(1000, "ms"),
            dtmin=Quantity(0.1, "ms"),
            global_damping=Quantity(0.5, "1/s"),
        )

        # Test dict serialization for YAML
        dict_data = analysis.model_dump()

        # Verify Quantity objects are serialized as strings
        assert dict_data["end_time"] == Quantity(1000, "ms")
        assert dict_data["dtmin"] == Quantity(0.1, "ms")
        assert dict_data["global_damping"] == Quantity(0.5, "1/s")

    def test_quantity_validator_from_string(self) -> None:
        """Test that Quantity validators properly deserialize from strings."""
        # Test creation from string representation
        analysis_data = {
            "end_time": "1000.0 millisecond",
            "dtmin": "0.1 millisecond",
            "dtmax": "10.0 millisecond",
            "global_damping": "0.5 / second",
        }

        analysis = Analysis(**analysis_data)

        # Verify proper Quantity objects are created
        assert isinstance(analysis.end_time, Quantity)
        assert analysis.end_time.magnitude == 1000.0
        assert str(analysis.end_time.units) == "millisecond"

        assert isinstance(analysis.global_damping, Quantity)
        assert analysis.global_damping.magnitude == 0.5
        assert str(analysis.global_damping.units) == "1 / second"

    def test_quantity_validator_from_string_invalid_value(self) -> None:
        """Test that the validator raises an error when a non-quantity string is received."""
        with pytest.raises(ValueError):
            Analysis(end_time="invalid unit string")

    def test_quantity_validator_from_quantity(self) -> None:
        """Test that Quantity validators pass through existing Quantity objects."""
        # Test creation from Quantity objects directly
        analysis = Analysis(
            end_time=Quantity(1000, "ms"),
            dtmin=Quantity(0.1, "ms"),
            global_damping=Quantity(0.5, "1/s"),
        )

        # Verify Quantity objects are preserved
        assert isinstance(analysis.end_time, Quantity)
        assert analysis.end_time.magnitude == 1000.0
        assert str(analysis.end_time.units) == "millisecond"

    @pytest.mark.xfail(reason="Default units not yet implemented in Pydantic validators")
    def test_quantity_validator_from_numeric_with_default_units(self) -> None:
        """Test that numeric values get default units when specified."""
        # Create stimulation with numeric values - should get default units
        stim_data = {
            "t_start": 10.0,  # Should become Quantity(10.0, "ms")
            "period": 800.0,  # Should become Quantity(800.0, "ms")
            "duration": 2.0,  # Should become Quantity(2.0, "ms")
        }

        stim = Stimulation(**stim_data)

        # Verify proper conversion to Quantity with default units
        assert isinstance(stim.t_start, Quantity)
        assert stim.t_start.magnitude == 10.0
        assert str(stim.t_start.units) == "millisecond"

    @pytest.mark.xfail(reason="Targets units not yet implemented in Pydantic validators")
    def test_quantity_validation_error_invalid_units(self) -> None:
        """Test that invalid unit strings raise proper validation errors."""
        with pytest.raises(ValueError, match="Unable to parse quantity"):
            Analysis(end_time="invalid_unit_string")

    @pytest.mark.xfail(reason="Targets dimensions not yet implemented in Pydantic validators")
    def test_quantity_validation_error_incompatible_dimensions(self) -> None:
        """Test validation of incompatible dimensions."""
        with pytest.raises(ValueError, match="incompatible dimensions"):
            # Trying to assign length unit to time field
            Analysis(end_time="100.0 meter")

    def test_nested_quantity_serialization(self) -> None:
        """Test serialization of nested models with Quantity fields."""
        mechanics = Mechanics()
        mechanics.analysis.end_time = Quantity(1000, "ms")
        mechanics.analysis.dtmin = Quantity(0.1, "ms")

        # Test JSON serialization of nested model
        json_data = mechanics.model_dump_json()
        parsed_data = json.loads(json_data)

        # Verify nested Quantity serialization
        assert parsed_data["analysis"]["end_time"] == "1000 millisecond"
        assert parsed_data["analysis"]["dtmin"] == "0.1 millisecond"

    def test_stimulation_node_ids_validation(self) -> None:
        """Test that Stimulation node_ids field is properly validated."""
        # Test with list of integers
        stim = Stimulation(node_ids=[1, 2, 3])
        assert stim.node_ids == [1, 2, 3]

        # Test with None
        stim_none = Stimulation(node_ids=None)
        assert stim_none.node_ids is None

        # Test with invalid type
        with pytest.raises(ValueError):
            Stimulation(node_ids="invalid_node_ids")

    def test_complex_nested_serialization_deserialization(self) -> None:
        """Test round-trip serialization/deserialization of complex nested structures."""
        # Create complex settings
        settings = SimulationSettings(
            mechanics=True,
            electrophysiology=True,
            fiber=False,
            purkinje=False,
            stress_free=False,
        )

        # Set some values
        settings.mechanics.analysis.end_time = Quantity(1000, "ms")
        settings.electrophysiology.analysis.end_time = Quantity(800, "ms")

        stim = Stimulation(
            node_ids=[1, 2, 3],
            t_start=Quantity(10, "ms"),
            amplitude=Quantity(50, "uF/mm^3"),
        )
        settings.electrophysiology.stimulation = {"apex": stim}

        # Serialize to JSON
        json_data = settings.mechanics.model_dump_json()

        # Deserialize back
        mechanics_dict = json.loads(json_data)
        reconstructed_mechanics = Mechanics(**mechanics_dict)

        # Verify values are preserved
        assert reconstructed_mechanics.analysis.end_time == Quantity(1000, "ms")
        assert isinstance(reconstructed_mechanics.analysis.end_time, Quantity)

    def test_yaml_round_trip_serialization(self) -> None:
        """Test YAML serialization/deserialization round trip."""
        analysis = Analysis(
            end_time=Quantity(1000, "ms"),
            dtmin=Quantity(0.1, "ms"),
            global_damping=Quantity(0.5, "1/s"),
        )

        # Serialize to dict for YAML
        data_dict = analysis.model_dump(mode="json")

        # Convert to YAML and back
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yml", delete=False) as f:
            yaml.dump(data_dict, f)
            temp_path = f.name

        try:
            with open(temp_path, "r") as f:
                loaded_dict = yaml.load(f, Loader=yaml.SafeLoader)

            # Reconstruct object
            reconstructed = Analysis(**loaded_dict)

            # Verify values match
            assert reconstructed.end_time == analysis.end_time
            assert reconstructed.dtmin == analysis.dtmin
            assert reconstructed.global_damping == analysis.global_damping

        finally:
            Path(temp_path).unlink()

    def test_unit_conversion_during_validation(self) -> None:
        """Test that units are properly converted during validation."""
        # Create analysis with different time units
        analysis = Analysis(
            end_time="1.0 second",  # Should be converted to milliseconds
            dtmin="0.1 second",
        )

        # Convert to consistent unit system
        analysis.to_consistent_unit_system()

        # Verify conversion to consistent units
        assert analysis.end_time.magnitude == 1000.0
        assert str(analysis.end_time.units) == "millisecond"
        assert analysis.dtmin.magnitude == 100.0
        assert str(analysis.dtmin.units) == "millisecond"

    def test_serialize_consistency(self) -> None:
        """Test the serialize method for consistency with modern serialization."""
        analysis = Analysis(
            end_time=Quantity(1000, "ms"),
            dtmin=Quantity(0.1, "ms"),
            global_damping=Quantity(0.5, "1/s"),
        )

        # Test modern serialization method
        data_modern = analysis.model_dump(mode="json")

        # Results should be identical (both convert to strings)
        assert data_modern["end_time"] == "1000 millisecond"
        assert data_modern["global_damping"] == "0.5 / second"

    def test_dimensionless_quantity_handling(self) -> None:
        """Test handling of dimensionless quantities."""
        from ansys.health.heart.settings.settings import Electrophysiology

        ep = Electrophysiology()
        ep.lambda_ratio = Quantity(0.2, "dimensionless")

        # Test serialization
        data = ep.model_dump()
        assert data["lambda_ratio"] == Quantity(0.2, "dimensionless")

        # Test deserialization
        reconstructed = Electrophysiology(**data)
        assert reconstructed.lambda_ratio == Quantity(0.2, "dimensionless")

    def test_none_quantity_fields(self) -> None:
        """Test handling of None values in optional Quantity fields."""
        # Create model with None values (if any fields allow it)
        stim = Stimulation(node_ids=None)  # node_ids can be None

        # Serialize and verify None is preserved
        data = stim.model_dump(mode="json")
        assert data["node_ids"] is None

        # Deserialize and verify None is preserved
        reconstructed = Stimulation(**data)
        assert reconstructed.node_ids is None

    def test_quantity_field_validation_edge_cases(self) -> None:
        """Test edge cases in Quantity field validation."""
        # Test zero values
        analysis = Analysis(end_time="0.0 millisecond")
        assert analysis.end_time.magnitude == 0.0

        # Test negative values where appropriate
        analysis = Analysis(global_damping="-0.1 / second")
        assert analysis.global_damping.magnitude == -0.1

        # Test very small values
        analysis = Analysis(dtmin="1e-6 millisecond")
        assert analysis.dtmin.magnitude == 1e-6

    def test_pydantic_validation_assignment(self) -> None:
        """Test that assignment validation works with Quantity fields."""
        analysis = Analysis()

        # Test assignment of string
        analysis.end_time = "500.0 millisecond"
        assert isinstance(analysis.end_time, Quantity)
        assert analysis.end_time.magnitude == 500.0

        # Test assignment of Quantity
        analysis.end_time = Quantity(1000, "ms")
        assert analysis.end_time.magnitude == 1000.0

        # Test invalid assignment
        with pytest.raises(ValueError):
            analysis.end_time = "invalid_quantity"
