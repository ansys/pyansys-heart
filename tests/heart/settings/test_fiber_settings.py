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

"""Simple tests for fiber method conditional settings."""

import pytest

from ansys.health.heart.settings.settings import (
    FibersBRBM,
    FibersDRBM,
    SimulationSettings,
)


class TestFiberMethodSelection:
    """Test fiber method selection and type safety."""

    def test_lsdyna_method_creates_correct_type(self):
        """Test that LSDYNA method creates FibersBRBM instance."""
        settings = SimulationSettings(fiber_method="LSDYNA")

        # Check method tracking
        assert settings._fiber_method == "LSDYNA"

        # Check correct type created
        assert isinstance(settings.fibers, FibersBRBM)
        assert not isinstance(settings.fibers, FibersDRBM)

        # Check LSDYNA-specific attributes exist
        assert hasattr(settings.fibers, "beta_endo_septum")
        assert hasattr(settings.fibers, "beta_epi_septum")

        # Check D-RBM-specific attributes don't exist
        assert not hasattr(settings.fibers, "left_ventricle")
        assert not hasattr(settings.fibers, "right_ventricle")

    def test_drbm_method_creates_correct_type(self):
        """Test that D-RBM method creates FibersDRBM instance."""
        settings = SimulationSettings(fiber_method="D-RBM")

        # Check method tracking
        assert settings._fiber_method == "D-RBM"

        # Check correct type created
        assert isinstance(settings.fibers, FibersDRBM)
        assert not isinstance(settings.fibers, FibersBRBM)

        # Check D-RBM-specific attributes exist
        assert hasattr(settings.fibers, "left_ventricle")
        assert hasattr(settings.fibers, "right_ventricle")
        assert hasattr(settings.fibers, "septal_fraction")

        # Check LSDYNA-specific attributes don't exist
        assert not hasattr(settings.fibers, "beta_endo_septum")
        assert not hasattr(settings.fibers, "beta_epi_septum")

    def test_invalid_method_raises_error(self):
        """Test that invalid fiber method raises ValueError."""
        with pytest.raises(ValueError, match="Invalid method to compute the fiber orientation"):
            SimulationSettings(fiber_method="INVALID")

    def test_fiber_disabled_no_method(self):
        """Test that disabling fibers sets method to None."""
        settings = SimulationSettings(fiber=False)

        assert settings._fiber_method is None
        assert not hasattr(settings, "fibers")


class TestFiberStructures:
    """Test fiber data structures for both methods."""

    def test_lsdyna_fiber_structure(self):
        """Test LSDYNA fiber structure has expected fields."""
        settings = SimulationSettings(fiber_method="LSDYNA")
        lsdyna_fibers = settings._get_fiber_config_lsdyna()

        # Check all required LSDYNA fields exist
        required_fields = [
            "alpha_endo",
            "alpha_epi",
            "beta_endo",
            "beta_epi",
            "beta_endo_septum",
            "beta_epi_septum",
        ]
        for field in required_fields:
            assert hasattr(lsdyna_fibers, field), f"Missing field: {field}"

    def test_drbm_fiber_structure(self):
        """Test D-RBM fiber structure has expected fields."""
        settings = SimulationSettings(fiber_method="D-RBM")
        drbm_fibers = settings._get_fiber_config_drbm()

        # Check nested ventricle structures exist
        assert hasattr(drbm_fibers, "left_ventricle")
        assert hasattr(drbm_fibers, "right_ventricle")

        # Check ventricle field structure
        for ventricle in [drbm_fibers.left_ventricle, drbm_fibers.right_ventricle]:
            for field in ["alpha_endo", "alpha_epi", "beta_endo", "beta_epi"]:
                assert hasattr(ventricle, field), f"Missing ventricle field: {field}"

    def test_drbm_get_rotation_dict(self):
        """Test D-RBM get_rotation_dict method."""
        settings = SimulationSettings(fiber_method="D-RBM")
        drbm_fibers = settings._get_fiber_config_drbm()

        rotation_dict = drbm_fibers._get_rotation_dict()

        # Check expected keys
        expected_keys = {
            "alpha_left",
            "alpha_right",
            "alpha_ot",
            "beta_left",
            "beta_right",
            "beta_ot",
        }
        assert set(rotation_dict.keys()) == expected_keys

        # Check ventricle values are lists of floats
        for key in ["alpha_left", "alpha_right", "beta_left", "beta_right"]:
            assert isinstance(rotation_dict[key], list)
            assert len(rotation_dict[key]) == 2
            assert all(isinstance(val, (int, float)) for val in rotation_dict[key])

        # Outflow tract values can be None or lists
        for key in ["alpha_ot", "beta_ot"]:
            value = rotation_dict[key]
            assert value is None or (isinstance(value, list) and len(value) == 2)


class TestFiberMethodConsistency:
    """Test that fiber method remains consistent throughout object lifecycle."""

    def test_method_consistency_lsdyna(self):
        """Test LSDYNA method consistency."""
        settings = SimulationSettings(fiber_method="LSDYNA")

        # Initial state
        assert settings._fiber_method == "LSDYNA"
        assert isinstance(settings.fibers, FibersBRBM)

        # After accessing fibers
        lsdyna_fibers = settings._get_fiber_config_lsdyna()
        assert settings._fiber_method == "LSDYNA"
        assert isinstance(settings.fibers, FibersBRBM)
        assert lsdyna_fibers is settings.fibers

    def test_method_consistency_drbm(self):
        """Test D-RBM method consistency."""
        settings = SimulationSettings(fiber_method="D-RBM")

        # Initial state
        assert settings._fiber_method == "D-RBM"
        assert isinstance(settings.fibers, FibersDRBM)

        # After accessing fibers
        drbm_fibers = settings._get_fiber_config_drbm()
        assert settings._fiber_method == "D-RBM"
        assert isinstance(settings.fibers, FibersDRBM)
        assert drbm_fibers is settings.fibers


class TestFiberDefaultsLoading:
    """Test that load_defaults() correctly initializes fiber settings for both methods."""

    def test_load_defaults_lsdyna_fibers(self):
        """Test load_defaults() correctly initializes B-RBM fiber settings."""
        from pint import Quantity

        # Create settings with LSDYNA fiber method
        settings = SimulationSettings(fiber_method="LSDYNA")

        # Load defaults
        settings.load_defaults()

        # Verify fiber method is maintained
        assert settings._fiber_method == "LSDYNA"
        assert isinstance(settings.fibers, FibersBRBM)

        # Verify B-RBM specific defaults are loaded from the defaults module
        lsdyna_fibers = settings._get_fiber_config_lsdyna()

        # Verify the specific default values are correct
        assert lsdyna_fibers.alpha_endo == Quantity(-60, "degree")
        assert lsdyna_fibers.alpha_epi == Quantity(60, "degree")
        assert lsdyna_fibers.beta_endo == Quantity(-65, "degree")
        assert lsdyna_fibers.beta_epi == Quantity(25, "degree")
        assert lsdyna_fibers.beta_endo_septum == Quantity(-65, "degree")
        assert lsdyna_fibers.beta_epi_septum == Quantity(25, "degree")

    def test_load_defaults_drbm_fibers(self):
        """Test load_defaults() correctly initializes D-RBM fiber settings."""
        from pint import Quantity

        # Create settings with D-RBM fiber method
        settings = SimulationSettings(fiber_method="D-RBM")

        # Verify initial state with Pydantic defaults
        assert settings._fiber_method == "D-RBM"
        assert isinstance(settings.fibers, FibersDRBM)

        # Load defaults (should not change D-RBM values since they use Pydantic defaults)
        settings.load_defaults()

        # Verify fiber method is maintained
        assert settings._fiber_method == "D-RBM"
        assert isinstance(settings.fibers, FibersDRBM)

        # Verify D-RBM structure and default values from class definition
        drbm_fibers = settings._get_fiber_config_drbm()

        # Verify that the structure is correct
        assert hasattr(drbm_fibers, "left_ventricle")
        assert hasattr(drbm_fibers, "right_ventricle")
        assert hasattr(drbm_fibers, "septal_fraction")

        # Verify left ventricle default values (from class definition)
        assert drbm_fibers.left_ventricle.alpha_endo == Quantity(60, "degree")
        assert drbm_fibers.left_ventricle.alpha_epi == Quantity(-60, "degree")
        assert drbm_fibers.left_ventricle.beta_endo == Quantity(-20, "degree")
        assert drbm_fibers.left_ventricle.beta_epi == Quantity(20, "degree")

        # Verify right ventricle default values (from class definition)
        assert drbm_fibers.right_ventricle.alpha_endo == Quantity(-90, "degree")
        assert drbm_fibers.right_ventricle.alpha_epi == Quantity(25, "degree")
        assert drbm_fibers.right_ventricle.beta_endo == Quantity(0, "degree")
        assert drbm_fibers.right_ventricle.beta_epi == Quantity(20, "degree")

        # Verify septal fraction default
        assert drbm_fibers.septal_fraction == 2.0 / 3.0

        # Verify outflow tract values are None by default
        assert drbm_fibers.alpha_outflow_tract is None
        assert drbm_fibers.beta_outflow_tract is None

        # Verify get_rotation_dict functionality works after defaults loading
        rotation_dict = drbm_fibers._get_rotation_dict()
        assert isinstance(rotation_dict, dict)
        expected_keys = {
            "alpha_left",
            "alpha_right",
            "alpha_ot",
            "beta_left",
            "beta_right",
            "beta_ot",
        }
        assert set(rotation_dict.keys()) == expected_keys

        # Verify correct values in rotation dict
        assert rotation_dict["alpha_left"] == [60.0, -60.0]
        assert rotation_dict["alpha_right"] == [-90.0, 25.0]
        assert rotation_dict["beta_left"] == [-20.0, 20.0]
        assert rotation_dict["beta_right"] == [0.0, 20.0]
        assert rotation_dict["alpha_ot"] is None
        assert rotation_dict["beta_ot"] is None

    def test_load_defaults_different_fiber_methods(self):
        """Test load_defaults() works correctly for both fiber methods independently."""
        from pint import Quantity

        # Test B-RBM first
        settings_lsdyna = SimulationSettings(fiber_method="LSDYNA")
        settings_lsdyna.load_defaults()

        assert settings_lsdyna._fiber_method == "LSDYNA"
        assert isinstance(settings_lsdyna.fibers, FibersBRBM)
        lsdyna_fibers = settings_lsdyna._get_fiber_config_lsdyna()
        assert lsdyna_fibers.alpha_endo == Quantity(-60, "degree")

        # Test D-RBM separately
        settings_drbm = SimulationSettings(fiber_method="D-RBM")
        settings_drbm.load_defaults()

        assert settings_drbm._fiber_method == "D-RBM"
        assert isinstance(settings_drbm.fibers, FibersDRBM)
        drbm_fibers = settings_drbm._get_fiber_config_drbm()
        assert hasattr(drbm_fibers, "left_ventricle")
        assert hasattr(drbm_fibers, "right_ventricle")

        # Verify they are independent and have different types
        assert settings_lsdyna._fiber_method != settings_drbm._fiber_method
        assert not isinstance(settings_lsdyna.fibers, type(settings_drbm.fibers))
