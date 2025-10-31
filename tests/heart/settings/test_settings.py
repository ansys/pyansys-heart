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

"""Collection of methods to test the settings module."""

import os
import tempfile

import numpy as np
from pint import Quantity
from pydantic import ValidationError
import pytest

from ansys.health.heart.settings.defaults import fibers as fibers_defaults
from ansys.health.heart.settings.settings import (
    Analysis,
    FibersBRBM,
    SimulationSettings,
    Stimulation,
    ZeroPressure,
    _get_consistent_units_str,
    _windows_to_wsl_path,
)
from tests.heart.conftest import compare_string_with_file

REF_STRING_SETTINGS_YML_MECHANICS = (
    "Simulation Settings:\n"
    "  mechanics:\n"
    "    analysis:\n"
    "      end_time: 1 second\n"
    "      dtmin: 2 second\n"
    "      dtmax: 3 second\n"
    "      dt_d3plot: 4 second\n"
    "      dt_icvout: 5 millisecond\n"
    "      global_damping: 0.33 / second\n"
    "      stiffness_damping: 0.1 second\n"
    "    boundary_conditions:\n"
    "      robin: null\n"
    "      valve: null\n"
    "      end_diastolic_cavity_pressure: null\n"
    "    system:\n"
    "      name: ConstantPreloadWindkesselAfterload\n"
    "      left_ventricle: null\n"
    "      right_ventricle: null\n"
)

REF_STRING_SETTINGS_YML_EP = (
    "Simulation Settings:\n"
    "  electrophysiology:\n"
    "    analysis:\n"
    "      end_time: 1 second\n"
    "      dtmin: 2 second\n"
    "      dtmax: 3 second\n"
    "      dt_d3plot: 4 second\n"
    "      dt_icvout: 5 millisecond\n"
    "      global_damping: 0 / second\n"
    "      stiffness_damping: 0 second\n"
    "      solvertype: Monodomain\n"
    "    stimulation:\n"
    "      stimdefaults:\n"
    "        node_ids: null\n"
    "        t_start: 0 millisecond\n"
    "        period: 800 millisecond\n"
    "        duration: 20 millisecond\n"
    "        amplitude: 50 microfarad / millimeter ** 3\n"
    "    layers:\n"
    "      percent_endo: 0.17 dimensionless\n"
    "      percent_mid: 0.41 dimensionless\n"
    "    lambda_ratio: 0.2 dimensionless\n"
)


def test_settings_save_001():
    """Test saving of settings to disk."""
    settings = SimulationSettings(
        mechanics=True,
        electrophysiology=False,
        fiber=False,
        purkinje=False,
        stress_free=False,
    )

    # fill some dummy data
    settings.mechanics.analysis.end_time = Quantity(1, "s")
    settings.mechanics.analysis.dtmin = Quantity(2, "s")
    settings.mechanics.analysis.dtmax = Quantity(3, "s")
    settings.mechanics.analysis.dt_d3plot = Quantity(4, "s")
    settings.mechanics.analysis.dt_icvout = Quantity(5, "ms")
    settings.mechanics.analysis.global_damping = Quantity(0.33, "s**-1")
    settings.mechanics.analysis.stiffness_damping = Quantity(0.1, "s")

    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        file_path = os.path.join(tempdir, "settings.yml")
        settings.save(file_path)
        compare_string_with_file(REF_STRING_SETTINGS_YML_MECHANICS, file_path)

    pass


def test_settings_save_002():
    """Test saving of EP settings to disk."""

    settings = SimulationSettings(
        mechanics=False,
        electrophysiology=True,
        fiber=False,
        purkinje=False,
        stress_free=False,
    )
    stim = Stimulation(
        t_start=Quantity(0, "ms"),
        period=Quantity(800, "ms"),
        duration=Quantity(20, "ms"),
        amplitude=Quantity(50, "uF/mm^3"),
    )

    settings.electrophysiology.stimulation = {"stimdefaults": stim}
    # fill some dummy data
    settings.electrophysiology.analysis.end_time = Quantity(1, "s")
    settings.electrophysiology.analysis.dtmin = Quantity(2, "s")
    settings.electrophysiology.analysis.dtmax = Quantity(3, "s")
    settings.electrophysiology.analysis.dt_d3plot = Quantity(4, "s")
    settings.electrophysiology.analysis.dt_icvout = Quantity(5, "ms")

    # settings.electrophysiology.material.beam["sigma"] = Quantity(1, "mS")

    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        file_path = os.path.join(tempdir, "settings_ep.yml")
        settings.save(file_path)

        compare_string_with_file(REF_STRING_SETTINGS_YML_EP, file_path)
    settings.load_defaults()
    stim2 = Stimulation(
        node_ids=[1, 2, 3],
        t_start=Quantity(10, "ms"),
        period=Quantity(100, "ms"),
        duration=Quantity(30, "ms"),
        amplitude=Quantity(40, "uF/mm^3"),
    )
    settings.electrophysiology.stimulation["stim2"] = stim2
    stim: Stimulation = settings.electrophysiology.stimulation["stim2"]

    assert stim.amplitude.m == 40
    assert stim.duration.m == 30
    assert stim.t_start.m == 10
    pass


def test_settings_load():
    """Test loading of settings from file."""
    # write file-to-load from reference string
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        file_path = os.path.join(tempdir, "settings.yml")

        with open(file_path, "w") as f:
            f.write(REF_STRING_SETTINGS_YML_MECHANICS)

        # load settings
        settings = SimulationSettings()
        settings.load(file_path)

        # assert the settings are properly populated
        assert settings.mechanics.analysis.end_time == Quantity(1, "s")
        assert settings.mechanics.analysis.dtmin == Quantity(2, "s")
        assert settings.mechanics.analysis.dtmax == Quantity(3, "s")
        assert settings.mechanics.analysis.dt_d3plot == Quantity(4, "s")
        assert settings.mechanics.analysis.dt_icvout == Quantity(5, "ms")
        assert settings.mechanics.analysis.global_damping == Quantity(0.33, "s**-1")
        assert settings.mechanics.analysis.stiffness_damping == Quantity(0.1, "s")
    pass


def test_get_consistent_units():
    """Test conversion to consistent unit system."""
    q = Quantity(30, "s")
    assert _get_consistent_units_str(q.dimensionality) == "ms**1"
    q = Quantity(30, "m")
    assert _get_consistent_units_str(q.dimensionality) == "mm**1"
    q = Quantity(30, "kg")
    assert _get_consistent_units_str(q.dimensionality) == "g**1"
    q = Quantity(30, "Pa")
    assert _get_consistent_units_str(q.dimensionality) == "MPa"
    q = Quantity(30, "kN")
    assert _get_consistent_units_str(q.dimensionality) == "N"


def test_convert_units_001():
    """Test consistent unit conversion."""
    settings = Analysis()

    # NOTE: use settings.end_time attribute as dummy.

    # s --> ms (*1e3)
    settings.end_time = Quantity(50, "s")
    settings.to_consistent_unit_system()
    assert abs(settings.end_time.m - 50 * 1e3) < 1e-15

    # Pa --> MPa (*1e-6)
    settings.end_time = Quantity(50, "Pa")
    settings.to_consistent_unit_system()
    assert abs(settings.end_time.m - 50 * 1e-6) < 1e-15

    # mm/s --> mm/ms (*1e-3)
    settings.end_time = Quantity(50, "mm/s")
    settings.to_consistent_unit_system()
    assert abs(settings.end_time.m - 50 * 1e-3) < 1e-15

    # kg / m^4 / s --> g / mm^4 / ms (*1e-12)
    settings.end_time = Quantity(1e12, "kg / m^4 / s")
    settings.to_consistent_unit_system()
    assert abs(settings.end_time.m - 1.0) < 1e-15


def test_convert_units_002():
    """Test consistent unit conversion."""
    stim = Stimulation()
    stim.amplitude = Quantity(50, "uF/mm^3")
    settings = FibersBRBM()
    settings.alpha_endo = Quantity(10, "radian")

    # uF/mm^3 --> uF/mm^3
    stim.to_consistent_unit_system()
    assert stim.amplitude.m == 50

    # radian --> degree
    settings.to_consistent_unit_system()
    assert abs(settings.alpha_endo.m - 10 * 180 / np.pi) < 1e-15

    settings.alpha_endo = Quantity(10, "dimensionless")
    # dimensionless --> dimensionless
    settings.to_consistent_unit_system()
    assert abs(settings.alpha_endo.m - 10) < 1e-15


def test_settings_set_defaults():
    """Check if defaults properly set using Pydantic model initialization."""
    # Create Fibers instance with defaults applied directly
    fibers_data = fibers_defaults.angles
    settings = FibersBRBM(**fibers_data)
    assert settings.alpha_endo.m == -60


@pytest.fixture(autouse=True)
def default_settings():
    settings = SimulationSettings(
        mechanics=True,
        electrophysiology=False,
        fiber=False,
        purkinje=False,
        stress_free=False,
    )
    settings.load_defaults()
    return settings


@pytest.fixture(autouse=True)
def default_allsettings():
    settings = SimulationSettings(
        mechanics=True,
        electrophysiology=True,
        fiber=True,
        purkinje=True,
        stress_free=True,
    )
    settings.load_defaults()
    return settings


def test_load_defaults(default_settings):
    default_settings.to_consistent_unit_system()

    assert default_settings.mechanics.analysis.end_time.m == 800.0


def test_purkinje_settings(default_allsettings: SimulationSettings):
    # check default default values
    assert default_allsettings.purkinje.node_id_origin_left is None
    assert default_allsettings.purkinje.node_id_origin_right is None

    node_origin_left = np.empty(0, dtype=int)
    node_origin_right = np.empty(0, dtype=int)
    if default_allsettings.purkinje.node_id_origin_left is None:
        node_origin_left = 9
    if default_allsettings.purkinje.node_id_origin_right is None:
        node_origin_right = 10
    assert node_origin_left == 9
    assert node_origin_right == 10
    default_allsettings.purkinje.node_id_origin_left = 1
    node_origin_left = default_allsettings.purkinje.node_id_origin_left
    assert node_origin_left == 1
    default_allsettings.purkinje.node_id_origin_right = 2
    node_origin_right = default_allsettings.purkinje.node_id_origin_right
    assert node_origin_right == 2


@pytest.mark.xfail(condition=os.name == "posix", reason="Windows-specific test")
def test_windows_path_to_wsl_path():
    """Test conversion of Windows path to WSL path."""
    assert _windows_to_wsl_path("C:\\Program Files\\LS-DYNA") == "/mnt/c/Program Files/LS-DYNA"
    assert _windows_to_wsl_path("D:\\") == "/mnt/d/"
    assert _windows_to_wsl_path("D:\\dev") == "/mnt/d/dev"

    assert (
        _windows_to_wsl_path("\\\\wsl.localhost\\Ubuntu\\home\\user\\project")
        == "/home/user/project"
    )


# ZeroPressure test reference strings
REF_STRING_ZERO_PRESSURE_YML = (
    "Simulation Settings:\n"
    "  stress_free:\n"
    "    analysis:\n"
    "      end_time: 500 millisecond\n"
    "      dtmin: 5 millisecond\n"
    "      dtmax: 50 millisecond\n"
    "      dt_d3plot: 25 millisecond\n"
    "      dt_icvout: 10 millisecond\n"
    "      global_damping: 0.1 / second\n"
    "      stiffness_damping: 0.05 second\n"
    "      dt_nodout: 15 millisecond\n"
    "      max_iters: 5\n"
    "      method: 1\n"
    "      tolerance: 1.0\n"
)


def test_zero_pressure_serialization_yaml():
    """Test YAML serialization of ZeroPressure settings."""
    settings = SimulationSettings(
        mechanics=False,
        electrophysiology=False,
        fiber=False,
        purkinje=False,
        stress_free=True,
    )

    # Set custom values
    settings.stress_free.analysis.end_time = Quantity(500, "ms")
    settings.stress_free.analysis.dtmin = Quantity(5, "ms")
    settings.stress_free.analysis.dtmax = Quantity(50, "ms")
    settings.stress_free.analysis.dt_d3plot = Quantity(25, "ms")
    settings.stress_free.analysis.dt_icvout = Quantity(10, "ms")
    settings.stress_free.analysis.dt_nodout = Quantity(15, "ms")
    settings.stress_free.analysis.global_damping = Quantity(0.1, "1/s")
    settings.stress_free.analysis.stiffness_damping = Quantity(0.05, "s")
    settings.stress_free.analysis.max_iters = 5
    settings.stress_free.analysis.method = 1
    settings.stress_free.analysis.tolerance = 1.0

    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        file_path = os.path.join(tempdir, "zero_pressure.yml")
        settings.save(file_path)

        # Read file contents and compare manually
        with open(file_path, "r") as f:
            content = f.read()

        # Verify key content exists
        assert "stress_free:" in content
        assert "analysis:" in content
        assert "end_time: 500 millisecond" in content
        assert "max_iters: 5" in content
        assert "method: 1" in content
        assert "tolerance: 1.0" in content


def test_zero_pressure_deserialization_yaml():
    """Test YAML deserialization of ZeroPressure settings."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        file_path = os.path.join(tempdir, "zero_pressure_load.yml")

        # Write reference string to file
        with open(file_path, "w") as f:
            f.write(REF_STRING_ZERO_PRESSURE_YML)

        # Load settings
        settings = SimulationSettings(
            mechanics=False,
            electrophysiology=False,
            fiber=False,
            purkinje=False,
            stress_free=True,
        )
        settings.load(file_path)

        # Verify loaded values
        assert settings.stress_free.analysis.end_time == Quantity(500, "ms")
        assert settings.stress_free.analysis.dtmin == Quantity(5, "ms")
        assert settings.stress_free.analysis.dtmax == Quantity(50, "ms")
        assert settings.stress_free.analysis.dt_d3plot == Quantity(25, "ms")
        assert settings.stress_free.analysis.dt_icvout == Quantity(10, "ms")
        assert settings.stress_free.analysis.dt_nodout == Quantity(15, "ms")
        assert settings.stress_free.analysis.global_damping == Quantity(0.1, "1/s")
        assert settings.stress_free.analysis.stiffness_damping == Quantity(0.05, "s")
        assert settings.stress_free.analysis.max_iters == 5
        assert settings.stress_free.analysis.method == 1
        assert settings.stress_free.analysis.tolerance == 1.0


def test_zero_pressure_serialization_json():
    """Test JSON serialization of ZeroPressure settings."""
    settings = SimulationSettings(
        mechanics=False,
        electrophysiology=False,
        fiber=False,
        purkinje=False,
        stress_free=True,
    )

    # Load defaults
    settings.load_defaults()

    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        file_path = os.path.join(tempdir, "zero_pressure.json")
        settings.save(file_path)

        # Read and parse JSON
        with open(file_path, "r") as f:
            content = f.read()
            import json

            data = json.loads(content)

            # Verify structure
            assert "Simulation Settings" in data
            assert "stress_free" in data["Simulation Settings"]
            assert "analysis" in data["Simulation Settings"]["stress_free"]

            analysis = data["Simulation Settings"]["stress_free"]["analysis"]
            assert analysis["end_time"] == "1000.0 millisecond"
            assert analysis["max_iters"] == 3
            assert analysis["method"] == 2
            assert analysis["tolerance"] == 5.0


def test_zero_pressure_roundtrip():
    """Test roundtrip serialization/deserialization of ZeroPressure."""
    # Create settings with custom values
    original_settings = SimulationSettings(
        mechanics=False,
        electrophysiology=False,
        fiber=False,
        purkinje=False,
        stress_free=True,
    )

    # Set specific values
    original_settings.stress_free.analysis.end_time = Quantity(750, "ms")
    original_settings.stress_free.analysis.max_iters = 7
    original_settings.stress_free.analysis.method = 3
    original_settings.stress_free.analysis.tolerance = 2.5
    original_settings.stress_free.analysis.dt_nodout = Quantity(100, "ms")

    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        file_path = os.path.join(tempdir, "roundtrip.yml")

        # Save settings
        original_settings.save(file_path)

        # Load settings into new object
        loaded_settings = SimulationSettings(
            mechanics=False,
            electrophysiology=False,
            fiber=False,
            purkinje=False,
            stress_free=True,
        )
        loaded_settings.load(file_path)

        # Verify roundtrip consistency
        assert (
            loaded_settings.stress_free.analysis.end_time
            == original_settings.stress_free.analysis.end_time
        )
        assert (
            loaded_settings.stress_free.analysis.max_iters
            == original_settings.stress_free.analysis.max_iters
        )
        assert (
            loaded_settings.stress_free.analysis.method
            == original_settings.stress_free.analysis.method
        )
        assert (
            loaded_settings.stress_free.analysis.tolerance
            == original_settings.stress_free.analysis.tolerance
        )
        assert (
            loaded_settings.stress_free.analysis.dt_nodout
            == original_settings.stress_free.analysis.dt_nodout
        )


def test_zero_pressure_unit_conversion():
    """Test unit conversion for ZeroPressure settings."""
    zero_pressure = ZeroPressure()

    # Set values with different units
    zero_pressure.analysis.end_time = Quantity(2, "s")  # Will convert to ms
    zero_pressure.analysis.dtmin = Quantity(0.01, "s")  # Will convert to ms
    zero_pressure.analysis.global_damping = Quantity(2, "1/s")  # Should stay 1/s

    # Apply unit conversion
    zero_pressure.to_consistent_unit_system()

    # Verify conversions
    assert zero_pressure.analysis.end_time.magnitude == 2000.0
    assert str(zero_pressure.analysis.end_time.units) == "millisecond"
    assert zero_pressure.analysis.dtmin.magnitude == 10.0
    assert str(zero_pressure.analysis.dtmin.units) == "millisecond"
    assert zero_pressure.analysis.global_damping.magnitude == 0.002
    assert str(zero_pressure.analysis.global_damping.units) == "1 / millisecond"


def test_zero_pressure_validation():
    """Test Pydantic validation for ZeroPressure fields."""
    zero_pressure = ZeroPressure()

    # Test invalid max_iters (should be int)
    with pytest.raises(ValidationError):
        zero_pressure.analysis.max_iters = "invalid"

    # Test invalid tolerance (should be float)
    with pytest.raises(ValidationError):
        zero_pressure.analysis.tolerance = "invalid"

    # Test valid updates
    zero_pressure.analysis.max_iters = 15
    zero_pressure.analysis.tolerance = 3.14
    zero_pressure.analysis.method = 5

    assert zero_pressure.analysis.max_iters == 15
    assert zero_pressure.analysis.tolerance == 3.14
    assert zero_pressure.analysis.method == 5


def test_zero_pressure_defaults_loading():
    """Test loading default values for ZeroPressure from defaults module."""
    settings = SimulationSettings(
        mechanics=False,
        electrophysiology=False,
        fiber=False,
        purkinje=False,
        stress_free=True,
    )

    # Load defaults
    settings.load_defaults()

    # Verify default values from zeropressure defaults module
    assert settings.stress_free.analysis.end_time == Quantity(1000, "ms")
    assert settings.stress_free.analysis.dtmin == Quantity(10, "ms")
    assert settings.stress_free.analysis.dtmax == Quantity(100, "ms")
    assert settings.stress_free.analysis.dt_d3plot == Quantity(100, "ms")
    assert settings.stress_free.analysis.dt_nodout == Quantity(200, "ms")

    # Verify base class defaults are preserved
    assert settings.stress_free.analysis.max_iters == 3
    assert settings.stress_free.analysis.method == 2
    assert settings.stress_free.analysis.tolerance == 5.0
