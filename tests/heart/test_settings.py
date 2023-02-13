"""Collection of methods to test the settings module."""

from ansys.heart.simulator.settings.settings1 import settings, _get_consistent_units_str
from pint import Quantity
from ansys.heart.simulator.settings.settings1 import SimulationSettings, Analysis
from .conftest import compare_string_with_file, get_workdir
import os


REF_STRING_SETTINGS_YML = (
    "Simulation Settings:\n"
    "  mechanics:\n"
    "    analysis:\n"
    "      end_time: 1 second\n"
    "      dtmin: 2 second\n"
    "      dtmax: 3 second\n"
    "      dt_d3plot: 4 second\n"
    "      dt_icvout: 5 millisecond\n"
    "      global_damping: 0.33 / second\n"
    "    material:\n"
    "      myocardium: null\n"
    "      atrium: null\n"
    "      cap: null\n"
    "    boundary_conditions:\n"
    "      pericardium: null\n"
    "      valve: null\n"
    "      end_diastolic_cavity_pressure: null\n"
    "    system:\n"
    "      left_ventricle: null\n"
    "      right_ventricle: null\n"
)


def test_settings_save():
    """Test saving of settings to disk."""
    settings = SimulationSettings(
        mechanics=True, electrophysiology=False, fiber=False, purkinje=False
    )
    # fill some dummy data
    settings.mechanics.analysis.end_time = Quantity(1, "s")
    settings.mechanics.analysis.dtmin = Quantity(2, "s")
    settings.mechanics.analysis.dtmax = Quantity(3, "s")
    settings.mechanics.analysis.dt_d3plot = Quantity(4, "s")
    settings.mechanics.analysis.dt_icvout = Quantity(5, "ms")
    settings.mechanics.analysis.global_damping = Quantity(0.33, "s**-1")

    file_path = os.path.join(get_workdir(), "settings.yml")
    settings.save(file_path)

    compare_string_with_file(REF_STRING_SETTINGS_YML, file_path)

    os.remove(file_path)
    pass


def test_settings_load():
    """Test loading of settings from file."""
    # write file-to-load from reference string
    file_path = os.path.join(get_workdir(), "settings.yml")
    with open(file_path, "w") as f:
        f.write(REF_STRING_SETTINGS_YML)

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

    os.remove(file_path)

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


def test_convert_units():
    """Test to consistent unit conversion."""
    settings = Analysis()

    # s --> ms
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
