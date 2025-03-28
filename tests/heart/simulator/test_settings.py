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
import pytest

from ansys.heart.simulator.settings.defaults import fibers as fibers_defaults
from ansys.heart.simulator.settings.settings import (
    Analysis,
    Fibers,
    SimulationSettings,
    Stimulation,
    _get_consistent_units_str,
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
    "    material:\n"
    "      myocardium: null\n"
    "      passive: null\n"
    "      cap: null\n"
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
    "    material:\n"
    "      myocardium: null\n"
    "      atrium: null\n"
    "      cap: null\n"
    "      beam: null\n"
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
    stim = Stimulation(t_start=0, period=800, duration=20, amplitude=50)

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
    stim2 = Stimulation(node_ids=[1, 2, 3], t_start=10, period=100, duration=30, amplitude=40)
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
    settings = Fibers()
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
    """Check if defaults properly set."""
    settings = Fibers()
    settings.set_values(fibers_defaults.angles)
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

    assert default_settings.mechanics.material.myocardium["isotropic"]["rho"].m == pytest.approx(
        0.001, 1e-9
    )


def test_get_meca_material(default_settings):
    default_settings.mechanics.material.myocardium["isotropic"]["rho"] = Quantity(8000, "kg/m^3")
    default_settings.to_consistent_unit_system()

    # test default value
    m1 = default_settings.get_mechanical_material("anisotropic")
    assert m1.active.actype == 1

    m1 = default_settings.get_mechanical_material("anisotropic", ep_coupled=True)
    assert m1.active.actype == 3

    m2 = default_settings.get_mechanical_material("isotropic")
    assert m2.iso.mu1 == pytest.approx(0.1, 1e-9)
    # test modified value
    assert m1.rho == pytest.approx(0.008, 1e-9)


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
