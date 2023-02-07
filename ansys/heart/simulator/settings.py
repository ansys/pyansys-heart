"""Module that stores settings."""
from dataclasses import asdict, dataclass
import json
import pathlib
import copy

import ansys.heart.simulator.settings.defaults as defaults
import yaml
from typing import List


from pint import Quantity, UnitRegistry




class Settings:
    """Generic settings class."""

    def set_defaults(self, defaults: dict):
        """Read default settings from dictionary."""
        for key, value in self.__dict__.items():
            # if isinstance(value, dict):
            #     # this is currently not working
            #     self.set_defaults(defaults)
            #     print(f"inner: {key}")
            if key in defaults.keys():
                setattr(self, key, defaults[key])


@dataclass
class AnalysisSettings(Settings):
    """Class for analysis settings."""

    end_time: float = 0
    """End time of the simulation."""
    dtmin: float = 0
    """Minimum time-step of simulation."""
    dtmax: float = 0
    """Maximum time-step of simulation."""
    dt_d3plot: float = 0
    """Time-step of d3plot export."""
    dt_icvout: float = 0
    """Time-step of icvout export."""
    global_damping: float = 0
    """Global damping constant."""


@dataclass
class MaterialSetting(Settings):
    """Class for storing material settings."""

    myocardium: dict = None
    """Myocardium material."""
    atrium: dict = None
    """Atrial material."""
    cap: dict = None
    """Cap material."""


@dataclass
class BoundaryConditions(Settings):
    """Stores settings/parameters for boundary conditions."""

    pericardium: dict = None
    """Parameters for pericardium spring b.c."""
    valve: dict = None
    """Parameters for valve spring b.c."""
    end_diastolic_cavity_pressure: dict = None
    """End-diastolic pressure."""


@dataclass
class SystemModel(Settings):
    """Stores settings/parameters for system model."""

    left_ventricle: dict = None
    """Parameters for left ventricle."""
    right_ventricle: dict = None
    """Parameters for right ventricle."""


@dataclass
class SimulationSettings:
    """Class for keeping track of settings."""

    analysis: AnalysisSettings = None
    """Generic analysis settings."""
    material: MaterialSetting = None
    """Material settings/configuration."""
    boundary_conditions: BoundaryConditions = None
    """Boundary condition specifications."""
    system: SystemModel = None
    """System model settings."""

    def save(self, filename: pathlib.Path):
        """Save simulation settings to disk.

        Parameters
        ----------
        filename : pathlib.Path
            Path to target .json or .yml file.
        """
        if not isinstance(filename, pathlib.Path):
            filename = pathlib.Path(filename)

        if filename.suffix not in [".yml", ".json"]:
            raise ValueError(f"Data format {filename.suffix} not supported")

        dictionary = asdict(self)

        # serialize Quantity objects
        dictionary = copy.deepcopy(dictionary)
        _serialize_quantity(dictionary)
        dictionary = {"Simulation Settings": dictionary}

        with open(filename, "w") as f:
            if filename.suffix == ".yml":
                yaml.dump(dictionary, f, sort_keys=False)
            elif filename.suffix == ".json":
                json.dump(dictionary, f, indent=4, sort_keys=False)

    def load(self, filename: pathlib.Path):
        """Load simulation parameters from disk.

        Parameters
        ----------
        filename : pathlib.Path
            Path to .json or .yml with settings.
        """
        if not isinstance(filename, pathlib.Path):
            filename = pathlib.Path(filename)
        ureg = UnitRegistry()
        # ureg()
        # loop over dictionary and try to convert to quantity.
        with open(filename, "r") as f:
            if filename.suffix == ".json":
                data = json.load(f)
            if filename.suffix == ".yml":
                data = yaml.load(f, Loader=yaml.SafeLoader)
        settings = data["Simulation Settings"]
        _deserialize_quantity(settings, ureg)

        # assign values to each respective attribute
        A = AnalysisSettings()
        A.set_defaults(settings["analysis"])
        M = MaterialSetting()
        M.set_defaults(settings["material"])
        BC = BoundaryConditions()
        BC.set_defaults(settings["boundary_conditions"])
        S = SystemModel()
        S.set_defaults(settings["system"])

        self.analysis = A
        self.material = M
        self.boundary_conditions = BC
        self.system = S

    def load_defaults(self):
        """Loads the default simulation settings."""

        # intialize parameters with defaults.
        A = AnalysisSettings()
        A.set_defaults(defaults.analysis)
        M = MaterialSetting()
        M.set_defaults(defaults.material)
        BC = BoundaryConditions()
        BC.set_defaults(defaults.boundary_conditions)
        S = SystemModel()
        S.set_defaults(defaults.system_model)

        self.analysis = A
        self.material = M
        self.boundary_conditions = BC
        self.system = S

def _remove_units_in_dictionary(d: dict):
    """Replace Quantity with value in a nested dictionary, that is, removes units."""
    for k, v in d.items():
        if isinstance(v, dict):
            _remove_units_in_dictionary(v)
        if isinstance(v, Quantity):
            d[k] = d[k].m
    return d

def _serialize_quantity(d: dict):
    """Serialize Quantity such that Quantity objects are replaced by <value> <units> string."""
    for k, v in d.items():
        if isinstance(v, dict):
            _serialize_quantity(v)
        if isinstance(v, Quantity):
            d[k] = str(d[k])
    return d


def _deserialize_quantity(d: dict, ureg: UnitRegistry):
    """Deserialize Quantity such that string <value> <units> is replaced by Quantity."""

    for k, v in d.items():
        if isinstance(v, dict):
            _deserialize_quantity(v, ureg)
        if isinstance(v, str):
            if isinstance(d[k], str):
                try:
                    float(d[k].split()[0])
                    q = ureg(d[k])
                except ValueError:
                    # failed to convert to quantity
                    continue
                d[k] = q
    return d


# some additional methods
def _get_dimensionality(d):
    dims = []
    for k, v in d.items():
        if isinstance(v, dict):
            dims += _get_dimensionality(v)
        elif isinstance(v, Quantity):
            dims.append(v.dimensionality)
    return dims


# get units
def _get_units(d):
    units = []
    for k, v in d.items():
        if isinstance(v, dict):
            units += _get_units(v)
        elif isinstance(v, Quantity):
            units.append(v.units)
    return units


settings = SimulationSettings()
settings.load_defaults()


# consistent_unit_system = MPa, mm, N, ms, g
#
# Length: mm
# Time: ms
# Mass: g
#
# Force: N [ Mass * Length / Time^2 ] --> [ kg*m / s^2]
# Pressure/Stress: MPa [ Mass / (Length * Time^2 )] --> [g / (mm * ms^2)]
#
# get dimensionality of consistent unit system
# _ureg = UnitRegistry()

# _unit_system = ["MPa", "mm", "N", "ms", "g"]

# get unit-system
_dimensions = _get_dimensionality(asdict(settings))
_units = _get_units(asdict(settings))

_units = list(set(_units))
_dimensions = list(set(_dimensions))
        