"""Module that stores settings."""

import copy
from dataclasses import asdict, dataclass, fields, replace
import json
import pathlib
from typing import List

import ansys.heart.simulator.settings.defaults as defaults
from ansys.heart.simulator.settings.defaults import mechanics as mech_defaults
from pint import Quantity, UnitRegistry
import yaml


def get_factorial(n: int):
    """recursively gets factorial"""
    if n < 2:
        return 1
    else:
        return n * get_factorial(n - 1)


class AttrDict(dict):
    """Dictionary subclass whose entries can be accessed by attributes
    (as well as normally).
    """

    def __init__(self, *args, **kwargs):
        def from_nested_dict(data):
            """Construct nested AttrDicts from nested dictionaries."""
            if not isinstance(data, dict):
                return data
            else:
                return AttrDict({key: from_nested_dict(data[key]) for key in data})

        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self
        for key in self.keys():
            self[key] = from_nested_dict(self[key])


class Settings:
    """Generic settings class."""

    def set_values(self, defaults: dict):
        """Read default settings from dictionary."""
        for key, value in self.__dict__.items():
            if key in defaults.keys():
                # set as AttrDict
                if isinstance(defaults[key], dict):
                    setattr(self, key, AttrDict(defaults[key]))
                else:
                    setattr(self, key, defaults[key])

    def serialize(self, remove_units: bool = False) -> dict:
        """Serialize the settings, that is formats the Quantity as str(<value> <unit>)"""
        dictionary = copy.deepcopy(asdict(self))
        _serialize_quantity(dictionary, remove_units)
        return dictionary

    def remove_units(self):
        def _remove_units(d):
            units = []
            if isinstance(d, Settings):
                d = d.__dict__
            for k, v in d.items():
                if isinstance(v, (dict, AttrDict, Settings)):
                    units += _remove_units(v)
                elif isinstance(v, Quantity):
                    # print(f"key: {k} | units {v.units}")
                    units.append(v.units)
                    d.update({k: v.m})
            return units

        removed_units = _remove_units(self)


@dataclass
class Analysis(Settings):
    """Class for analysis settings."""

    end_time: Quantity = 0
    """End time of the simulation."""
    dtmin: Quantity = 0
    """Minimum time-step of simulation."""
    dtmax: Quantity = 0
    """Maximum time-step of simulation."""
    dt_d3plot: Quantity = 0
    """Time-step of d3plot export."""
    dt_icvout: Quantity = 0
    """Time-step of icvout export."""
    global_damping: Quantity = 0
    """Global damping constant."""


@dataclass
class Material(Settings):
    """Class for storing material settings."""

    myocardium: AttrDict = None
    """Myocardium material."""
    atrium: AttrDict = None
    """Atrial material."""
    cap: AttrDict = None
    """Cap material."""


@dataclass
class BoundaryConditions(Settings):
    """Stores settings/parameters for boundary conditions."""

    pericardium: AttrDict = None
    """Parameters for pericardium spring b.c."""
    valve: AttrDict = None
    """Parameters for valve spring b.c."""
    end_diastolic_cavity_pressure: AttrDict = None
    """End-diastolic pressure."""


@dataclass
class SystemModel(Settings):
    """Stores settings/parameters for system model."""

    left_ventricle: AttrDict = None
    """Parameters for left ventricle."""
    right_ventricle: AttrDict = None
    """Parameters for right ventricle."""


@dataclass
class Mechanics(Settings):
    """Class for keeping track of settings."""

    analysis: Analysis = Analysis()
    """Generic analysis settings."""
    material: Material = Material()
    """Material settings/configuration."""
    boundary_conditions: BoundaryConditions = BoundaryConditions()
    """Boundary condition specifications."""
    system: SystemModel = SystemModel()
    """System model settings."""


@dataclass
class Electrophysiology(Settings):
    """Class for keeping track of electrophysiology settings."""

    analysis: Analysis = Analysis()
    """Generic analysis settings."""


@dataclass
class Fibers(Settings):
    """Class for keeping track of fiber settings."""

    # what parameters to expose?


@dataclass
class Purkinje(Settings):
    """Class for keeping track of purkinje settings."""


class SimulationSettings:
    """Class for keeping track of settings."""

    def __init__(
        self,
        mechanics: bool = True,
        electrophysiology: bool = True,
        fiber: bool = True,
        purkinje: bool = True,
    ) -> None:
        if mechanics:
            self.mechanics: Mechanics = Mechanics()
            """Settings for mechanical simulation."""

        if electrophysiology:
            self.electrophysiology: Electrophysiology = Electrophysiology()
            """Settings for electrophysiology simulation."""

        if fiber:
            self.fibers: Fibers = Fibers()
            """Settings for fiber generation."""

        if purkinje:
            self.purkinje: Purkinje = Purkinje()
            """Settings for Purkinje generation"""

        pass

    def save(self, filename: pathlib.Path, remove_units: bool = False):
        """Save simulation settings to disk.

        Parameters
        ----------
        filename : pathlib.Path
            Path to target .json or .yml file
        remove_units : bool, optional
            Flag indicating whether to remove units before writing, by default False
        """
        if not isinstance(filename, pathlib.Path):
            filename = pathlib.Path(filename)

        if filename.suffix not in [".yml", ".json"]:
            raise ValueError(f"Data format {filename.suffix} not supported")

        # serialize each of the settings.
        serialized_settings = {}
        for attribute_name in self.__dict__.keys():
            if not isinstance(getattr(self, attribute_name), Settings):
                continue
            else:
                setting: Settings = getattr(self, attribute_name)
                serialized_settings[attribute_name] = setting.serialize(remove_units=remove_units)

        serialized_settings = {"Simulation Settings": serialized_settings}

        with open(filename, "w") as f:
            if filename.suffix == ".yml":
                # note: this supresses writing of tags from AttrDict
                yaml.dump(json.loads(json.dumps(serialized_settings)), f, sort_keys=False)

            elif filename.suffix == ".json":
                json.dump(serialized_settings, f, indent=4, sort_keys=False)

    def load(self, filename: pathlib.Path):
        """Load simulation parameters from disk.

        Parameters
        ----------
        filename : pathlib.Path
            Path to .json or .yml with settings.
        """
        if not isinstance(filename, pathlib.Path):
            filename = pathlib.Path(filename)

        with open(filename, "r") as f:
            if filename.suffix == ".json":
                data = json.load(f)
            if filename.suffix == ".yml":
                data = yaml.load(f, Loader=yaml.SafeLoader)
        settings = data["Simulation Settings"]

        # unit registry to convert back to Quantity object
        ureg = UnitRegistry()

        try:
            attribute_name = "mechanics"
            _deserialize_quantity(settings[attribute_name], ureg)
            # assign values to each respective attribute
            A = Analysis()
            A.set_values(settings[attribute_name]["analysis"])
            M = Material()
            M.set_values(settings[attribute_name]["material"])
            BC = BoundaryConditions()
            BC.set_values(settings[attribute_name]["boundary_conditions"])
            S = SystemModel()
            S.set_values(settings[attribute_name]["system"])
            self.mechanics.analysis = A
            self.mechanics.material = M
            self.mechanics.boundary_conditions = BC
            self.mechanics.system = S
        except:
            print("Failed to load mechanics settings.")

    def load_defaults(self):
        """Load the default simulation settings."""
        # TODO move to Settings class
        if hasattr(self, "mechanics"):
            A = Analysis()
            A.set_values(mech_defaults.analysis)
            M = Material()
            M.set_values(mech_defaults.material)
            BC = BoundaryConditions()
            BC.set_values(mech_defaults.boundary_conditions)
            S = SystemModel()
            S.set_values(mech_defaults.system_model)

            self.mechanics.analysis = A
            self.mechanics.material = M
            self.mechanics.boundary_conditions = BC
            self.mechanics.system = S


def _remove_units_from_dictionary(d: dict):
    """Replace Quantity with value in a nested dictionary, that is, removes units."""
    for k, v in d.items():
        if isinstance(v, (dict, AttrDict)):
            _remove_units_from_dictionary(v)
        if isinstance(v, Quantity):
            d[k] = d[k].m
    return d


def _serialize_quantity(d: dict, remove_units: bool = False):
    """Serialize Quantity such that Quantity objects are replaced by <value> <units> string."""
    for k, v in d.items():
        # if isinstance(v, AttrDict):
        #     v = dict(v)  # cast to regular dict
        if isinstance(v, (dict, AttrDict)):
            _serialize_quantity(v, remove_units=remove_units)
        if isinstance(v, Quantity):
            if remove_units:
                d[k] = str(d[k].m)
            else:
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
        if isinstance(v, (dict, AttrDict)):
            dims += _get_dimensionality(v)
        elif isinstance(v, Quantity):
            dims.append(v.dimensionality)
    return dims


# get units
def _get_units(d):
    units = []
    for k, v in d.items():
        if isinstance(v, (dict, AttrDict)):
            units += _get_units(v)
        elif isinstance(v, Quantity):
            units.append(v.units)
    return units


settings = SimulationSettings()
settings.load_defaults()

settings.mechanics.remove_units()

settings.save("test.yml")
settings.save("test.json")

# settings.mechanics.analysis.end_time = Quantity(30, "ms")
# settings.save("test.yml")

settings.load("test.yml")

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
_dimensions = _get_dimensionality(asdict(settings.mechanics))
_units = _get_units(asdict(settings.mechanics))

_units = list(set(_units))
_dimensions = list(set(_dimensions))
