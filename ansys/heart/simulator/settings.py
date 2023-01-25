"""Module that stores settings."""
from dataclasses import asdict, dataclass
import json
import pathlib

import ansys.heart.simulator.settings.defaults as defaults
from pint import Quantity
import yaml


def remove_units_in_dictionary(d: dict):
    """Replace Quantity with value in a nested dictionary, that is, removes units."""
    for k, v in d.items():
        if isinstance(v, dict):
            remove_units_in_dictionary(v)
        if isinstance(v, Quantity):
            d[k] = d[k].m
    return d


class Settings:
    """Generic settings class."""

    def set_defaults(self, defaults: dict):
        """Read default settings."""
        for key, value in self.__dict__.items():
            if isinstance(value, dict):
                self.set_defaults(defaults)
            if key in defaults.keys():
                setattr(self, key, defaults[key])


@dataclass
class AnalysisSettings(Settings):
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

    analysis: AnalysisSettings
    """Generic analysis settings."""
    material: MaterialSetting
    """Material settings/configuration."""
    boundary_conditions: BoundaryConditions
    """Boundary condition specifications."""
    system: SystemModel
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
        # dictionary = remove_units_in_dictionary(dictionary)
        dictionary = {"Simulation Settings": dictionary}

        with open(filename, "w") as f:
            if filename.suffix == ".yml":
                yaml.dump(dictionary, f)
            elif filename.suffix == ".json":
                json.dump(dictionary, f, indent=4)


A = AnalysisSettings()
A.set_defaults(defaults.analysis)
M = MaterialSetting()
M.set_defaults(defaults.material)
BC = BoundaryConditions()
BC.set_defaults(defaults.boundary_conditions)
S = SystemModel()
S.set_defaults(defaults.system_model)

SETTINGS = SimulationSettings(
    analysis=A,
    material=M,
    boundary_conditions=BC,
    system=S,
)

# get dimensionality
def _get_dimensionality(d):
    dims = []
    for k, v in d.items():
        if isinstance(v, dict):
            dims += _get_dimensionality(v)
        elif isinstance(v, Quantity):
            print(k, ":", v)
            dims.append(v.dimensionality)
    return dims


# get units
def _get_units(d):
    units = []
    for k, v in d.items():
        if isinstance(v, dict):
            units += _get_units(v)
        elif isinstance(v, Quantity):
            print(k, ":", v)
            units.append(v.units)
    return units


# get unit-system

_values = _get_dimensionality(asdict(SETTINGS))
_units = _get_units(asdict(SETTINGS))

_all_units = list(set(_units))
