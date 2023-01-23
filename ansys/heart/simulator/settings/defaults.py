"""Module contains default values for the simulations."""
from pint import Quantity
import copy as copy
from dataclasses import dataclass, asdict
import yaml, json
import pathlib


"""Generic analysis settings."""
analysis = {
    "end_time": Quantity(3000.0, "ms"),
    "dtmin": Quantity(10.0, "ms"),
    "dtmax": Quantity(10.0, "ms"),
    "dt_d3plot": Quantity(50.0, "ms"),
    "dt_icvout": Quantity(1.0, "ms"),
    "global_damping": Quantity(0.5, "1/ms"),
}

"""Material settings."""
material = {
    "myocardium": {
        "isotropic": {
            "rho": Quantity(0.001, "g/mm^3"),
            "k1": Quantity(0.0023599999999999997, "N/mm"),
            "k2": Quantity(1.75),
        },
        "anisotropic": {"k1": Quantity(0.00049, "N/mm"), "k2": Quantity(9.01, "dimensionless")},
        "active": {
            "actype": 1,
            "tmax": Quantity(0.125, "MPa"),
            "ca2ionm": Quantity(4.35, "umol/L"),
            "prefill": Quantity(1000.0, "ms"),
        },
    },
    "atrium": {
        "rho": Quantity(0.001, "g/mm^3"),
        "itype": -1,
        "mu1": 0.0349,
        "alpha1": 2,
        "Comment": "Should be equivalent with MAT_077_H",
    },
    "cap": {
        "rho": Quantity(0.001, "g/mm^3"),
        "nu": 0.499,
        "c10": 1.0,
        "thickness": Quantity(5.0, "mm"),
    },
}

"""Boundary condition settings."""
boundary_conditions = {
    "pericardium": {"penalty_function": [0.65, 25], "spring_stiffness": Quantity(0.05, "MPa")},
    "valve": {
        "biventricle": Quantity(0.002, "MPa"),
        "fourchamber": Quantity(0.02, "MPa"),
        "scale_factor": {"normal": 0.5, "radial": 1.0},
    },
    "end_diastolic_cavity_pressure": {
        "left_ventricle": Quantity(0.002, "MPa"),
        "right_ventricle": Quantity(0.0005333299999999999, "MPa"),
    },
}

"""System model parameters."""
system_model = {
    "name": "ConstantPreloadWindkesselAfterload",
    "left_ventricle": {
        "constants": {
            "Rv": Quantity(6.65e-06, "MPa/(mm^3/ms)"),
            "Ra": Quantity(1.729e-05, "MPa/(mm^3/ms)"),
            "Rp": Quantity(0.000766, "MPa/(mm^3/ms)"),
            "Ca": Quantity(6390000.0, "mm^3/MPa"),
            "Pven": Quantity(0.002, "MPa"),
        },
        "initial_value": {"part": Quantity(0.00931, "MPa")},
    },
    "right_ventricle": {
        "constants": {
            "Rv": Quantity(6.65e-06, "MPa/(mm^3/ms)"),
            "Ra": Quantity(6.0514999999999996e-06, "MPa/(mm^3/ms)"),
            "Rp": Quantity(9.575e-05, "MPa/(mm^3/ms)"),
            "Ca": Quantity(28755000.0, "mm^3/MPa"),
            "Pven": Quantity(0.0005333299999999999, "MPa"),
        },
        "initial_value": {"part": Quantity(0.0019950000000000002, "MPa")},
    },
}

parameters = {
    "analysis": analysis,
    "material": material,
    "boundary_conditions": boundary_conditions,
    "system_model": system_model,
}


def remove_units_in_dictionary(d: dict):
    """Method to replace Quantity with value in a nested dictionary, that is, removes units."""
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
    """Stores settings/parmaters for boundary conditions."""

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
        dictionary = remove_units_in_dictionary(dictionary)
        dictionary = {"Simulation Settings": dictionary}

        with open(filename, "w") as f:
            if filename.suffix == ".yml":
                yaml.dump(dictionary, f)
            elif filename.suffix == ".json":
                json.dump(dictionary, f, indent=4)


A = AnalysisSettings()
A.set_defaults(analysis)
M = MaterialSetting()
M.set_defaults(material)
BC = BoundaryConditions()
BC.set_defaults(boundary_conditions)
S = SystemModel()
S.set_defaults(system_model)

SETTINGS = SimulationSettings(
    analysis=A,
    material=M,
    boundary_conditions=BC,
    system=S,
)
