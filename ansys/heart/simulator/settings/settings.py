"""Module that defines some classes for settings."""

import copy
from dataclasses import asdict, dataclass
import json
import pathlib

from ansys.heart import LOG as LOGGER
from ansys.heart.simulator.settings.defaults import mechanics as mech_defaults
from ansys.heart.simulator.settings.defaults import zeropressure as zero_pressure_defaults
from pint import Quantity, UnitRegistry
import yaml


class AttrDict(dict):
    """Dict subclass whose entries can be accessed by attributes as well as normally."""

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

    def __repr__(self):
        """Represent object in as dictionary in YAML style."""
        d = self.serialize()
        d = {self.__class__.__name__: d}
        return yaml.dump(json.loads(json.dumps(d)), sort_keys=False)

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
        """Serialize the settings, that is formats the Quantity as str(<value> <unit>)."""
        dictionary = copy.deepcopy(asdict(self))
        _serialize_quantity(dictionary, remove_units)
        return dictionary

    def to_consistent_unit_system(self):
        """Convert units to consistent unit system.

        Note
        ----
        Currently the only supported unit system is ["MPa", "mm", "N", "ms", "g"]
        For instance:
        Quantity(10, "mm/s") --> Quantity(0.01, "mm/ms")
        """

        def _to_consitent_units(d):
            if isinstance(d, Settings):
                d = d.__dict__
            for k, v in d.items():
                if isinstance(v, (dict, AttrDict, Settings)):
                    _to_consitent_units(v)
                elif isinstance(v, Quantity):
                    # print(f"key: {k} | units {v.units}")
                    if "[substance]" in list(v.dimensionality):
                        LOGGER.warning("Not converting [substance] / [length]^3")
                        continue
                    d.update({k: v.to(_get_consistent_units_str(v.dimensionality))})
            return

        _to_consitent_units(self)
        return

    def _remove_units(self):
        """Remove all units from Quantity objects."""

        def __remove_units(d):
            units = []
            if isinstance(d, Settings):
                d = d.__dict__
            for k, v in d.items():
                if isinstance(v, (dict, AttrDict, Settings)):
                    units += __remove_units(v)
                elif isinstance(v, Quantity):
                    # LOGGER.debug(f"key: {k} | units {v.units}")
                    units.append(v.units)
                    d.update({k: v.m})
            return units

        removed_units = __remove_units(self)
        return removed_units


@dataclass(repr=False)
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


@dataclass(repr=False)
class Material(Settings):
    """Class for storing material settings."""

    myocardium: AttrDict = None
    """Myocardium material."""
    atrium: AttrDict = None
    """Atrial material."""
    cap: AttrDict = None
    """Cap material."""


@dataclass(repr=False)
class BoundaryConditions(Settings):
    """Stores settings/parameters for boundary conditions."""

    pericardium: AttrDict = None
    """Parameters for pericardium spring b.c."""
    valve: AttrDict = None
    """Parameters for valve spring b.c."""
    end_diastolic_cavity_pressure: AttrDict = None
    """End-diastolic pressure."""


@dataclass(repr=False)
class SystemModel(Settings):
    """Stores settings/parameters for system model."""

    name: str = "ConstantPreloadWindkesselAfterload"
    """Name of the system model."""

    left_ventricle: AttrDict = None
    """Parameters for left ventricle."""
    right_ventricle: AttrDict = None
    """Parameters for right ventricle."""


@dataclass(repr=False)
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


@dataclass(repr=False)
class AnalysisZeroPressure(Analysis):
    """Class for keeping track of zero pressure analysis settings."""

    dt_nodout: Quantity = 0
    """Time interval of nodeout export."""

    max_iters: int = 3
    """Maximum iterations for stress-free-configuration algorithm."""
    method: int = 2
    """Method to use."""
    tolerance: float = 5.0
    """Tolerance to use for iterative algorithm."""


@dataclass(repr=False)
class ZeroPressure(Settings):
    """Class for keeping track of settings for stress-free-configuration computation."""

    analysis: AnalysisZeroPressure = AnalysisZeroPressure()
    """Generic analysis settings."""


@dataclass(repr=False)
class Electrophysiology(Settings):
    """Class for keeping track of electrophysiology settings."""

    analysis: Analysis = Analysis()
    """Generic analysis settings."""


@dataclass(repr=False)
class Fibers(Settings):
    """Class for keeping track of fiber settings."""

    # what parameters to expose?


@dataclass(repr=False)
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
        stress_free: bool = True,
    ) -> None:
        """Initialize Simulation Settings.

        Parameters
        ----------
        mechanics : bool, optional
            Flag indicating whether to add settings for mechanics, by default True
        electrophysiology : bool, optional
            Flag indicating whether to add settings for electrophysiology, by default True
        fiber : bool, optional
            Flag indicating whether to add settings for fiber generation, by default True
        purkinje : bool, optional
            Flag indicating whether to add settings for purkinje generation, by default True
        stress_free : bool, optional
            Flag indicating whether to add settings for the stress free
            configuration computation, by default True

        Example
        -------
        Instantiate settings and load defaults

        >>> from ansys.heart.simulator.settings.settings import SimulationSettings
        >>> settings = SimulationSettings()
        >>> settings.load_defaults()
        >>> print(settings)
        SimulationSettings
          mechanics
          electrophysiology
          fibers
          purkinje

        >>> print(settings.mechanics.analysis)
        Analysis:
          end_time: 3000.0 millisecond
          dtmin: 10.0 millisecond
          dtmax: 10.0 millisecond
          dt_d3plot: 50.0 millisecond
          dt_icvout: 1.0 millisecond
          global_damping: 0.5 / millisecond

        """
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
            """Settings for Purkinje generation."""

        if stress_free:
            self.stress_free: ZeroPressure = ZeroPressure()
            """Settings for stress free configuration simulation."""

        pass

    def __repr__(self):
        """Represent object as list of relevant attribute names."""
        repr_str = "\n  ".join(
            [attr for attr in self.__dict__ if isinstance(getattr(self, attr), Settings)]
        )
        repr_str = self.__class__.__name__ + "\n  " + repr_str
        return repr_str

    def save(self, filename: pathlib.Path, remove_units: bool = False):
        """Save simulation settings to disk.

        Parameters
        ----------
        filename : pathlib.Path
            Path to target .json or .yml file
        remove_units : bool, optional
            Flag indicating whether to remove units before writing, by default False

        Example
        -------
        Create examples settings with default values.

        >>> from ansys.heart.simulator.settings.settings import SimulationSettings
        >>> settings = SimulationSettings()
        >>> settings.load_defaults()
        >>> settings.save("my_settings.yml")

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
                # NOTE: this suppress writing of tags from AttrDict
                yaml.dump(json.loads(json.dumps(serialized_settings)), f, sort_keys=False)

            elif filename.suffix == ".json":
                json.dump(serialized_settings, f, indent=4, sort_keys=False)

    def load(self, filename: pathlib.Path):
        """Load simulation settings.

        Parameters
        ----------
        filename : pathlib.Path
            Path to yaml or json file.

        Example
        -------
        Create examples settings with default values.

        >>> from ansys.heart.simulator.settings.settings import SimulationSettings
        >>> settings = SimulationSettings()
        >>> settings.load_defaults()
        >>> settings.save("my_settings.yml")

        Load settings in second SimulationSettings object.

        >>> settings1 = SimulationSettings()
        >>> settings1.load("my_settings.yml")
        >>> assert settings.mechanics.analysis == settings1.mechanics.analysis
        True

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

            attribute_name = "stress_free"
            _deserialize_quantity(settings[attribute_name], ureg)
            A = AnalysisZeroPressure()
            A.set_values(settings[attribute_name]["analysis"])
            self.stress_free.analysis = A

        except:
            LOGGER.error("Failed to load mechanics settings.")

    def load_defaults(self):
        """Load the default simulation settings.

        Example
        -------
        Create examples settings with default values.

        Load module
        >>> from ansys.heart.simulator.settings.settings import SimulationSettings

        Instantiate settings object.

        >>> settings = SimulationSettings()
        >>> settings.load_defaults()
        >>> settings.mechanics.analysis
        Analysis:
          end_time: 3000.0 millisecond
          dtmin: 10.0 millisecond
          dtmax: 10.0 millisecond
          dt_d3plot: 50.0 millisecond
          dt_icvout: 1.0 millisecond
          global_damping: 0.5 / millisecond

        """
        # TODO move to Settings class
        for attr in self.__dict__:
            if isinstance(getattr(self, attr), Mechanics):
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

            if isinstance(getattr(self, attr), ZeroPressure):
                A = AnalysisZeroPressure()
                A.set_values(zero_pressure_defaults.analysis)
                self.stress_free.analysis = A

            elif isinstance(getattr(self, attr), (Electrophysiology, Fibers, Purkinje)):
                LOGGER.warning(
                    "Reading EP, Fiber, ZeroPressure, and Purkinje settings not yet supported."
                )

    def to_consistent_unit_system(self):
        """Convert all settings to consistent unit-system ["MPa", "mm", "N", "ms", "g"].

        Example
        -------
        Convert to the consistent unit system ["MPa", "mm", "N", "ms", "g"].

        Import necessary modules
        >>> from ansys.heart.simulator.settings.settings import SimulationSettings
        >>> from pint import Quantity

        Instantiate settings
        >>> settings = SimulationSettings()
        >>> settings.mechanics.analysis.end_time = Quantity(1, "s")
        >>> settings.to_consistent_unit_system()
        >>> settings.mechanics.analysis.end_time
        <Quantity(1000.0, 'millisecond')>

        """
        attributes = [
            getattr(self, attr)
            for attr in self.__dict__
            if isinstance(getattr(self, attr), Settings)
        ]

        for attr in attributes:
            if isinstance(attr, Settings):
                attr.to_consistent_unit_system()
        return


def _remove_units_from_dictionary(d: dict):
    """Replace Quantity with value in a nested dictionary (removes units)."""
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
    """Deserialize string such that "<value> <units>" is replaced by Quantity(value, units)."""
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


# desired consistent unit system is:
# ["MPa", "mm", "N", "ms", "g"]
# Time: ms
# Length: mm
# Mass: g
# Pressure: MPa
# Force: N
# base_quantitiy / unit mapping

_base_quantity_unit_mapper = {
    "[time]": "ms",
    "[length]": "mm",
    "[mass]": "g",
    "[substance]": "umol",
}
# these are derived quantities:
_derived = [
    [Quantity(30, "MPa").dimensionality, Quantity(30, "N").dimensionality],
    ["MPa", "N"],
]


def _get_consistent_units_str(dimensions: set):
    """Get consistent units formatted as string."""
    if dimensions in _derived[0]:
        _to_units = _derived[1][_derived[0].index(dimensions)]
        return _to_units

    _to_units = []
    for quantity in dimensions:
        _to_units.append(
            "{:s}**{:d}".format(_base_quantity_unit_mapper[quantity], dimensions[quantity])
        )
    return "*".join(_to_units)


if __name__ == "__main__":
    LOGGER.debug("protected")
