"""Module that defines some classes for settings."""

import copy
from dataclasses import asdict, dataclass, field
import json
import os
import pathlib
from typing import List, Literal

from ansys.heart import LOG as LOGGER
from ansys.heart.simulator.settings.defaults import electrophysiology as ep_defaults
from ansys.heart.simulator.settings.defaults import fibers as fibers_defaults
from ansys.heart.simulator.settings.defaults import mechanics as mech_defaults
from ansys.heart.simulator.settings.defaults import purkinje as purkinje_defaults
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
                elif isinstance(v, Quantity) and not (v.unitless):
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

    end_time: Quantity = Quantity(0, "s")
    """End time of the simulation."""
    dtmin: Quantity = Quantity(0, "s")
    """Minimum time-step of simulation."""
    dtmax: Quantity = Quantity(0, "s")
    """Maximum time-step of simulation."""
    dt_d3plot: Quantity = Quantity(0, "s")
    """Time-step of d3plot export."""
    dt_icvout: Quantity = Quantity(0, "s")
    """Time-step of icvout export."""
    global_damping: Quantity = Quantity(0, "1/s")
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
class EpMaterial(Settings):
    """Class for storing ep material settings."""

    myocardium: AttrDict = None
    """Myocardium material."""
    atrium: AttrDict = None
    """Atrial material."""
    cap: AttrDict = None
    """Cap material."""
    beam: AttrDict = None
    """beam material."""
    # TODO consider 'other', e.g passive conductor, soft tissue...?


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

    analysis: Analysis = field(default_factory=lambda: Analysis())
    """Generic analysis settings."""
    material: Material = field(default_factory=lambda: Material())
    """Material settings/configuration."""
    boundary_conditions: BoundaryConditions = field(default_factory=lambda: BoundaryConditions())
    """Boundary condition specifications."""
    system: SystemModel = field(default_factory=lambda: SystemModel())
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

    analysis: AnalysisZeroPressure = field(default_factory=lambda: AnalysisZeroPressure())
    """Generic analysis settings."""


@dataclass(repr=False)
class Electrophysiology(Settings):
    """Class for keeping track of electrophysiology settings."""

    material: EpMaterial = field(default_factory=lambda: EpMaterial())
    """Material settings/configuration."""
    analysis: Analysis = field(default_factory=lambda: Analysis())
    """Generic analysis settings."""


@dataclass(repr=False)
class Fibers(Settings):
    """Class for keeping track of fiber settings."""

    alpha_endo: Quantity = 0
    "Helical angle in endocardium."
    alpha_epi: Quantity = 0
    "Helical angle in epicardium."
    beta_endo: Quantity = 0
    "Angle to the outward transmural axis of the heart in endocardium."
    beta_epi: Quantity = 0
    "Angle to the outward transmural axis of the heart in epicardium."
    beta_endo_septum: Quantity = 0
    "Angle to the outward transmural axis of the heart in left septum."
    beta_epi_septum: Quantity = 0
    "Angle to the outward transmural axis of the heart in right septum."


@dataclass(repr=False)
class Purkinje(Settings):
    """Class for keeping track of purkinje settings."""

    edgelen: Quantity = 0
    """Edge length."""
    ngen: Quantity = 0
    """Number of generations."""
    nbrinit: Quantity = 0
    """Number of beams from origin point."""
    nsplit: Quantity = 0
    """Number of splits at each leaf."""
    pmjtype: Quantity = 0
    """Purkinje muscle junction type."""
    pmjradius: Quantity = 0
    """Purkinje muscle junction radius."""


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

            if isinstance(getattr(self, attr), Electrophysiology):
                A = Analysis()
                A.set_values(ep_defaults.analysis)
                M = EpMaterial()
                M.set_values(ep_defaults.material)
                self.electrophysiology.analysis = A
                self.electrophysiology.material = M
                # TODO add stim params, monodomain/bidomain/eikonal,cellmodel
            # TODO add settings for purkinje  fibers and epmecha
            if isinstance(getattr(self, attr), Fibers):
                self.fibers.set_values(fibers_defaults.angles)
            if isinstance(getattr(self, attr), Purkinje):
                self.purkinje.set_values(purkinje_defaults.build)

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
    "[current]": "mA",
}
# these are derived quantities:
_derived = [
    [
        Quantity(30, "MPa").dimensionality,
        Quantity(30, "N").dimensionality,
        Quantity(30, "mS/mm").dimensionality,
        Quantity(30, "uF/mm^2").dimensionality,
        Quantity(30, "1/mS").dimensionality,
        Quantity(30, "degree").dimensionality,
    ],
    ["MPa", "N", "mS/mm", "uF/mm^2", "1/mS", "degree"],
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


class DynaSettings:
    """Class for collecting, managing and validating LS-DYNA settings."""

    def __init__(
        self,
        lsdyna_path: pathlib.Path,
        dynatype: Literal["smp", "intelmpi", "platformmpi", "msmpi"] = "intelmpi",
        num_cpus: int = 1,
        platform: Literal["windows", "wsl", "linux"] = "windows",
        dyna_options: str = "",
        mpi_options: str = "",
    ):
        """Initialize Dyna settings.

        Parameters
        ----------
        lsdyna_path : Path
            Path to LS-DYNA
        dynatype : Literal[&quot;smp&quot;, &quot;intelmpi&quot;, &quot;platformmpi&quot;]
            Type of LS-DYNA executable. Shared Memory Parallel or Massively Parallel Processing
        num_cpus : int, optional
            Number of CPU's requested, by default 1
        platform : Literal[&quot;windows&quot;, &quot;wsl&quot;, &quot;linux&quot;], optional
            Platform, by default "windows"
        dyna_options : str, optional
            Additional command line options, by default ''
        mpi_options : str, optional
            Additional mpi options, by default ''
        """
        self.lsdyna_path: pathlib.Path = lsdyna_path
        """Path to LS-DYNA executable."""
        self.dynatype: str = dynatype
        """Type of dyna executable."""
        self.num_cpus: int = num_cpus
        """Number of CPU's requested."""
        self.platform: str = platform
        """Platform on which dyna is executed."""

        self.dyna_options = dyna_options
        """Additional command line options for dyna."""

        if dynatype in ["intelmpi", "platformmpi", "msmpi"]:
            self.mpi_options = mpi_options
            """additional mpi options."""
        elif dynatype == "smp":
            self.mpi_options = ""

        return

    def get_commands(self, path_to_input: pathlib.Path) -> List[str]:
        """Get command line arguments from the defined settings.

        Parameters
        ----------
        path_to_input : pathlib.Path
            Path to the LS-DYNA input file.

        Returns
        -------
        List[str]
            List of strings of each of the commands.
        """
        import subprocess

        lsdyna_path = self.lsdyna_path

        if self.platform == "windows" or self.platform == "linux":
            if self.dynatype in ["intelmpi", "platformmpi"]:
                commands = [
                    "mpirun",
                    self.mpi_options,
                    "-np",
                    str(self.num_cpus),
                    lsdyna_path,
                    "i=" + path_to_input,
                    self.dyna_options,
                ]
            elif self.dynatype in ["smp"]:
                commands = [
                    lsdyna_path,
                    "i=" + path_to_input,
                    "ncpu=" + str(self.num_cpus),
                    self.dyna_options,
                ]
        if self.platform == "windows" and self.dynatype == "msmpi":
            commands = [
                "mpiexec",
                self.mpi_options,
                "-np",
                str(self.num_cpus),
                lsdyna_path,
                "i=" + path_to_input,
                self.dyna_options,
            ]

        elif self.platform == "wsl":
            path_to_input_wsl = (
                subprocess.run(
                    ["wsl", "wslpath", os.path.basename(path_to_input)], capture_output=1
                )
                .stdout.decode()
                .strip()
            )
            # redefines LS-DYNA path.
            lsdyna_path = (
                subprocess.run(
                    ["wsl", "wslpath", str(lsdyna_path).replace("\\", "/")], capture_output=1
                )
                .stdout.decode()
                .strip()
            )

            if self.dynatype in ["intelmpi", "platformmpi", "msmpi"]:
                commands = [
                    "mpirun",
                    self.mpi_options,
                    "-np",
                    str(self.num_cpus),
                    lsdyna_path,
                    "i=" + path_to_input_wsl,
                    self.dyna_options,
                ]
            elif self.dynatype in ["smp"]:
                commands = [
                    lsdyna_path,
                    "i=" + path_to_input_wsl,
                    "ncpu=" + str(self.num_cpus),
                    self.dyna_options,
                ]

            path_to_run_script = os.path.join(pathlib.Path(path_to_input).parent, "run_lsdyna.sh")
            with open(path_to_run_script, "w", newline="\n") as f:
                f.write("#!/usr/bin/env sh\n")
                f.write("echo start lsdyna in wsl...\n")
                f.write(" ".join([i.strip() for i in commands]))

            commands = ["powershell", "-Command", "wsl", "-e", "bash", "-lic", "./run_lsdyna.sh"]

        # remove empty strings from commands
        commands = [c for c in commands if c != ""]

        # expand any environment variables if any
        commands = [os.path.expandvars(c) for c in commands]

        return commands

    def _set_env_variables(self):
        r"""Try to set environment variables for MPI run using Ansys installation root directories.

        Notes
        -----
        This is based on lsdynaintelvar.bat and lsdynamsvar.bat in
        ANSYS Inc\v231\ansys\bin\winx64\lsprepost48\LS-Run 1.0
        and requires you to install the MPI libraries with the Ansys installer.
        """
        # get latest installed Ansys version
        env_variables = list(os.environ.keys())
        ansys_env_roots = sorted([k for k in env_variables if "AWP_ROOT" in k])
        ansys_env_roots.reverse()

        if self.platform == "windows":
            platform = "winx64"
        elif self.platform == "linux":
            platform = "linx64"

        if self.dynatype in "intelmpi":
            INTEL_CMP_REV = "2019.5.281"
            INTEL_MKL_REV = "2020.0.166"
            INTEL_MPI_REV = "2018.3.210"

            if os.getenv("MPI_ROOT"):
                LOGGER.warning(
                    "MPI_ROOT already defined. Not trying to automatically set MPI env variables."
                )
                return

            for root in ansys_env_roots:
                mpi_root = os.path.join(
                    os.getenv(root), "commonfiles", "MPI", "Intel", INTEL_MPI_REV, platform
                )
                mpi_path = os.path.join(mpi_root, "bin")
                mkl_path = os.path.join(os.getenv(root), "tp", "IntelMKL", INTEL_MKL_REV, platform)
                cmp_path = os.path.join(
                    os.getenv(root), "tp", "IntelCompiler", INTEL_CMP_REV, platform
                )
                for p in [mpi_root, mpi_path, mkl_path, cmp_path]:
                    if not os.path.isdir(p):
                        LOGGER.debug("Failed to set env variables with %s " % root)
                        continue

                LOGGER.info("Setting MPI environment variable MPI_ROOT to %s" % mpi_root)
                LOGGER.info("Adding %s to PATH" % [mpi_path, mkl_path, cmp_path])

                os.environ["MPI_ROOT"] = mpi_root
                os.environ["PATH"] = (
                    ";".join([mpi_path, mkl_path, cmp_path]) + ";" + os.environ["PATH"]
                )
                os.environ["I_MPI_AUTH_METHOD"] = "delegate"
                os.environ["KMP_AFFINITY"] = "verbose"
                break

        elif self.dynatype == "msmpi" and self.platform == "windows":
            INTEL_CMP_REV = "2019.5.281"
            INTEL_MKL_REV = "2020.0.166"
            MS_MPI_REV = "10.1.12498.18"

            if os.getenv("MPI_ROOT"):
                LOGGER.warning(
                    "MPI_ROOT already defined. Not trying to automatically set MPI env variables."
                )
                return

            for root in ansys_env_roots:
                mpi_root = os.path.join(
                    os.getenv(root), "commonfiles", "MPI", "Microsoft", MS_MPI_REV, platform
                )
                mpi_path = os.path.join(mpi_root, "bin")
                mkl_path = os.path.join(os.getenv(root), "tp", "IntelMKL", INTEL_MKL_REV, platform)
                cmp_path = os.path.join(
                    os.getenv(root), "tp", "IntelCompiler", INTEL_CMP_REV, platform
                )
                for p in [mpi_root, mpi_path, mkl_path, cmp_path]:
                    if not os.path.isdir(p):
                        LOGGER.debug("Failed to set env variables with %s " % root)
                        continue

                LOGGER.info("Setting MPI environment variable MPI_ROOT to %s" % mpi_root)
                LOGGER.info("Adding %s to PATH" % [mpi_path, mkl_path, cmp_path])

                os.environ["MPI_ROOT"] = mpi_root
                os.environ["PATH"] = (
                    ";".join([mpi_path, mkl_path, cmp_path]) + ";" + os.environ["PATH"]
                )
                break

        elif self.dynatype == "platformmpi":
            LOGGER.error("Automatically setting env variables for platform mpi not implemented yet")
            return

        return

    def __repr__(self):
        """Represent self as string."""
        return yaml.dump(vars(self), allow_unicode=True, default_flow_style=False)


if __name__ == "__main__":
    LOGGER.debug("protected")
