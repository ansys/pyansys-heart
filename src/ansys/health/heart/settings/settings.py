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

"""Module that defines classes that hold settings relevant for PyAnsys-Heart."""

import copy
from dataclasses import asdict, dataclass, field
import json
import os
import pathlib
from pathlib import Path
import shutil
from typing import List, Literal

from pint import Quantity, UnitRegistry
import yaml

from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.exceptions import (
    LSDYNANotFoundError,
    MPIProgamNotFoundError,
    WSLNotFoundError,
)
from ansys.health.heart.settings.defaults import (
    electrophysiology as ep_defaults,
    fibers as fibers_defaults,
    mechanics as mech_defaults,
    purkinje as purkinje_defaults,
    zeropressure as zero_pressure_defaults,
)
from ansys.health.heart.settings.material.curve import (
    ActiveCurve,
    constant_ca2,
)
from ansys.health.heart.settings.material.material import (
    ACTIVE,
    ANISO,
    ISO,
    ActiveModel,
    Mat295,
)


class AttrDict(dict):
    """Dictionary subclass whose entries can be accessed by attributes as well as normally."""

    def __init__(self, *args, **kwargs):
        """Construct nested AttrDicts from nested dictionaries."""

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
        """Represent object in dictionary in YAML style."""
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
        """Convert units to a consistent unit system.

        Notes
        -----
        Currently the only supported unit system is ["MPa", "mm", "N", "ms", "g"]
        For instance:
        Quantity(10, "mm/s") --> Quantity(0.01, "mm/ms")
        """

        def _to_consitent_units(d):
            """Convert units to a consistent unit system."""
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
    """End time of simulation."""
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
    stiffness_damping: Quantity = Quantity(0, "s")
    """Stiffness damping constant."""


@dataclass(repr=False)
class EPAnalysis(Analysis):
    """Class for EP analysis settings."""

    solvertype: Literal["Monodomain", "Eikonal", "ReactionEikonal"] = "Monodomain"


@dataclass(repr=False)
class Material(Settings):
    """Class for storing material settings."""

    myocardium: AttrDict = None
    """Myocardium material."""
    passive: AttrDict = None
    """Passive material. For example, the vessel wall."""
    cap: AttrDict = None
    """Cap material."""


@dataclass(repr=False)
class EpMaterial(Settings):
    """Class for storing EP material settings."""

    myocardium: AttrDict = None
    """Myocardium material."""
    atrium: AttrDict = None
    """Atrial material."""
    cap: AttrDict = None
    """Cap material."""
    beam: AttrDict = None
    """Beam material."""
    # TODO: consider 'other', e.g passive conductor, soft tissue...?


@dataclass(repr=False)
class BoundaryConditions(Settings):
    """Stores settings/parameters for boundary conditions."""

    robin: AttrDict = None
    """Parameters for pericardium spring/damper b.c."""
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
    """Parameters for the left ventricle."""
    right_ventricle: AttrDict = None
    """Parameters for the right ventricle."""


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
    """Class for keeping track of zero-pressure analysis settings."""

    dt_nodout: Quantity = 0
    """Time interval of nodeout export."""

    max_iters: int = 3
    """Maximum iterations for stress-free-configuration algorithm."""
    method: int = 2
    """Method to use."""
    # TODO: this should be a Quantity type
    tolerance: float = 5.0
    """Tolerance to use for iterative algorithm."""


@dataclass(repr=False)
class ZeroPressure(Settings):
    """Class for keeping track of settings for stress-free-configuration computation."""

    analysis: AnalysisZeroPressure = field(default_factory=lambda: AnalysisZeroPressure())
    """Generic analysis settings."""


@dataclass
class Stimulation(Settings):
    """Stimulation settings."""

    node_ids: List[int] = None
    t_start: Quantity = Quantity(0.0, "ms")
    period: Quantity = Quantity(800, "ms")
    duration: Quantity = Quantity(2, "ms")
    amplitude: Quantity = Quantity(50, "uF/mm^3")

    def __setattr__(self, __name: str, __value) -> None:
        """Set attributes.

        Parameters
        ----------
        __name : str
            Attribute name.
        __value : _type_
            Attribute value.

        """
        if __name == "node_ids":
            if isinstance(__value, list):
                try:
                    __value = [int(x) for x in __value]
                except ValueError:
                    print("Failed to cast node_ids to list of integers.")

            return super().__setattr__(__name, __value)
        elif __name == "t_start" or __name == "period" or __name == "duration":
            return super().__setattr__(__name, Quantity(__value, "ms"))
        elif __name == "amplitude":
            return super().__setattr__(__name, Quantity(__value, "uF/mm^3"))


@dataclass(repr=False)
class Electrophysiology(Settings):
    """Class for keeping track of EP settings."""

    material: EpMaterial = field(default_factory=lambda: EpMaterial())
    """Material settings/configuration."""
    analysis: EPAnalysis = field(default_factory=lambda: EPAnalysis())
    """Generic analysis settings."""
    stimulation: AttrDict[str, Stimulation] = None
    """Stimulation settings."""


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
class AtrialFiber(Settings):
    """
    Class for keeping track of atrial fiber settings.

    Default parameters are from doi.org/10.1016/j.cma.2020.113468 for idealized geometry.
    """

    tau_mv: float = 0
    tau_lpv: float = 0
    tau_rpv: float = 0

    tau_tv: float = 0
    tau_raw: float = 0
    tau_ct_minus: float = 0
    tau_ct_plus: float = 0
    tau_icv: float = 0
    tau_scv: float = 0
    tau_ib: float = 0
    tau_ras: float = 0


@dataclass(repr=False)
class Purkinje(Settings):
    """Class for keeping track of Purkinje settings."""

    node_id_origin_left: int = None
    """Left Purkinje origin ID."""
    node_id_origin_right: int = None
    """Right Purkinje origin id."""
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

        Examples
        --------
        Instantiate settings and load defaults

        >>> from ansys.health.heart.settings.settings import SimulationSettings
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
            self.atrial_fibers: AtrialFiber = AtrialFiber()
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

        Examples
        --------
        Create examples settings with default values.

        >>> from ansys.health.heart.settings.settings import SimulationSettings
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

        Examples
        --------
        Create examples settings with default values.

        >>> from ansys.health.heart.settings.settings import SimulationSettings
        >>> settings = SimulationSettings()
        >>> settings.load_defaults()
        >>> settings.save("my_settings.yml")

        Load settings in second SimulationSettings object.

        >>> settings1 = SimulationSettings()
        >>> settings1.load("my_settings.yml")
        >>> print(
        ...     "True" if settings.mechanics.analysis == settings1.mechanics.analysis else "False"
        ... )
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
            analysis = Analysis()
            analysis.set_values(settings[attribute_name]["analysis"])
            material = Material()
            material.set_values(settings[attribute_name]["material"])
            boundary_conditions = BoundaryConditions()
            boundary_conditions.set_values(settings[attribute_name]["boundary_conditions"])
            system_model = SystemModel()
            system_model.set_values(settings[attribute_name]["system"])
            self.mechanics.analysis = analysis
            self.mechanics.material = material
            self.mechanics.boundary_conditions = boundary_conditions
            self.mechanics.system = system_model

            attribute_name = "stress_free"
            _deserialize_quantity(settings[attribute_name], ureg)
            analysis = AnalysisZeroPressure()
            analysis.set_values(settings[attribute_name]["analysis"])
            self.stress_free.analysis = analysis

        except KeyError:
            LOGGER.error("Failed to load mechanics settings.")

    def load_defaults(self):
        """Load default simulation settings.

        Examples
        --------
        Create examples settings with default values.

        Load module
        >>> from ansys.health.heart.settings.settings import SimulationSettings

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
        # TODO: move to Settings class
        for attr in self.__dict__:
            if isinstance(getattr(self, attr), Mechanics):
                analysis = Analysis()
                analysis.set_values(mech_defaults.analysis)
                material = Material()
                material.set_values(mech_defaults.material)
                boundary_conditions = BoundaryConditions()
                boundary_conditions.set_values(mech_defaults.boundary_conditions)
                system_model = SystemModel()
                system_model.set_values(mech_defaults.system_model)

                self.mechanics.analysis = analysis
                self.mechanics.material = material
                self.mechanics.boundary_conditions = boundary_conditions
                self.mechanics.system = system_model

            if isinstance(getattr(self, attr), ZeroPressure):
                analysis = AnalysisZeroPressure()
                analysis.set_values(zero_pressure_defaults.analysis)
                self.stress_free.analysis = analysis

            if isinstance(getattr(self, attr), Electrophysiology):
                analysis = EPAnalysis()
                analysis.set_values(ep_defaults.analysis)
                material = EpMaterial()
                material.set_values(ep_defaults.material)

                self.electrophysiology.analysis = analysis
                self.electrophysiology.material = material
                self.electrophysiology.stimulation: AttrDict[str, Stimulation] = AttrDict()
                for key in ep_defaults.stimulation.keys():
                    system_model = Stimulation()
                    system_model.set_values(ep_defaults.stimulation[key])
                    self.electrophysiology.stimulation[key] = system_model
                # TODO: add stim params, monodomain/bidomain/eikonal,cellmodel
            # TODO: add settings for purkinje  fibers and epmecha
            if isinstance(getattr(self, attr), Fibers):
                self.fibers.set_values(fibers_defaults.angles)
            if isinstance(getattr(self, attr), Purkinje):
                self.purkinje.set_values(purkinje_defaults.build)
            if isinstance(getattr(self, attr), AtrialFiber):
                self.atrial_fibers.set_values(fibers_defaults.la_bundle)
                self.atrial_fibers.set_values(fibers_defaults.ra_bundle)

    def to_consistent_unit_system(self):
        """Convert all settings to consistent unit-system ["MPa", "mm", "N", "ms", "g"].

        Examples
        --------
        Convert to the consistent unit system ["MPa", "mm", "N", "ms", "g"].

        Import necessary modules
        >>> from ansys.health.heart.settings.settings import SimulationSettings
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

    def get_mechanical_material(
        self, required_type: Literal["isotropic", "anisotropic"], ep_coupled=False
    ) -> Mat295:
        """Load mechanical materials from settings.

        Parameters
        ----------
        required_type : Literal[&#39;isotropic&#39;,&#39;anisotropic&#39;]
            Type of required maerial
        ep_coupled : bool, optional
            If MAT295 is coupled with EP simulation, by default False

        Returns
        -------
        MAT295
            material with parameters in settings
        """
        if required_type == "anisotropic":
            material = _read_myocardium_property(
                self.mechanics.material.myocardium, coupled=ep_coupled
            )
        elif required_type == "isotropic":
            material = _read_passive_property(self.mechanics.material.passive)

        return material

    def get_ventricle_fiber_rotation(self, method: Literal["LSDYNA", "D-RBM"]) -> dict:
        """Get rotation angles from settings.

        Parameters
        ----------
        method : Literal[&quot;LSDYNA&quot;, &quot;D
            Fiber rule based methods

        Returns
        -------
        dict
            rotation angles alpha and beta
        """
        if method == "LSDYNA":
            rotation = {
                "alpha": [
                    self.fibers.alpha_endo.m,
                    self.fibers.alpha_epi.m,
                ],
                "beta": [
                    self.fibers.beta_endo.m,
                    self.fibers.beta_epi.m,
                ],
                "beta_septum": [
                    self.fibers.beta_endo_septum.m,
                    self.fibers.beta_epi_septum.m,
                ],
            }
        elif method == "D-RBM":
            rotation = {
                "alpha_left": [
                    self.fibers.alpha_endo.m,
                    self.fibers.alpha_epi.m,
                ],
                "alpha_right": [
                    self.fibers.alpha_endo.m,
                    self.fibers.alpha_epi.m,
                ],
                "alpha_ot": None,
                "beta_left": [
                    self.fibers.beta_endo.m,
                    self.fibers.beta_epi.m,
                ],
                "beta_right": [
                    self.fibers.beta_endo.m,
                    self.fibers.beta_epi.m,
                ],
                "beta_ot": None,
            }
        return rotation


def _read_passive_property(passive: AttrDict) -> Mat295:
    """Read passive property from settings."""
    passive = Mat295(
        rho=passive["rho"].m,
        iso=ISO(
            itype=passive["itype"],
            beta=2,
            kappa=passive["kappa"].m,
            mu1=passive["mu1"].m,
            alpha1=passive["alpha1"],
        ),
    )
    return passive


def _read_myocardium_property(mat: AttrDict, coupled=False) -> Mat295:
    """Read myocardium property from settings."""
    rho = mat["isotropic"]["rho"].m

    iso = ISO(
        kappa=mat["isotropic"]["kappa"].m,
        k1=mat["isotropic"]["k1"].m,
        k2=mat["isotropic"]["k2"].m,
        beta=2,
    )

    fibers = [ANISO.HGOFiber(k1=mat["anisotropic"]["k1f"].m, k2=mat["anisotropic"]["k2f"].m)]

    if "k1s" in mat["anisotropic"]:
        sheet = ANISO.HGOFiber(k1=mat["anisotropic"]["k1s"].m, k2=mat["anisotropic"]["k2s"].m)
        fibers.append(sheet)

    if "k1fs" in mat["anisotropic"]:
        k1fs, k2fs = mat["anisotropic"]["k1fs"].m, mat["anisotropic"]["k2fs"].m
    else:
        k1fs, k2fs = None, None
    aniso = ANISO(fibers=fibers, k1fs=k1fs, k2fs=k2fs)

    max = mat["active"]["taumax"].m
    bt = mat["active"]["beat_time"].m
    ss = mat["active"]["ss"]
    sn = mat["active"]["sn"]

    if not coupled:
        ac_mdoel = ActiveModel.Model1(taumax=max)  # use default field in Model1 except taumax
        curve = ActiveCurve(constant_ca2(tb=bt), threshold=0.1, type="ca2")
        active = ACTIVE(
            ss=ss,
            sn=sn,
            model=ac_mdoel,
            ca2_curve=curve,
        )
    else:
        ac_mdoel = ActiveModel.Model3(
            ca2ion50=0.001,
            n=2,
            f=0.0,
            l=1.9,
            eta=1.45,
            sigmax=max,  # MPa
        )

        active = ACTIVE(
            ss=ss,
            sn=sn,
            acthr=0.0002,
            model=ac_mdoel,
            ca2_curve=None,
        )

    return Mat295(rho=rho, iso=iso, aniso=aniso, active=active)


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
        Quantity(30, "uF/mm^3").dimensionality,
    ],
    ["MPa", "N", "mS/mm", "uF/mm^2", "1/mS", "degree", "uF/mm^3"],
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


def _windows_to_wsl_path(windows_path: str):
    """Convert Windows to WSL path."""
    win_path = Path(windows_path)
    if isinstance(win_path, pathlib.PosixPath):
        return None

    if "\\\\wsl.localhost" in str(win_path):
        wsl_path = Path(*win_path.parts[1:])
        wsl_path = "/" + wsl_path.as_posix()
        return wsl_path

    elif win_path.drive != "":
        wsl_mount = ("/mnt/" + win_path.drive.replace(":", "")).lower()
        wsl_path = win_path.as_posix().replace(win_path.drive, wsl_mount)

    return wsl_path


class DynaSettings:
    """Class for collecting, managing, and validating LS-DYNA settings."""

    @staticmethod
    def _get_available_mpi_exe():
        """Find whether mpiexec or mpirun are available."""
        # preference for mpirun if it is added to PATH. mpiexec is the fallback option.
        if shutil.which("mpirun"):
            return shutil.which("mpirun")
        elif shutil.which("mpiexec"):
            LOGGER.debug("mpirun not found. Using mpiexec.")
            return shutil.which("mpiexec")
        else:
            raise MPIProgamNotFoundError("mpirun or mpiexec not found. Please configure MPI.")

    def __init__(
        self,
        lsdyna_path: pathlib.Path = "lsdyna.exe",
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
        """Type of LS-DYNA executable."""
        self.num_cpus: int = num_cpus
        """Number of CPU's requested."""
        self.platform: str = platform
        """Platform LS-DYNA is executed on."""

        self.dyna_options = dyna_options
        """Additional command line options for dyna."""

        if dynatype in ["intelmpi", "platformmpi", "msmpi"]:
            self.mpi_options = mpi_options
            """additional mpi options."""
        elif dynatype == "smp":
            self.mpi_options = ""

        self._modify_from_global_settings()
        LOGGER.info("LS-DYNA Configuration:")
        LOGGER.info(
            f"path: {self.lsdyna_path} | type: {self.dynatype} | platform: {self.platform} | cpus: {self.num_cpus}"  # noqa: E501
        )

        # Ensure path to LS-DYNA executable is absolute
        ls_dyna_abs_path = shutil.which(self.lsdyna_path)

        if self.platform == "wsl":
            ls_dyna_abs_path = str(Path(self.lsdyna_path).resolve())

        if ls_dyna_abs_path is None or not Path(ls_dyna_abs_path).is_file():
            raise LSDYNANotFoundError(
                f"LS-DYNA executable not found at {ls_dyna_abs_path}. Please check the path."
            )

        self.lsdyna_path: pathlib.Path = ls_dyna_abs_path

        if self.platform == "wsl" and os.name != "nt":
            raise WSLNotFoundError(f"""WSL is not supported on {os.name}.""")

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
        if self.platform == "wsl":
            mpi_exe = "mpirun"
        elif self.dynatype in ["msmpi", "intelmpi", "platformmpi"]:
            mpi_exe = self._get_available_mpi_exe()

        lsdyna_path = self.lsdyna_path

        if self.platform == "windows" or self.platform == "linux":
            if self.dynatype in ["intelmpi", "platformmpi"]:
                commands = [
                    mpi_exe,
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
            wsl_exe_path = shutil.which("wsl.exe")
            if wsl_exe_path is None:
                raise WSLNotFoundError("wsl.exe not found. Please install WSL.")

            # Convert paths to WSL compatible paths.
            path_to_input_wsl = _windows_to_wsl_path(path_to_input)
            lsdyna_path = _windows_to_wsl_path(self.lsdyna_path)

            if self.dynatype in ["intelmpi", "platformmpi", "msmpi"]:
                commands = [
                    mpi_exe,
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

            commands = [
                "powershell",
                "-Command",
                wsl_exe_path,
                "-e",
                "bash",
                "-lic",
                "./run_lsdyna.sh",
            ]

        # remove empty strings from commands
        commands = [c for c in commands if c != ""]

        # expand any environment variables if any
        commands = [os.path.expandvars(c) for c in commands]

        return commands

    def _modify_from_global_settings(self):
        """Set DynaSettings based on globally defined settings for PyAnsys-Heart."""
        keys = [key for key in os.environ.keys() if "PYANSYS_HEART" in key]
        LOGGER.debug(f"PYANSYS_HEART Environment variables: {keys}")
        self.lsdyna_path = os.getenv("PYANSYS_HEART_LSDYNA_PATH", self.lsdyna_path)
        self.platform = os.getenv("PYANSYS_HEART_LSDYNA_PLATFORM", self.platform)
        self.dynatype = os.getenv("PYANSYS_HEART_LSDYNA_TYPE", self.dynatype)
        self.num_cpus = int(os.getenv("PYANSYS_HEART_NUM_CPU", self.num_cpus))
        return

    def __repr__(self):
        """Represent self as string."""
        return yaml.dump(vars(self), allow_unicode=True, default_flow_style=False)
