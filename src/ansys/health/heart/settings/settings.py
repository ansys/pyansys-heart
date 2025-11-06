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

"""Module that defines classes that hold settings relevant for PyAnsys-Heart.

Examples
--------
Create and configure simulation settings:

>>> from ansys.health.heart.settings.settings import SimulationSettings
>>> settings = SimulationSettings()
>>> settings.load_defaults()
>>> settings.mechanics.analysis.end_time = Quantity(1000, "ms")
>>> settings.save("config.yml")

Load existing configuration:

>>> settings = SimulationSettings()
>>> settings.load("config.yml")
>>> print(settings.mechanics.analysis.end_time)
1000.0 millisecond
"""

import json
import os
import pathlib
from pathlib import Path
import shutil
from typing import Any, Literal

from pint import Quantity, UnitRegistry
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    ValidationError,
    field_serializer,
    field_validator,
)
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

ureg = UnitRegistry()


class BaseSettings(BaseModel):
    """Base class for all settings with Pydantic validation and serialization.

    Features
    --------
    - Automatic validation of types and values
    - Built-in JSON/YAML serialization with Pint Quantity support
    - Unit conversion to consistent unit system ["MPa", "mm", "N", "ms", "g"]
    - Nested model validation and type safety

    Examples
    --------
    >>> settings = BaseSettings()
    >>> settings.to_consistent_unit_system()
    >>> data = settings.model_dump_json()
    """

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        validate_assignment=True,
        extra="forbid",
        str_strip_whitespace=True,
    )

    def __repr__(self) -> str:
        """Represent object in YAML-style format using Pydantic v2 serialization.

        Returns
        -------
        str
            YAML-formatted string representation of the object.
        """
        data = self.model_dump(mode="json", exclude_none=True)
        data = {self.__class__.__name__: data}
        return yaml.dump(json.loads(json.dumps(data)), sort_keys=False)

    @field_serializer("*", when_used="json")
    def serialize_quantities_for_json(self, value: Any, _info) -> str | float | Any:
        """Serialize Quantity objects for JSON output.

        This serializer handles Quantity objects during JSON serialization,
        providing string representation.
        Handles nested Quantity objects in dictionaries and lists.

        Parameters
        ----------
        value : Any
            The field value to serialize.
        _info : SerializationInfo
            Pydantic serialization context (unused but required by signature).

        Returns
        -------
        str | float | Any
            String representation if Quantity, otherwise unchanged.
        """

        def _serialize_recursive(obj: Any) -> Any:
            """Recursively serialize Quantity objects in nested structures."""
            if isinstance(obj, Quantity):
                return str(obj)
            elif isinstance(obj, dict):
                return {key: _serialize_recursive(val) for key, val in obj.items()}
            elif isinstance(obj, (list, tuple)):
                return [_serialize_recursive(item) for item in obj]
            return obj

        return _serialize_recursive(value)

    @field_validator("*", mode="before")
    def parse_quantity(cls, v, info):  # noqa D102
        """Parse string values to Quantity objects for fields annotated as Quantity.

        This validator applies to all fields and attempts to parse string values
        as Quantity objects when the field is annotated with Quantity type.
        For nested models, it ensures proper Quantity parsing across all levels.

        Parameters
        ----------
        v : Any
            The value to validate and potentially convert to a Quantity.
        info : ValidationInfo
            Pydantic validation context containing field information.

        Returns
        -------
        Any
            Quantity object if conversion successful, otherwise the original value.
        """
        # If it's already a Quantity, return as-is
        if isinstance(v, Quantity):
            return v

        # Only attempt parsing for string values
        if not isinstance(v, str):
            return v

        # Get field annotation if available
        field_name = getattr(info, "field_name", None)
        if not field_name:
            return v

        # Check if this field should be a Quantity based on annotation
        if hasattr(cls, "__annotations__") and field_name in cls.__annotations__:
            field_annotation = cls.__annotations__[field_name]

            # Check if the field is annotated as Quantity
            if field_annotation == Quantity:
                try:
                    return ureg(v)
                except Exception as e:
                    LOGGER.warning(
                        f"Failed to parse quantity from string '{v}' for field '{field_name}': {e}"
                    )
                    return v

            # Handle generic aliases (for Python 3.9+ compatibility)
            if hasattr(field_annotation, "__origin__") and field_annotation.__origin__ == Quantity:
                try:
                    return ureg(v)
                except Exception as e:
                    LOGGER.warning(
                        f"Failed to parse quantity from string '{v}' for field '{field_name}': {e}"
                    )
                    return v

        # For non-Quantity fields, pass through unchanged
        return v

    def to_consistent_unit_system(self) -> None:
        """Convert units to consistent system ["MPa", "mm", "N", "ms", "g"].

        This method converts all Quantity objects to use the PyAnsys Heart
        standard unit system for cardiac simulations.

        Examples
        --------
        >>> from pint import Quantity
        >>> settings = BaseSettings()
        >>> # Assuming a Quantity field exists
        >>> settings.to_consistent_unit_system()
        """

        def _to_consistent_units(obj: Any) -> None:
            """Convert units recursively."""
            if isinstance(obj, BaseSettings):
                obj_dict = obj.__dict__
            elif isinstance(obj, dict):
                obj_dict = obj
            else:
                return

            for key, value in obj_dict.items():
                if isinstance(value, (dict, BaseSettings)):
                    _to_consistent_units(value)
                elif isinstance(value, Quantity) and not value.unitless:
                    if "[substance]" in list(value.dimensionality):
                        LOGGER.warning("Not converting [substance] / [length]^3")
                        continue
                    new_quantity = value.to(_get_consistent_units_str(value.dimensionality))
                    if isinstance(obj, BaseSettings):
                        setattr(obj, key, new_quantity)
                    else:
                        obj[key] = new_quantity

        _to_consistent_units(self)


class Analysis(BaseSettings):
    """Class for analysis settings.

    Defines core simulation analysis parameters including time stepping,
    output intervals, and damping parameters for cardiac simulations.

    Parameters
    ----------
    end_time : Quantity, default: Quantity(0, "s")
        End time of simulation in time units.
    dtmin : Quantity, default: Quantity(0, "s")
        Minimum time-step of simulation in time units.
    dtmax : Quantity, default: Quantity(0, "s")
        Maximum time-step of simulation in time units.
    dt_d3plot : Quantity, default: Quantity(0, "s")
        Time-step of d3plot export in time units.
    dt_icvout : Quantity, default: Quantity(0, "s")
        Time-step of icvout export in time units.
    global_damping : Quantity, default: Quantity(0, "1/s")
        Global damping constant in 1/time units.
    stiffness_damping : Quantity, default: Quantity(0, "s")
        Stiffness damping constant in time units.

    Examples
    --------
    >>> from pint import Quantity
    >>> analysis = Analysis(
    ...     end_time=Quantity(1000, "ms"), dtmin=Quantity(0.1, "ms"), dtmax=Quantity(10, "ms")
    ... )
    >>> analysis.to_consistent_unit_system()
    >>> print(analysis.end_time)
    1000.0 millisecond
    """

    end_time: Quantity = Field(default=Quantity(0, "s"), description="End time of simulation")
    dtmin: Quantity = Field(default=Quantity(0, "s"), description="Minimum time-step of simulation")
    dtmax: Quantity = Field(default=Quantity(0, "s"), description="Maximum time-step of simulation")
    dt_d3plot: Quantity = Field(default=Quantity(0, "s"), description="Time-step of d3plot export")
    dt_icvout: Quantity = Field(default=Quantity(0, "s"), description="Time-step of icvout export")
    global_damping: Quantity = Field(
        default=Quantity(0, "1/s"), description="Global damping constant"
    )
    stiffness_damping: Quantity = Field(
        default=Quantity(0, "s"), description="Stiffness damping constant"
    )


class EPAnalysis(Analysis):
    """Class for EP analysis settings.

    Extends Analysis with electrophysiology-specific solver configuration.
    Supports different EP solver types for cardiac electrical simulation.

    Parameters
    ----------
    solvertype : Literal["Monodomain", "Eikonal", "ReactionEikonal"], default: "Monodomain"
        Type of electrophysiology solver to use.

    Examples
    --------
    >>> ep_analysis = EPAnalysis(solvertype="Monodomain")
    >>> print(ep_analysis.solvertype)
    Monodomain
    """

    solvertype: Literal["Monodomain", "Eikonal", "ReactionEikonal"] = Field(
        default="Monodomain", description="Type of electrophysiology solver"
    )


class BoundaryConditions(BaseSettings):
    """Stores settings/parameters for boundary conditions.

    Manages boundary condition parameters for cardiac simulation including
    pericardium constraints, valve mechanics, and pressure loading.

    Parameters
    ----------
    robin : dict[str, Any] | None, default: None
        Parameters for pericardium spring/damper boundary conditions.
    valve : dict[str, Any] | None, default: None
        Parameters for valve spring boundary conditions.
    end_diastolic_cavity_pressure : dict[str, Any] | None, default: None
        End-diastolic pressure configuration.

    Examples
    --------
    >>> bc = BoundaryConditions(
    ...     robin={"stiffness": 1.0, "damping": 0.1}, valve={"spring_constant": 100.0}
    ... )
    >>> print(bc.robin["stiffness"])
    1.0
    """

    robin: dict[str, Any] | None = Field(
        default=None, description="Parameters for pericardium spring/damper b.c."
    )
    valve: dict[str, Any] | None = Field(
        default=None, description="Parameters for valve spring b.c."
    )
    end_diastolic_cavity_pressure: dict[str, Any] | None = Field(
        default=None, description="End-diastolic pressure"
    )


class SystemModel(BaseSettings):
    """Stores settings/parameters for the system model.

    Manages system-level model configuration including circulatory system
    models and ventricular-specific parameters.

    Parameters
    ----------
    name : str, default: "ConstantPreloadWindkesselAfterload"
        Name of the system model implementation.
    left_ventricle : dict[str, Any] | None, default: None
        Parameters specific to left ventricle modeling.
    right_ventricle : dict[str, Any] | None, default: None
        Parameters specific to right ventricle modeling.

    Examples
    --------
    >>> system = SystemModel(
    ...     name="ConstantPreloadWindkesselAfterload",
    ...     left_ventricle={"volume": 150.0},
    ...     right_ventricle={"volume": 120.0},
    ... )
    >>> print(system.name)
    ConstantPreloadWindkesselAfterload
    """

    name: str = Field(
        default="ConstantPreloadWindkesselAfterload", description="Name of the system model"
    )
    left_ventricle: dict[str, Any] | None = Field(
        default=None, description="Parameters for the left ventricle"
    )
    right_ventricle: dict[str, Any] | None = Field(
        default=None, description="Parameters for the right ventricle"
    )


class Mechanics(BaseSettings):
    """Class for keeping track of mechanical simulation settings.

    Complete mechanical simulation configuration including analysis parameters,
    boundary conditions, and system model settings.

    Parameters
    ----------
    analysis : Analysis, default: Analysis()
        Generic analysis settings for time stepping and output.
    boundary_conditions : BoundaryConditions, default: BoundaryConditions()
        Boundary condition specifications and parameters.
    system : SystemModel, default: SystemModel()
        System model settings and configurations.

    Examples
    --------
    >>> mechanics = Mechanics()
    >>> mechanics.analysis.end_time = Quantity(1000, "ms")
    >>> mechanics.boundary_conditions.robin = {"stiffness": 1.0}
    >>> print(mechanics.analysis.end_time)
    1000.0 millisecond
    """

    analysis: Analysis = Field(default_factory=Analysis, description="Generic analysis settings")
    boundary_conditions: BoundaryConditions = Field(
        default_factory=BoundaryConditions, description="Boundary condition specifications"
    )
    system: SystemModel = Field(default_factory=SystemModel, description="System model settings")


class AnalysisZeroPressure(Analysis):
    """Class for keeping track of zero-pressure analysis settings.

    Extends Analysis with specific settings for stress-free configuration
    computation, including iterative solver parameters.

    Parameters
    ----------
    dt_nodout : Quantity, default: Quantity(0, "s")
        Time interval of nodeout export.
    max_iters : int, default: 3
        Maximum iterations for stress-free-configuration algorithm.
    method : int, default: 2
        Method identifier to use for computation.
    tolerance : float, default: 5.0
        Tolerance for iterative algorithm convergence.

    Examples
    --------
    >>> zero_p = AnalysisZeroPressure(max_iters=5, tolerance=1.0, method=2)
    >>> print(zero_p.max_iters)
    5
    """

    dt_nodout: Quantity = Field(
        default=Quantity(0, "s"), description="Time interval of nodeout export"
    )
    max_iters: int = Field(
        default=3, description="Maximum iterations for stress-free-configuration algorithm"
    )
    method: int = Field(default=2, description="Method to use")
    tolerance: float = Field(default=5.0, description="Tolerance to use for iterative algorithm")


class ZeroPressure(BaseSettings):
    """Class for keeping track of settings for stress-free-configuration computation.

    Configuration for computing the stress-free (unloaded) configuration
    of cardiac geometry, essential for accurate mechanical simulations.

    Parameters
    ----------
    analysis : AnalysisZeroPressure, default: AnalysisZeroPressure()
        Analysis settings specific to zero-pressure computation.

    Examples
    --------
    >>> zero_pressure = ZeroPressure()
    >>> zero_pressure.analysis.max_iters = 5
    >>> zero_pressure.analysis.tolerance = 1.0
    >>> print(zero_pressure.analysis.method)
    2
    """

    analysis: AnalysisZeroPressure = Field(
        default_factory=AnalysisZeroPressure, description="Generic analysis settings"
    )


class Stimulation(BaseSettings):
    """Stimulation settings for electrophysiology simulations.

    Defines electrical stimulation parameters including timing, location,
    and amplitude for cardiac electrophysiology simulations.

    Parameters
    ----------
    node_ids : list[int] | None, default: None
        List of node IDs where stimulation is applied.
    t_start : Quantity, default: Quantity(0.0, "ms")
        Start time of stimulation.
    period : Quantity, default: Quantity(800, "ms")
        Period between stimulation cycles.
    duration : Quantity, default: Quantity(2, "ms")
        Duration of each stimulation pulse.
    amplitude : Quantity, default: Quantity(50, "uF/mm^3")
        Stimulation amplitude.

    Examples
    --------
    >>> from pint import Quantity
    >>> stim = Stimulation(
    ...     node_ids=[1, 2, 3],
    ...     t_start=Quantity(10, "ms"),
    ...     period=Quantity(800, "ms"),
    ...     duration=Quantity(2, "ms"),
    ...     amplitude=Quantity(50, "uF/mm^3"),
    ... )
    >>> print(stim.node_ids)
    [1, 2, 3]
    """

    node_ids: list[int] | None = Field(default=None, description="Node IDs for stimulation")
    t_start: Quantity = Field(default=Quantity(0.0, "ms"), description="Start time of stimulation")
    period: Quantity = Field(default=Quantity(800, "ms"), description="Period between cycles")
    duration: Quantity = Field(default=Quantity(2, "ms"), description="Duration of pulse")
    amplitude: Quantity = Field(
        default=Quantity(50, "uF/mm^3"), description="Stimulation amplitude"
    )

    @field_validator("node_ids")
    @classmethod
    def validate_node_ids(cls, v: Any) -> list[int] | None:
        """Validate and convert node_ids to list of integers.

        Parameters
        ----------
        v : Any
            Input value to validate.

        Returns
        -------
        list[int] | None
            Validated list of integer node IDs or None.

        Raises
        ------
        ValueError
            If node_ids cannot be converted to list of integers.
        """
        if v is None:
            return None
        if isinstance(v, list):
            try:
                return [int(x) for x in v]
            except (ValueError, TypeError) as e:
                raise ValueError("Failed to cast node_ids to list of integers") from e
        raise ValueError("node_ids must be a list of integers or None")


class Electrophysiology(BaseSettings):
    """Class for keeping track of EP settings.

    Complete electrophysiology simulation configuration including analysis settings,
    stimulation protocols, layer definitions, and conductivity parameters.

    Parameters
    ----------
    analysis : EPAnalysis, default EPAnalysis()
        Generic analysis settings for EP simulation.
    stimulation : dict[str, Stimulation] | None, default None
        Dictionary of stimulation settings by name.
    layers : dict[str, Quantity]
        Layer definitions for material assignment of myocardium.
        Default: {"percent_endo": Quantity(0.17, "dimensionless"),
        "percent_mid": Quantity(0.41, "dimensionless")}
    lambda_ratio : Quantity, default Quantity(0.2, "dimensionless")
        Intra to extracellular conductivity ratio for EP solve.

    Examples
    --------
    >>> ep = Electrophysiology()
    >>> ep.analysis.solvertype = "Monodomain"
    >>> ep.stimulation = {"apex": Stimulation(node_ids=[1, 2, 3])}
    >>> print(ep.analysis.solvertype)
    Monodomain
    """

    analysis: EPAnalysis = Field(
        default_factory=EPAnalysis, description="Generic analysis settings"
    )
    stimulation: dict[str, Stimulation] | None = Field(
        default=None, description="Stimulation settings"
    )

    layers: dict[str, Quantity] = Field(
        default_factory=lambda: {
            "percent_endo": Quantity(0.17, "dimensionless"),  # thickness of endocardial layer
            "percent_mid": Quantity(0.41, "dimensionless"),  # thickness of midmyocardial layer
        },
        description="Layers for material assignment of the myocardium",
    )

    lambda_ratio: Quantity = Field(
        default=Quantity(0.2, "dimensionless"),
        description="Intra to extracellular conductivity ratio",
    )


class BaseFiberSettings(BaseSettings):
    """Base class for keeping track of fiber orientation settings.

    Defines fundamental fiber orientation parameters for cardiac muscle fiber
    modeling. These settings control the helical and transverse angles that
    define the spatial orientation of cardiac muscle fibers across the
    myocardial wall from endocardium to epicardium.

    Parameters
    ----------
    alpha_endo : Quantity, default: Quantity(-60, "degree")
        Helical angle in endocardium (inner heart wall surface) in degrees.
        Positive values indicate right-handed helix orientation.
    alpha_epi : Quantity, default: Quantity(60, "degree")
        Helical angle in epicardium (outer heart wall surface) in degrees.
        Typically opposite sign to alpha_endo for transmural rotation.
    beta_endo : Quantity, default: Quantity(-65, "degree")
        Angle to the outward transmural axis in endocardium in degrees.
        Controls fiber inclination relative to heart wall thickness direction.
    beta_epi : Quantity, default: Quantity(25, "degree")
        Angle to the outward transmural axis in epicardium in degrees.
        Defines fiber inclination at the outer wall surface.

    Notes
    -----
    The fiber orientation by rotating a the local coordinate system:
    - Longitudinal direction (e_l): apex to base
    - Transmural direction (e_t): endocardium to epicardium
    - Circumferential direction (e_c): orthogonal to both (right-hand rule)

    Alpha angles define rotation of the circumferential direction around the
    transmural axis (helical angle).
    Beta angles define rotation of the circumferential direction around the longitudinal
    axis (inclination angle).

    Default values are based on Bayer et al. 2012 http://dx.doi.org/10.1007/s10439-012-0593-5. Note
    that for the Bayer method the transmural direction points inward from epicardium, hence the
    negative sign of alpha_endo and beta_endo.

    Examples
    --------
    Create basic fiber settings with default Bayer et al. values:

    >>> fiber_settings = BaseFiberSettings()
    >>> print(fiber_settings.alpha_endo)
    -60 degree
    >>> print(fiber_settings.alpha_epi)
    60 degree

    Create custom fiber orientation:

    >>> from pint import Quantity
    >>> custom_fibers = BaseFiberSettings(
    ...     alpha_endo=Quantity(-45, "degree"),
    ...     alpha_epi=Quantity(45, "degree"),
    ...     beta_endo=Quantity(-30, "degree"),
    ...     beta_epi=Quantity(30, "degree"),
    ... )
    """

    alpha_endo: Quantity = Field(
        default=Quantity(-60, "degree"), description="Helical angle in endocardium"
    )
    alpha_epi: Quantity = Field(
        default=Quantity(60, "degree"), description="Helical angle in epicardium"
    )
    beta_endo: Quantity = Field(
        default=Quantity(-65, "degree"),
        description="Angle to the outward transmural axis in endocardium",
    )
    beta_epi: Quantity = Field(
        default=Quantity(25, "degree"),
        description="Angle to the outward transmural axis in epicardium",
    )


class FibersBRBM(BaseFiberSettings):
    """Class for keeping track of fiber settings for the Bayer et al rule-based method.

    Extends BaseFiberSettings with additional septum-specific fiber orientations
    required for the LS-DYNA/B-RBM fiber generation method.

    The B-RBM method uses distinct angle specifications for the septum region
    to account for the unique fiber architecture in the interventricular septum.

    Parameters
    ----------
    alpha_endo : Quantity, default: Quantity(-60, "degree")
        Helical angle in endocardium (inherited from BaseFiberSettings).
    alpha_epi : Quantity, default: Quantity(60, "degree")
        Helical angle in epicardium (inherited from BaseFiberSettings).
    beta_endo : Quantity, default: Quantity(-65, "degree")
        Angle to the outward transmural axis in endocardium (inherited).
    beta_epi : Quantity, default: Quantity(25, "degree")
        Angle to the outward transmural axis in epicardium (inherited).
    beta_endo_septum : Quantity, default: Quantity(-65, "degree")
        Angle to the outward transmural axis on the left septum endocardium.
        Specific to the septal region fiber orientation.
    beta_epi_septum : Quantity, default Quantity(25, "degree")
        Angle to the outward transmural axis in the septum epicardium.
        Controls septal fiber inclination at the epicardial surface.

    Notes
    -----
    Based on Bayer et al. https://doi.org/10.1007/s10439-012-0593-5.

    The B-RBM method distinguishes between:
    - Free wall fiber orientations (using base class angles)
    - Septal fiber orientations (using septum-specific beta angles)

    This allows for realistic modeling of the complex fiber architecture
    in the interventricular septum. Note that for the Bayer method the transmural
    direction points inward from epicardium to endocardium, and hence positive rotation
    has a different meaning that for the D-RBM method.

    Examples
    --------
    Create B-RBM fiber settings with default values:

    >>> brbm_fibers = FibersBRBM()
    >>> print(brbm_fibers.beta_endo_septum)
    -65 degree
    >>> print(brbm_fibers.beta_epi_septum)
    25 degree

    Create custom B-RBM settings:

    >>> from pint import Quantity
    >>> custom_brbm = FibersBRBM(
    ...     alpha_endo=Quantity(-50, "degree"),
    ...     alpha_epi=Quantity(50, "degree"),
    ...     beta_endo_septum=Quantity(-70, "degree"),
    ...     beta_epi_septum=Quantity(30, "degree"),
    ... )
    """

    beta_endo_septum: Quantity = Field(
        default=Quantity(-65, "degree"),
        description="Angle to the outward transmural axis on the left septum",
    )
    beta_epi_septum: Quantity = Field(
        default=Quantity(25, "degree"),
        description="Angle to the outward transmural axis in the septum",
    )

    def _get_rotation_dict(self) -> dict[str, list[float]]:
        """Get B-RBM rotation angles formatted for legacy compute_fibers method.

        Converts the Pydantic model data to the dictionary format expected by
        the compute_fibers(method="LSDYNA") method for B-RBM fiber generation.

        Returns
        -------
        dict[str, list[float]]
            Dictionary with keys "alpha", "beta", "beta_septum" containing
            [endo, epi] angle pairs in degrees for B-RBM method.
        """
        return {
            "alpha": [
                self.alpha_endo.to("degree").magnitude,
                self.alpha_epi.to("degree").magnitude,
            ],
            "beta": [
                self.beta_endo.to("degree").magnitude,
                self.beta_epi.to("degree").magnitude,
            ],
            "beta_septum": [
                self.beta_endo_septum.to("degree").magnitude,
                self.beta_epi_septum.to("degree").magnitude,
            ],
        }


class FibersDRBM(BaseSettings):
    """Class for storing settings for the Doste et al rule-based method.

    Implements the D-RBM (Doste Rule-Based Method) for fiber orientation
    generation in biventricular cardiac models. This method uses separate
    fiber orientation settings for left and right ventricles, with optional
    outflow tract specifications. Moreover, it includes a septal fraction
    parameter to define the portion of the septum assigned to the left ventricle.

    The D-RBM method provides ventricle-specific fiber orientations to
    better represent the distinct fiber architectures in each ventricle,
    particularly important for accurate mechanical and electrical modeling.

    Parameters
    ----------
    left_ventricle : BaseFiberSettings, default: BaseFiberSettings(alpha_endo=Quantity(60, "degree")
        ,alpha_epi=Quantity(-60, "degree"), beta_endo=Quantity(-20, "degree"),
        beta_epi=Quantity(20, "degree"))
        Fiber orientation settings specific to the left ventricle.
        Contains alpha and beta angles for left ventricular wall.
    right_ventricle : BaseFiberSettings: default: BaseFiberSettings(alpha_endo=
        Quantity(-90, "degree"),alpha_epi=Quantity(25, "degree"), beta_endo=Quantity(0, "degree"),
        beta_epi=Quantity(20, "degree"))
        Fiber orientation settings specific to the right ventricle.
        Contains alpha and beta angles for right ventricular wall.
    alpha_outflow_tract : Quantity | None, default: None
        Helical angle for the outflow tract region in degrees.
        Set to None if outflow tract fiber orientation not specified.
    beta_outflow_tract : Quantity | None, default: None
        Inclination angle for the outflow tract region in degrees.
        Set to None if outflow tract fiber orientation not specified.
    septal_fraction : float, default: 2.0/3.0
        The fraction of the septum that belongs to the left ventricle.
        Typically 2/3 (0.667) based on anatomical measurements.

    Notes
    -----
    Based on Doste et al. https://doi.org/10.1002/cnm.3185.

    Coordinate system definition in PyAnsys-Heart:
    - Longitudinal direction: apex to base
    - Transmural direction: endo to epicardium
    - Circumferential direction: e_c = e_l × e_t (right-hand rule)

    Alpha defines rotation around transmural axis (helical angle).
    Beta defines rotation around longitudinal axis (inclination angle).

    Examples
    --------
    Create D-RBM fiber settings with default values:

    >>> drbm_fibers = FibersDRBM()
    >>> print(drbm_fibers.left_ventricle.alpha_endo)
    60 degree
    >>> print(drbm_fibers.right_ventricle.alpha_endo)
    -90 degree
    >>> print(drbm_fibers.septal_fraction)
    0.6666666666666666

    Create custom D-RBM settings:

    >>> from pint import Quantity
    >>> custom_drbm = FibersDRBM(
    ...     left_ventricle=BaseFiberSettings(
    ...         alpha_endo=Quantity(45, "degree"), alpha_epi=Quantity(-45, "degree")
    ...     ),
    ...     septal_fraction=0.7,
    ... )
    """

    left_ventricle: BaseFiberSettings = Field(
        default_factory=lambda: BaseFiberSettings(
            alpha_endo=Quantity(60, "degree"),
            alpha_epi=Quantity(-60, "degree"),
            beta_endo=Quantity(-20, "degree"),
            beta_epi=Quantity(20, "degree"),
        ),
        description="Fiber orientation settings for the left ventricle",
    )
    right_ventricle: BaseFiberSettings = Field(
        default_factory=lambda: BaseFiberSettings(
            alpha_endo=Quantity(-90, "degree"),
            alpha_epi=Quantity(25, "degree"),
            beta_endo=Quantity(0, "degree"),
            beta_epi=Quantity(20, "degree"),
        ),
        description="Fiber orientation settings for the right ventricle",
    )
    alpha_outflow_tract: Quantity | None = Field(
        default=None, description="Helical angle for outflow tract region (None if not specified)"
    )
    beta_outflow_tract: Quantity | None = Field(
        default=None,
        description="Inclination angle for outflow tract region (None if not specified)",
    )
    septal_fraction: float = Field(
        default=2.0 / 3.0,
        description="Fraction of septum belonging to left ventricle (typically 2/3)",
    )

    def _get_rotation_dict(self) -> dict[str, list[float] | None]:
        """Get D-RBM rotation angles formatted for legacy compute_fibers method.

        Converts the Pydantic model data to the dictionary format expected by
        the compute_fibers(method="D-RBM") method.

        Returns
        -------
        dict[str, list[float] | None]
            Dictionary with keys "alpha_left", "alpha_right", "alpha_ot",
            "beta_left", "beta_right", "beta_ot" containing [endo, epi] angle
            pairs in degrees, or None for outflow tract if not configured.

        Examples
        --------
        >>> settings = SimulationSettings(fiber_method="D-RBM")
        >>> drbm_fibers = settings.get_fibers_drbm()
        >>> rotation_angles = drbm_fibers.get_rotation_dict()
        >>> print(rotation_angles["alpha_left"])
        [60.0, -60.0]
        """
        return {
            "alpha_left": [
                self.left_ventricle.alpha_endo.to("degree").magnitude,
                self.left_ventricle.alpha_epi.to("degree").magnitude,
            ],
            "alpha_right": [
                self.right_ventricle.alpha_endo.to("degree").magnitude,
                self.right_ventricle.alpha_epi.to("degree").magnitude,
            ],
            "alpha_ot": (
                [self.alpha_outflow_tract.to("degree").magnitude] * 2
                if self.alpha_outflow_tract is not None
                else None
            ),
            "beta_left": [
                self.left_ventricle.beta_endo.to("degree").magnitude,
                self.left_ventricle.beta_epi.to("degree").magnitude,
            ],
            "beta_right": [
                self.right_ventricle.beta_endo.to("degree").magnitude,
                self.right_ventricle.beta_epi.to("degree").magnitude,
            ],
            "beta_ot": (
                [self.beta_outflow_tract.to("degree").magnitude] * 2
                if self.beta_outflow_tract is not None
                else None
            ),
        }


class AtrialFiber(BaseSettings):
    """Class for keeping track of atrial fiber settings.

    Default parameters are from doi.org/10.1016/j.cma.2020.113468 for idealized geometry.
    Defines atrial fiber bundle parameters and orientations.

    Parameters
    ----------
    tau_mv : float, default: 0.0
        Mitral valve parameter.
    tau_lpv : float, default: 0.0
        Left pulmonary vein parameter.
    tau_rpv : float, default: 0.0
        Right pulmonary vein parameter.
    tau_tv : float, default: 0.0
        Tricuspid valve parameter.
    tau_raw : float, default: 0.0
        Right atrial wall parameter.
    tau_ct_minus : float, default: 0.0
        Crista terminalis minus parameter.
    tau_ct_plus : float, default: 0.0
        Crista terminalis plus parameter.
    tau_icv : float, default: 0.0
        Inferior vena cava parameter.
    tau_scv : float, default: 0.0
        Superior vena cava parameter.
    tau_ib : float, default: 0.0
        Isthmus bundle parameter.
    tau_ras : float, default: 0.0
        Right atrial septum parameter.

    Examples
    --------
    >>> atrial = AtrialFiber(tau_mv=0.5, tau_tv=0.3)
    >>> print(atrial.tau_mv)
    0.5
    """

    tau_mv: float = Field(default=0.0, description="Mitral valve parameter")
    tau_lpv: float = Field(default=0.0, description="Left pulmonary vein parameter")
    tau_rpv: float = Field(default=0.0, description="Right pulmonary vein parameter")
    tau_tv: float = Field(default=0.0, description="Tricuspid valve parameter")
    tau_raw: float = Field(default=0.0, description="Right atrial wall parameter")
    tau_ct_minus: float = Field(default=0.0, description="Crista terminalis minus parameter")
    tau_ct_plus: float = Field(default=0.0, description="Crista terminalis plus parameter")
    tau_icv: float = Field(default=0.0, description="Inferior vena cava parameter")
    tau_scv: float = Field(default=0.0, description="Superior vena cava parameter")
    tau_ib: float = Field(default=0.0, description="Isthmus bundle parameter")
    tau_ras: float = Field(default=0.0, description="Right atrial septum parameter")


class Purkinje(BaseSettings):
    """Class for keeping track of Purkinje settings.

    Defines parameters for Purkinje network generation and electrical
    properties including geometry, branching, and junction characteristics.

    Parameters
    ----------
    node_id_origin_left : int | None, default: None
        Left Purkinje origin node ID.
    node_id_origin_right : int | None, default: None
        Right Purkinje origin node ID.
    edgelen : Quantity, default: Quantity(0, "mm")
        Edge length for Purkinje segments.
    ngen : Quantity, default: Quantity(0, "dimensionless")
        Number of generations in the network.
    nbrinit : Quantity, default: Quantity(0, "dimensionless")
        Number of initial branches from origin.
    nsplit : Quantity, default: Quantity(0, "dimensionless")
        Number of splits at each leaf.
    pmjtype : Quantity, default: Quantity(0, "dimensionless")
        Purkinje muscle junction type identifier.
    pmjradius : Quantity, default: Quantity(0, "mm")
        Purkinje muscle junction radius.
    pmjrestype : Quantity, default: Quantity(1, "dimensionless")
        Purkinje muscle junction resistance type.
    pmjres : Quantity, default: Quantity(0.001, "1/mS")
        Purkinje muscle junction resistance value.

    Examples
    --------
    >>> from pint import Quantity
    >>> purkinje = Purkinje(
    ...     node_id_origin_left=1,
    ...     node_id_origin_right=2,
    ...     edgelen=Quantity(1.0, "mm"),
    ...     ngen=Quantity(5, "dimensionless"),
    ... )
    >>> print(purkinje.node_id_origin_left)
    1
    """

    node_id_origin_left: int | None = Field(default=None, description="Left Purkinje origin ID")
    node_id_origin_right: int | None = Field(default=None, description="Right Purkinje origin id")
    edgelen: Quantity = Field(default=Quantity(0, "mm"), description="Edge length")
    ngen: Quantity = Field(
        default=Quantity(0, "dimensionless"), description="Number of generations"
    )
    nbrinit: Quantity = Field(
        default=Quantity(0, "dimensionless"), description="Number of beams from origin point"
    )
    nsplit: Quantity = Field(
        default=Quantity(0, "dimensionless"), description="Number of splits at each leaf"
    )
    pmjtype: Quantity = Field(
        default=Quantity(0, "dimensionless"), description="Purkinje muscle junction type"
    )
    pmjradius: Quantity = Field(
        default=Quantity(0, "mm"), description="Purkinje muscle junction radius"
    )
    pmjrestype: Quantity = Field(
        default=Quantity(1, "dimensionless"), description="Purkinje muscle junction resistance type"
    )
    pmjres: Quantity = Field(
        default=Quantity(0.001, "1/mS"), description="Purkinje muscle junction resistance"
    )


class SimulationSettings:
    """Class for keeping track of settings.

    Attributes are conditionally created based on initialization parameters.
    All parameters default to True, so all attributes exist by default.
    """

    # Type annotations for conditionally created attributes
    # Note: These attributes will only exist if the corresponding boolean parameter is True
    # All parameters default to True, so these attributes exist in the default case
    mechanics: Mechanics  # Exists when mechanics=True (default)
    electrophysiology: Electrophysiology  # Exists when electrophysiology=True (default)
    fibers: FibersBRBM | FibersDRBM  # Exists when fiber=True (default)
    atrial_fibers: AtrialFiber  # Exists when fiber=True (default)
    purkinje: Purkinje  # Exists when purkinje=True (default)
    stress_free: ZeroPressure  # Exists when stress_free=True (default)

    def __init__(
        self,
        mechanics: bool = True,
        electrophysiology: bool = True,
        fiber: bool = True,
        fiber_method: Literal["LSDYNA", "D-RBM"] = "LSDYNA",
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
            if fiber_method == "LSDYNA":
                self.fibers: FibersBRBM = FibersBRBM()
                """Fiber settings for the LSDYNA rule-based method."""
            elif fiber_method == "D-RBM":
                self.fibers: FibersDRBM = FibersDRBM()
                """Fiber settings for the D-RBM method."""
            else:
                raise ValueError(
                    "Invalid method to compute the fiber orientation. "
                    "Valid methods include: [LSDYNA, D-RBM]"
                )

            # Store the fiber method for later validation and loading
            self._fiber_method = fiber_method

            self.atrial_fibers: AtrialFiber = AtrialFiber()
            """Settings for atrial fiber generation."""
        else:
            self._fiber_method = None

        if purkinje:
            self.purkinje: Purkinje = Purkinje()
            """Settings for Purkinje generation."""

        if stress_free:
            self.stress_free: ZeroPressure = ZeroPressure()
            """Settings for stress free configuration simulation."""

        return

    def _get_fiber_config_lsdyna(self) -> FibersBRBM:
        """Get LSDYNA fiber settings with type safety.

        Returns
        -------
        FibersBRBM
            LSDYNA fiber settings instance.

        Raises
        ------
        ValueError
            If fiber method is not LSDYNA or fibers not configured.
        """
        if not hasattr(self, "fibers"):
            raise ValueError("Fiber settings not configured")
        if self._fiber_method != "LSDYNA":
            raise ValueError(f"Fiber method is {self._fiber_method}, not LSDYNA")
        return self.fibers  # type: ignore (we know it's FibersBRBM from the check)

    def _get_fiber_config_drbm(self) -> FibersDRBM:
        """Get D-RBM fiber settings with type safety.

        Returns
        -------
        FibersDRBM
            D-RBM fiber settings instance.

        Raises
        ------
        ValueError
            If fiber method is not D-RBM or fibers not configured.
        """
        if not hasattr(self, "fibers"):
            raise ValueError("Fiber settings not configured")
        if self._fiber_method != "D-RBM":
            raise ValueError(f"Fiber method is {self._fiber_method}, not D-RBM")
        return self.fibers  # type: ignore (we know it's FibersDRBM from the check)

    def __repr__(self):
        """Represent object as list of relevant attribute names.

        Returns
        -------
        str
            String representation showing the class name and active
            settings attribute names.

        Examples
        --------
        >>> settings = SimulationSettings()
        >>> print(repr(settings))
        SimulationSettings
          mechanics
          electrophysiology
          fibers
          atrial_fibers
          purkinje
          stress_free
        """
        repr_str = "\n  ".join(
            [attr for attr in self.__dict__ if isinstance(getattr(self, attr), BaseSettings)]
        )
        repr_str = self.__class__.__name__ + "\n  " + repr_str
        return repr_str

    def save(self, filename: pathlib.Path):
        """Save simulation settings to disk.

        Parameters
        ----------
        filename : pathlib.Path
            Path to target .json or .yml file

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

        # Serialize each of the settings using Pydantic v2's enhanced model_dump
        serialized_settings = {}
        for attribute_name in self.__dict__.keys():
            if not isinstance(getattr(self, attribute_name), BaseSettings):
                continue
            else:
                setting: BaseSettings = getattr(self, attribute_name)
                # Use the simplified model dump method (no unit removal)
                serialized_settings[attribute_name] = setting.model_dump(
                    mode="json", exclude_none=False
                )

        serialized_settings = {"Simulation Settings": serialized_settings}

        with open(filename, "w") as f:
            if filename.suffix == ".yml":
                # Serialize settings using modern Pydantic serialization
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

        # Load file data with proper error handling
        try:
            with open(filename, "r", encoding="utf-8") as f:
                if filename.suffix == ".json":
                    data = json.load(f)
                elif filename.suffix == ".yml":
                    data = yaml.load(f, Loader=yaml.SafeLoader)
                else:
                    raise ValueError(f"Unsupported file format: {filename.suffix}")
        except FileNotFoundError as e:
            LOGGER.error(f"Settings file not found: {filename}")
            raise FileNotFoundError(f"Settings file not found: {filename}") from e
        except (json.JSONDecodeError, yaml.YAMLError) as e:
            LOGGER.error(f"Failed to parse settings file {filename}: {e}")
            raise ValueError(f"Invalid file format in {filename}: {e}") from e

        settings_data = data.get("Simulation Settings", {})
        if not settings_data:
            LOGGER.warning("No 'Simulation Settings' found in file")
            return

        try:
            # Use streamlined approach - Pydantic handles all validation automatically
            self._load_settings_section("mechanics", settings_data, Mechanics)
            self._load_settings_section("stress_free", settings_data, ZeroPressure)
            self._load_settings_section("electrophysiology", settings_data, Electrophysiology)

            # Handle fiber settings conditionally based on detected method
            self._load_fiber_settings(settings_data)

            self._load_settings_section("atrial_fibers", settings_data, AtrialFiber)
            self._load_settings_section("purkinje", settings_data, Purkinje)

        except ValidationError as e:
            LOGGER.error(f"Validation error while loading settings: {e}")
            raise ValueError(f"Invalid settings data: {e}") from e
        except Exception as e:
            LOGGER.error(f"Unexpected error loading settings: {e}")
            raise RuntimeError(f"Failed to load settings: {e}") from e

    def _load_settings_section(
        self, section_name: str, settings_data: dict[str, Any], model_class: type[BaseSettings]
    ) -> None:
        """Load a specific settings section using Pydantic v2 validation.

        This helper method streamlines the loading process by using Pydantic's
        automatic validation and type conversion. It pre-processes nested data
        to convert string quantities to Quantity objects before validation.

        Parameters
        ----------
        section_name : str
            Name of the settings section to load.
        settings_data : dict[str, Any]
            Complete settings data dictionary.
        model_class : type[BaseSettings]
            Pydantic model class to validate against.
        """
        if section_name in settings_data and hasattr(self, section_name):
            section_data = settings_data[section_name].copy()

            # Pre-process nested data to convert string quantities
            section_data = self._convert_quantities_recursive(section_data)

            # Let Pydantic handle all validation and type conversion automatically
            validated_model = model_class.model_validate(section_data)
            setattr(self, section_name, validated_model)

    def _load_fiber_settings(self, settings_data: dict[str, Any]) -> None:
        """Load fiber settings based on detected method or current configuration.

        This method detects the fiber method from the loaded data
        structure and uses the appropriate Pydantic model for validation.

        Parameters
        ----------
        settings_data : dict[str, Any]
            Settings data dictionary.
        """
        if "fibers" not in settings_data or not hasattr(self, "fibers"):
            return

        primary_method = "LSDYNA"
        fallback_method = "D-RBM"

        try:
            # Try loading with the primary method first.
            self._load_settings_section("fibers", settings_data, FibersBRBM)
            self._fiber_method = primary_method
            LOGGER.info(f"Successfully loaded fiber settings using {self._fiber_method} method")

            return self.fibers

        except (ValidationError, ValueError) as e:
            LOGGER.debug(f"Failed to load fiber settings {primary_method}. {e}")

        try:
            # Try the alternative method.
            self._load_settings_section("fibers", settings_data, FibersDRBM)
            self._fiber_method = fallback_method
            LOGGER.info(f"Successfully loaded fiber settings using {self._fiber_method} method")

            return self.fibers

        except (ValidationError, ValueError) as error:
            # Both methods failed - provide helpful error message
            raise ValueError(
                f"Failed to load fiber settings with both {primary_method} and "
                f"{fallback_method} methods. {error}."
            )

    def _convert_quantities_recursive(self, data: Any) -> Any:
        """Recursively convert string quantities to Quantity objects in nested data.

        This helper method processes nested dictionaries and lists to convert
        string representations of quantities to actual Quantity objects before
        Pydantic validation. This ensures proper handling of nested models.

        Parameters
        ----------
        data : Any
            Data structure to process (dict, list, or primitive value).

        Returns
        -------
        Any
            Processed data with string quantities converted to Quantity objects.
        """
        if isinstance(data, dict):
            return {key: self._convert_quantities_recursive(value) for key, value in data.items()}
        elif isinstance(data, list):
            return [self._convert_quantities_recursive(item) for item in data]
        elif isinstance(data, str):
            # Try to parse as a quantity if it looks like one
            if self._looks_like_quantity(data):
                try:
                    return ureg(data)
                except Exception:
                    # If parsing fails, return the original string
                    return data
            return data
        else:
            # Return primitive values unchanged
            return data

    def _looks_like_quantity(self, value: str) -> bool:
        """Check if a string looks like a quantity that can be parsed.

        Parameters
        ----------
        value : str
            String value to check.

        Returns
        -------
        bool
            True if the string appears to be a quantity representation.
        """
        # Simple heuristic: contains a space and has numeric part
        if " " not in value:
            return False

        parts = value.split()
        if len(parts) < 2:
            return False

        # Check if first part is numeric
        try:
            float(parts[0])
            return True
        except ValueError:
            return False

    def load_defaults(self):
        """Load default simulation settings using Pydantic model initialization.

        This method properly initializes all settings with default values using
        Pydantic's built-in validation and type conversion capabilities.

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
          end_time: 800.0 millisecond
          dtmin: 5.0 millisecond
          dtmax: 5.0 millisecond
          dt_d3plot: 20.0 millisecond
          dt_icvout: 5.0 millisecond
          global_damping: 0.1 / millisecond

        """
        try:
            # Load mechanics defaults using Pydantic model initialization
            if hasattr(self, "mechanics") and isinstance(self.mechanics, Mechanics):
                self.mechanics.analysis = Analysis(**mech_defaults.analysis)
                self.mechanics.boundary_conditions = BoundaryConditions(
                    **mech_defaults.boundary_conditions
                )
                self.mechanics.system = SystemModel(**mech_defaults.system_model)

            # Load zero pressure defaults
            if hasattr(self, "stress_free") and isinstance(self.stress_free, ZeroPressure):
                self.stress_free.analysis = AnalysisZeroPressure(**zero_pressure_defaults.analysis)

            # Load electrophysiology defaults
            if hasattr(self, "electrophysiology") and isinstance(
                self.electrophysiology, Electrophysiology
            ):
                self.electrophysiology.analysis = EPAnalysis(**ep_defaults.analysis)

                # Create stimulation dictionary with Pydantic validation
                stimulation_dict = {}
                for key, stim_data in ep_defaults.stimulation.items():
                    stimulation_dict[key] = Stimulation(**stim_data)
                self.electrophysiology.stimulation = stimulation_dict

            # Load Purkinje defaults
            if hasattr(self, "purkinje") and isinstance(self.purkinje, Purkinje):
                # Update Purkinje with defaults - handle properly based on Purkinje model structure
                for field_name, value in purkinje_defaults.build.items():
                    if hasattr(self.purkinje, field_name):
                        setattr(self.purkinje, field_name, value)

            # Load atrial fiber defaults
            if hasattr(self, "atrial_fibers") and isinstance(self.atrial_fibers, AtrialFiber):
                # Update atrial fibers with defaults - handle both la_bundle and ra_bundle
                for field_name, value in fibers_defaults.la_bundle.items():
                    if hasattr(self.atrial_fibers, field_name):
                        setattr(self.atrial_fibers, field_name, value)
                for field_name, value in fibers_defaults.ra_bundle.items():
                    if hasattr(self.atrial_fibers, field_name):
                        setattr(self.atrial_fibers, field_name, value)

        except Exception as e:
            LOGGER.error(f"Failed to load default settings: {e}")
            raise RuntimeError(f"Failed to initialize settings with defaults: {e}") from e

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
            if isinstance(getattr(self, attr), BaseSettings)
        ]

        for attr in attributes:
            if isinstance(attr, BaseSettings):
                attr.to_consistent_unit_system()
        return


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
    """Get consistent units formatted as string.

    Converts dimensionality to the PyAnsys Heart consistent unit system string
    representation based on the defined base quantities and derived units.

    Parameters
    ----------
    dimensions : set
        Set of dimensions from a Quantity object.

    Returns
    -------
    str
        String representation of consistent units for the given dimensions.
    """
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
    r"""Convert Windows path to WSL-compatible path format.

    Handles conversion from Windows drive paths and WSL localhost paths
    to proper Unix-style paths for use within Windows Subsystem for Linux.

    Parameters
    ----------
    windows_path : str
        Windows path to convert.

    Returns
    -------
    str | None
        WSL-compatible path string, or None if conversion not applicable.

    Examples
    --------
    >>> _windows_to_wsl_path(r"C:\Users\example")
    '/mnt/c/Users/example'
    >>> _windows_to_wsl_path(r"\\wsl.localhost\Ubuntu\home")
    '/Ubuntu/home'
    """
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
    """Class for collecting, managing, and validating LS-DYNA settings.

    This class provides configuration management for LS-DYNA simulations,
    including executable paths, parallelization settings, platform-specific
    configurations, and command-line argument generation.

    Parameters
    ----------
    lsdyna_path : pathlib.Path
        Path to LS-DYNA executable.
    dynatype : str
        Type of LS-DYNA executable (smp, intelmpi, platformmpi, msmpi).
    num_cpus : int
        Number of CPUs requested for parallel execution.
    platform : str
        Platform for LS-DYNA execution (windows, wsl, linux).
    dyna_options : str
        Additional command line options for LS-DYNA.
    mpi_options : str
        Additional MPI options for parallel execution.

    Examples
    --------
    >>> dyna_settings = DynaSettings(
    ...     lsdyna_path="lsdyna.exe", dynatype="intelmpi", num_cpus=4, platform="windows"
    ... )
    >>> commands = dyna_settings.get_commands("input.k")
    """

    @staticmethod
    def _get_available_mpi_exe():
        """Find whether mpiexec or mpirun are available.

        Searches for MPI executables in PATH, preferring mpirun over mpiexec.

        Returns
        -------
        str
            Path to available MPI executable.

        Raises
        ------
        MPIProgamNotFoundError
            If neither mpirun nor mpiexec are found in PATH.
        """
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
        num_cpus : int, default: 1
            Number of CPUs requested.
        platform : Literal["windows", "wsl", "linux"], default: "windows"
            Platform.
        dyna_options : str, default: ""
            Additional command line options.
        mpi_options : str, default: ""
            Additional MPI options.
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

    def get_commands(self, path_to_input: pathlib.Path) -> list[str]:
        """Get command line arguments from the defined settings.

        Builds platform-specific command line arguments for running LS-DYNA
        with the configured settings including MPI and parallelization options.

        Parameters
        ----------
        path_to_input : pathlib.Path
            Path to the LS-DYNA input file.

        Returns
        -------
        list[str]
            List of command line arguments for executing LS-DYNA.

        Raises
        ------
        WSLNotFoundError
            If WSL platform is specified but wsl.exe is not found.

        Examples
        --------
        >>> dyna_settings = DynaSettings(dynatype="smp", num_cpus=4)
        >>> commands = dyna_settings.get_commands(Path("input.k"))
        >>> print(commands[0])  # LS-DYNA executable path
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
        """Set DynaSettings based on globally defined settings for PyAnsys-Heart.

        Checks for PYANSYS_HEART environment variables and updates settings
        accordingly. Supported environment variables:
        - PYANSYS_HEART_LSDYNA_PATH: Path to LS-DYNA executable
        - PYANSYS_HEART_LSDYNA_PLATFORM: Execution platform
        - PYANSYS_HEART_LSDYNA_TYPE: LS-DYNA executable type
        - PYANSYS_HEART_NUM_CPU: Number of CPUs for parallel execution
        """
        keys = [key for key in os.environ.keys() if "PYANSYS_HEART" in key]
        LOGGER.debug(f"PYANSYS_HEART Environment variables: {keys}")
        self.lsdyna_path = os.getenv("PYANSYS_HEART_LSDYNA_PATH", self.lsdyna_path)
        self.platform = os.getenv("PYANSYS_HEART_LSDYNA_PLATFORM", self.platform)
        self.dynatype = os.getenv("PYANSYS_HEART_LSDYNA_TYPE", self.dynatype)
        self.num_cpus = int(os.getenv("PYANSYS_HEART_NUM_CPU", self.num_cpus))
        return

    def __repr__(self):
        """Represent self as YAML-formatted string.

        Returns
        -------
        str
            YAML representation of the DynaSettings object attributes.
        """
        return yaml.dump(vars(self), allow_unicode=True, default_flow_style=False)
