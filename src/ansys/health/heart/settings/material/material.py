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

"""Material module."""

from typing import List, Literal, Optional

from deprecated import deprecated
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.settings.material.curve import ActiveCurve, constant_ca2


class ISO(BaseModel):
    """Isotropic module of MAT_295."""

    itype: Literal[-3, -1, 1, 3] = -3
    """Isotropic material type."""
    """+/-3: HGO model, +/-1: Ogden model."""
    beta: float = 0.0
    """Volumetric response coefficient."""
    nu: float = Field(default=0.499, le=0.5, ge=0.0)
    """Possion's ratio."""
    k1: Optional[float] = None
    """k1 for HGO model."""
    k2: Optional[float] = None
    """k2 for HGO model."""
    mu1: Optional[float] = None
    """mu1 for Ogden model."""
    alpha1: Optional[float] = None
    """alpha1 for Ogden model."""
    kappa: Optional[float] = None
    """Bulk modulus."""

    @model_validator(mode="after")
    def validate_material_parameters(self):
        """Validate that required parameters are provided based on itype."""
        # Check HGO model requirements (itype = ±3)
        if abs(self.itype) == 3:
            if self.k1 is None or self.k2 is None:
                raise ValueError(
                    f"For HGO model (itype = ±3), both k1 and k2 must be specified. "
                    f"Got k1={self.k1}, k2={self.k2}"
                )
            # Ensure Ogden parameters are not provided
            if self.mu1 is not None or self.alpha1 is not None:
                raise ValueError(
                    f"For HGO model (itype = {self.itype}), mu1 and alpha1 should not be specified."
                    f"Got mu1={self.mu1}, alpha1={self.alpha1}"
                )

        # Check Ogden model requirements (itype = ±1)
        elif abs(self.itype) == 1:
            if self.mu1 is None or self.alpha1 is None:
                raise ValueError(
                    f"For Ogden model (itype = ±1), both mu1 and alpha1 must be specified. "
                    f"Got mu1={self.mu1}, alpha1={self.alpha1}"
                )
            # Ensure HGO parameters are not provided
            if self.k1 is not None or self.k2 is not None:
                raise ValueError(
                    f"For Ogden model (itype = {self.itype}), k1 and k2 should not be specified. "
                    f"Got k1={self.k1}, k2={self.k2}"
                )

        # Update Poisson's ratio if kappa is provided
        if self.kappa is not None:
            mu = self.k1 if abs(self.itype) == 3 else self.mu1
            if mu is None:
                raise ValueError(
                    "Cannot calculate nu from kappa: material stiffness parameter is None"
                )
            self.nu = (3 * self.kappa - 2 * mu) / (6 * self.kappa + 2 * mu)

        # Warning for low Poisson's ratio
        if self.nu < 0.49:
            LOGGER.warning("Poisson's ratio lower than 0.49 is not recommended.")

        return self


class HGOFiber(BaseModel):
    """Define HGO type fiber from k1 and k2."""

    k1: Optional[float] = None
    """k1 for HGO model."""
    k2: Optional[float] = None
    """k2 for HGO model."""
    a: float = 0.0
    """Fiber dispersion tensor parameter."""
    b: float = 1.0
    """Fiber dispersion tensor parameter."""
    _theta: float = None
    """0 for fiber 1, 90 for fiber 2. Don't change it."""
    _ftype: int = 1
    """Fiber type, 1 for HGO."""
    _fcid: int = 0
    """Not used yet."""


class ANISO(BaseModel):
    """Anisotropic module of MAT_295."""

    atype: int = -1
    """Type of anisotropic model."""
    fibers: List[HGOFiber] = Field(default_factory=lambda: [HGOFiber()])
    """List of fibers."""

    k1fs: Optional[float] = None
    """k1 for HGO model for coupling between fibers."""
    k2fs: Optional[float] = None
    """k2 for HGO model for coupling between fibers."""

    vec_a: tuple = (1.0, 0.0, 0.0)
    """Component of fiber direction, don't change it."""
    vec_d: tuple = (0.0, 1.0, 0.0)
    """Component of sheet direction, don't change it."""

    _nf = None
    """Number of fibers."""
    _intype = None
    """Interaction type."""

    @model_validator(mode="after")
    def validate_material_parameters(self):
        """Check and deduce other parameters from input."""
        # check if legal
        if len(self.fibers) != 1 and len(self.fibers) != 2:
            LOGGER.error("No. of fiber must be 1 or 2.")
            raise ValueError("No. of fibers must be 1 or 2.")

        # deduce input
        self._nf = len(self.fibers)

        if self.k1fs is not None and self.k2fs is not None:
            if len(self.fibers) == 2:
                self._intype = 1
            else:
                LOGGER.error("One fiber cannot have an interaction term.")
                raise ValueError("One fiber cannot have an interaction term.")
        else:
            self._intype = 0

        self.fibers[0]._theta = 0.0
        if self._nf > 1:
            self.fibers[1]._theta = 90.0

        return self

    def __repr__(self):
        """Make sure print contains field in __post_init__."""
        attrs = ", ".join(f"{attr}={getattr(self, attr)}" for attr in self.__annotations__)
        attrs += f", nf={self._nf}, intype={self._intype}"
        return f"{self.__class__.__name__}({attrs})"


class ActiveModel(BaseModel):
    """Active model base class."""


class ActiveModel1(ActiveModel):
    """Hold data for active model 1, check manual for details."""

    t0: float = None
    ca2ion: float = None
    ca2ionm: float = 4.35
    n: int = 2
    taumax: float = 0.125
    stf: float = 0.0
    b: float = 4.75
    l0: float = 1.58
    l: float = 1.85  # noqa: E741
    dtmax: float = 150
    mr: float = 1048.9
    tr: float = -1629.0


class ActiveModel3(ActiveModel):
    """Hold data for active model 3, check manual for details."""

    t0: float = None
    ca2ion50: float = 1.0
    n: float = 1.0
    f: float = 0.0
    l: float = 1.0  # no effect if eta=0 #noqa: E741
    eta: float = 0.0
    sigmax: float = None


class ACTIVE(BaseModel):
    """Active module of MAT_295."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    acid: Optional[int] = None  # empty for ep_coupled, or curve ID from writer
    """Do not define it, it will be assigned with an ID of Ca2+ curve
    for mechanical problem or empty for ep-coupled problem."""
    actype: Optional[int] = None  # defined in __post_init__
    """Type of active model, will be deduced the model validator."""
    acthr: Optional[float] = None
    """Ca2+ threshold for active stress, need to be defined for ep-coupled,
    for mechanics it's defined in ActiveCurve."""
    acdir: int = 1  # always act in fiber direction
    """Direction of active stress, don't change it."""
    sf: float = 1.0  # always 1.0 and controls contractility in ActiveModel
    """Scaling factor on fiber direction."""
    ss: float = 0.0
    """Scaling factor on sheet direction."""
    sn: float = 0.0
    """Scaling factor on normal direction."""
    model: ActiveModel = ActiveModel1()
    """Active model."""

    ca2_curve: Optional[ActiveCurve] = ActiveCurve(func=constant_ca2(), threshold=0.1, type="ca2")
    """Ca2+ curve for mechanical problem."""

    @model_validator(mode="after")
    def validate_material_parameters(self):  # noqa: D102
        if isinstance(self.model, ActiveModel1):
            self.actype = 1
        elif isinstance(self.model, ActiveModel3):
            self.actype = 3
        else:
            LOGGER.error("Unknown actype.")
        if self.ca2_curve is not None:
            self.acthr = self.ca2_curve.threshold

        return self


class MechanicalMaterialModel(BaseModel):
    """Base class for all mechanical material model."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    pass


class Mat295(MechanicalMaterialModel):
    """Hold data for MAT_ANISOTROPIC_HYPERELASTIC (MAT_295)."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    rho: float
    """Density of the material."""
    iso: ISO
    """Isotropic module."""
    aopt: float = 2.0
    """Matrerial axe option, don't change it."""
    aniso: Optional[ANISO] = None
    """Anisotropic module."""
    active: Optional[ACTIVE] = None
    """Active module."""

    @field_validator("aniso", mode="before")
    def check_aniso(cls, v):  # noqa D102
        """Check if aniso is provided."""
        if v is None or isinstance(v, (ANISO, dict)):
            return v
        else:
            raise TypeError("aniso must be of type ANISO or dict.")


@deprecated(reason="Use *MAT_295 with the ISO module instead.")
class NeoHookean(MechanicalMaterialModel):
    """Passive isotropic material with MAT_77H."""

    rho: float
    """Density of the material."""
    c10: float  # mu/2
    """c10."""
    kappa: Optional[float] = None
    """Bulk modulus."""
    nu: Optional[float] = None
    """Poisson's ratio."""

    @model_validator(mode="after")
    def validate_material_parameters(self):  # noqa: D102
        """Deduce Poisson's ratio if not given."""
        if self.kappa is not None:
            # replace Poisson's coefficient
            mu = self.c10 * 2
            self.nu = (3 * self.kappa - 2 * mu) / (6 * self.kappa + 2 * mu)
        if self.nu is not None and self.nu < 0.49:
            LOGGER.warning("Poisson's ratio lower than 0.49 is not recommended.")
