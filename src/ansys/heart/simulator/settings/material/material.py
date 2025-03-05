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

from dataclasses import dataclass, field
from typing import List, Optional

from ansys.heart.core import LOG as LOGGER
from ansys.heart.simulator.settings.material.curve import ActiveCurve, constant_ca2


@dataclass
class ISO:
    """Isotropic module of MAT_295."""

    itype: int = -3
    """Type of isotropic model."""
    beta: float = 0.0
    """Volumetric response function coefficient."""
    nu: float = 0.499
    """Poisson's coefficient."""
    k1: float = None
    """Holzapfel-Ogden modulus."""
    k2: float = None
    """Holzapfel-Ogden constant."""
    mu1: float = None
    """Ogden modulus."""
    alpha1: float = None
    """Odgen constant."""
    kappa: float = None
    """Bulk modulus."""

    def __post_init__(self):
        """Test inputs."""
        if self.k1 is not None and self.k2 is not None:
            if not abs(self.itype) == 3:
                # must be HGO model
                raise ValueError("ISO input is invalid, abs(itype)!=3.")
        elif self.mu1 is not None and self.alpha1 is not None:
            if not abs(self.itype) == 1:
                # must be Odgen model
                raise ValueError("ISO input is invalid, abs(itype)!=1.")
        else:
            raise ValueError("ISO input is not supported.")

        if self.kappa is not None:
            # replace Poisson's coefficient
            mu = self.k1 if abs(self.itype) == 3 else self.mu1
            self.nu = (3 * self.kappa - 2 * mu) / (6 * self.kappa + 2 * mu)


@dataclass
class ANISO:
    """Anisotropic module of MAT_295."""

    @dataclass
    class HGOFiber:
        """Define HGO type fiber from k1 and k2."""

        k1: float = None
        """Holzapfel-Ogden modulus of fiber."""
        k2: float = None
        """Holzapfel-Ogden constant of fiber."""
        a: float = 0.0
        """First dispersion parameter."""
        b: float = 1.0
        """Second dispersion parameter."""
        _theta: float = None
        """Second fiber rotation angle, 90."""
        _ftype: int = 1
        """Fiber type."""
        _fcid: int = 0
        """Curve ID, keep it 0."""

    atype: int = -1
    "Type of anisotropic model."
    fibers: List[HGOFiber] = None
    """List of fibers."""

    k1fs: Optional[float] = None
    """Coupling modulus between the fiber and sheet directions"""
    k2fs: Optional[float] = None
    """Coupling constant between the fiber and sheet directions"""

    vec_a: tuple = (1.0, 0.0, 0.0)
    """Components of vector a, don't change."""
    vec_d: tuple = (0.0, 1.0, 0.0)
    """Components of vector d, don't change."""

    def __post_init__(self):
        """Check and deduce other parameters from input."""
        # create a default one if not given
        if self.fibers is None:
            self.fibers = [self.HGOFiber()]

        # check if legal
        if len(self.fibers) != 1 and len(self.fibers) != 2:
            LOGGER.error("The number of fibers must be 1 or 2.")
            exit()

        # deduce input
        self.nf = len(self.fibers)

        if self.k1fs is not None and self.k2fs is not None:
            if len(self.fibers) == 2:
                self.intype = 1
            else:
                LOGGER.error("A single fiber cannot have interaction.")
                exit()
        else:
            self.intype = 0

        self.fibers[0]._theta = 0.0
        if self.nf > 1:
            self.fibers[1]._theta = 90.0

    def __repr__(self):
        """Ensure the print output contains fields from __post_init__."""
        attrs = ", ".join(f"{attr}={getattr(self, attr)}" for attr in self.__annotations__)
        attrs += f", nf={self.nf}, intype={self.intype}"
        return f"{self.__class__.__name__}({attrs})"


@dataclass
class ActiveModel:
    """Abstract class for different active models."""

    pass

    @dataclass
    class Model1:
        """Hold data for active model 1."""

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

    @dataclass
    class Model3:
        """Hold data for active model 3."""

        t0: float = None
        ca2ion50: float = 1.0
        n: float = 1.0
        f: float = 0.0
        l: float = 1.0  # no effect if eta=0 #noqa: E741
        eta: float = 0.0
        sigmax: float = None


@dataclass
class ACTIVE:
    """Active module of MAT_295."""

    acid: int = None  # empty for ep_coupled, or curve ID from writer
    actype: int = None  # defined in __post_init__
    acthr: float = (
        None  # need to be defined for ep_coupled, for mechanics it's defined in ActiveCurve
    )
    acdir: int = 1  # always act in fiber direction
    sf: float = 1.0  # always 1.0 and controls contractility in ActiveModel
    ss: float = 0.0
    sn: float = 0.0
    model: ActiveModel = field(default_factory=ActiveModel.Model1)
    ca2_curve: ActiveCurve = field(
        default_factory=lambda: ActiveCurve(constant_ca2(), threshold=0.1, type="ca2")
    )

    def __post_init__(self):
        """Deduce actype."""
        if isinstance(self.model, ActiveModel.Model1):
            self.actype = 1
        elif isinstance(self.model, ActiveModel.Model3):
            self.actype = 3
        else:
            LOGGER.error("Unknown actype.")
        if self.ca2_curve is not None:
            self.acthr = self.ca2_curve.threshold


@dataclass
class MechanicalMaterialModel:
    """Base class for mechanical material models."""

    pass

    @dataclass
    class DummyMaterial:
        """Just for initialization."""

        pass

        def __repr__(self):
            """Print a message."""
            return "Material is empty."


@dataclass
class MAT295(MechanicalMaterialModel):
    """MAT_295, check LS-DYNA manul for more details."""

    rho: float
    """Mass density."""
    iso: ISO
    """Isotropic module."""
    aopt: float = 2.0
    """Material axes option, dont' change."""
    aniso: Optional[ANISO] = None
    """Anisotropic module."""
    active: Optional[ACTIVE] = None
    """Active module."""


@dataclass
class NeoHookean(MechanicalMaterialModel):
    """NeoHookean model, passive isotropic material."""

    rho: float
    """Mass density."""
    c10: float  # mu/2
    """Neohookean parameter, is half of shear modulus."""
    kappa: float
    """Bulk modulus."""
    nu: float = None
    """Poisson's coefficient."""

    def __post_init__(self):
        """Deduce Poisson's coefficient."""
        if self.kappa is not None:
            # replace Poisson's coefficient
            mu = self.c10 / 2.0
            self.nu = (3 * self.kappa - 2 * mu) / (6 * self.kappa + 2 * mu)
