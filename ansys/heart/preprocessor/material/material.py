"""Material module."""

import logging

LOGGER = logging.getLogger("pyheart_global.material")

from dataclasses import dataclass
from typing import List, Optional

from ansys.heart.preprocessor.material.curve import ActiveCurve, constant_ca2


@dataclass
class ISO:
    """Isotropic module of *mat295."""

    itype: int = -3
    beta: float = 0.0
    nu: float = 0.499
    k1: float = None
    k2: float = None


@dataclass
class ANISO:
    """Anisotropic module of *mat295."""

    @dataclass
    class HGO_Fiber:
        """Define HGO type fiber from k1 and k2."""

        k1: float = None
        k2: float = None
        a: float = 0.0
        b: float = 1.0
        _theta: float = None
        _ftype: int = 1
        _fcid: int = 0

    atype: int = -1
    fibers: List[HGO_Fiber] = None

    k1fs: Optional[float] = None
    k2fs: Optional[float] = None

    vec_a: tuple = (1.0, 0.0, 0.0)
    vec_d: tuple = (0.0, 1.0, 0.0)

    def __post_init__(self):
        """Check and deduce other parameters from input."""
        # create a default one if not given
        if self.fibers is None:
            self.fibers = [self.HGO_Fiber()]

        # check if legal
        if len(self.fibers) != 1 and len(self.fibers) != 2:
            LOGGER.error("No. of fiber must be 1 or 2.")
            exit()

        # deduce input
        self.nf = len(self.fibers)

        if self.k1fs is not None and self.k2fs is not None:
            if len(self.fibers) == 2:
                self.intype = 1
            else:
                LOGGER.error("One fiber cannot have interaction.")
                exit()
        else:
            self.intype = 0

        self.fibers[0]._theta = 0.0
        if self.nf > 1:
            self.fibers[1]._theta = 90.0

    def __repr__(self):
        """Make sure print contains field in __post_init__."""
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
        l: float = 1.85
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
        l: float = 1.0  # no effect if eta=0
        eta: float = 0.0
        sigmax: float = None


@dataclass
class ACTIVE:
    """Active module of *mat295."""

    acid: int = None  # not used
    actype: int = None  # defined by curve
    acthr: float = None  # defined by curve
    acdir: int = 1
    sf: float = 1.0
    ss: float = 0.0
    sn: float = 0.0
    model: ActiveModel = ActiveModel.Model1()
    ca2_curve: ActiveCurve = ActiveCurve(constant_ca2(), threshold=0.1, type="ca2")

    def __post_init__(self):
        """Deduce actype."""
        if isinstance(self.model, ActiveModel.Model1):
            self.actype = 1
        elif isinstance(self.model, ActiveModel.Model3):
            self.actype = 3
        else:
            LOGGER.error("Unknown actype.")

        self.acthr = self.ca2_curve.threshold


@dataclass
class MechaMaterialModel:
    """Base class for all mechanical material model."""

    pass

    @dataclass
    class DummyMaterial:
        """Just for initialization."""

        pass

        def __repr__(self):
            """Print a message."""
            return "Material is empty."


@dataclass
class MAT295(MechaMaterialModel):
    """Hold data for *mat295."""

    rho: float
    iso: ISO
    aopt: float = 2.0
    aniso: Optional[ANISO] = None
    active: Optional[ACTIVE] = None


@dataclass
class NeoHookean(MechaMaterialModel):
    """Passive isotropic material."""

    rho: float
    c10: float  # mu/2
    nu: float = 0.499
