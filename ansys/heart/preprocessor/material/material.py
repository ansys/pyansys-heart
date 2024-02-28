"""Material module."""

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
    fiber: List[HGO_Fiber] = None

    k1fs: Optional[float] = None
    k2fs: Optional[float] = None

    vec_a: tuple = (1.0, 0.0, 0.0)
    vec_d: tuple = (0.0, 1.0, 0.0)

    def __post_init__(self):
        """Check and deduce other parameters from input."""
        # create a default one if not given
        if self.fiber is None:
            self.fiber = [self.HGO_Fiber()]

        # check if legal
        if len(self.fiber) != 1 and len(self.fiber) != 2:
            print("No. of fiber must be 1 or 2.")
            exit()

        # deduce input
        self.nf = len(self.fiber)

        if self.k1fs is not None and self.k2fs is not None:
            if len(self.fiber) == 2:
                self.intype = 1
            else:
                print("One fiber cannot have interaction.")
                exit()
        else:
            self.intype = 0

        self.fiber[0]._theta = 0.0
        if self.nf > 1:
            self.fiber[1]._theta = 90.0

    def __repr__(self):
        """Make sure print contains field in __post_init__."""
        attrs = ", ".join(f"{attr}={getattr(self, attr)}" for attr in self.__annotations__)
        attrs += f", nf={self.nf}, intype={self.intype}"
        return f"{self.__class__.__name__}({attrs})"


@dataclass
class ActiveModel:
    """Abstract class for active models."""

    pass


@dataclass
class ActiveModel1(ActiveModel):
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
class ActiveModel3(ActiveModel):
    """Hold data for active model 3."""

    # TODO
    t0: float = None


@dataclass
class ACTIVE:
    """Active module of *mat295."""

    acdir: int = 1
    acid: int = None
    acthr: float = None
    sf: float = 1.0
    ss: float = 0.0
    sn: float = 0.0
    model: ActiveModel = ActiveModel1()
    ca2_curve: ActiveCurve = ActiveCurve(constant_ca2(), threshold=0.1, type="ca2")

    def __post_init__(self):
        """Deduce actype."""
        if isinstance(self.model, ActiveModel1):
            self.actype = 1
        elif isinstance(self.model, ActiveModel3):
            self.actype = 3
        else:
            print("Unknown actype.")

        self.acthr = self.ca2_curve.threshold


@dataclass
class MAT295:
    """Hold data for *mat295."""

    rho: float = 0.001
    aopt: float = 2.0

    iso: ISO = ISO()
    aniso: Optional[ANISO] = ANISO()
    active: Optional[ACTIVE] = ACTIVE()


class NeoHookean:
    """Passive isotropic material."""

    rho: float
    c10: float  # mu/2
    nnu: float = 0.499


import dataclasses

from ansys.dyna.keywords import keywords
import pandas as pd


class _MaterialHGOMyocardium(keywords.Mat295):
    def __init__(self, id: int, mat: MAT295, ignore_active: bool = False):
        # 1st line
        super().__init__(mid=id)
        setattr(self, "rho", mat.rho)
        setattr(self, "aopt", mat.aopt)

        # iso
        for field in dataclasses.fields(mat.iso):
            value = getattr(mat.iso, field.name)
            setattr(self, field.name, value)

        # aniso
        if mat.aniso is not None:
            self.atype = mat.aniso.atype
            self.intype = mat.aniso.intype
            self.nf = mat.aniso.nf
            self.ftype = mat.aniso.fiber[0]._ftype  # not used but must be defined

            self.a1 = mat.aniso.vec_a[0]
            self.a2 = mat.aniso.vec_a[1]
            self.a3 = mat.aniso.vec_a[2]

            self.d1 = mat.aniso.vec_d[0]
            self.d2 = mat.aniso.vec_d[1]
            self.d3 = mat.aniso.vec_d[2]

            fiber_sheet = []
            for i in range(len(mat.aniso.fiber)):
                dct = {
                    "theta": mat.aniso.fiber[i]._theta,
                    "a": mat.aniso.fiber[i].a,
                    "b": mat.aniso.fiber[i].b,
                    "ftype": mat.aniso.fiber[i]._ftype,
                    "fcid": mat.aniso.fiber[i]._fcid,
                    "k1": mat.aniso.fiber[i].k1,
                    "k2": mat.aniso.fiber[i].k2,
                }
                fiber_sheet.append(dct)
            self.anisotropic_settings = pd.DataFrame(fiber_sheet)

            if mat.aniso.intype == 1:
                self.coupling_k1 = mat.aniso.k1fs
                self.coupling_k2 = mat.aniso.k2fs

        # active
        if not ignore_active and mat.active is not None:
            setattr(self, "actype", mat.active.actype)

            for field in dataclasses.fields(mat.active):
                if field.type is ActiveModel:  # nested data of active model
                    for nested_f in dataclasses.fields(mat.active.model):
                        name = nested_f.name
                        value = getattr(mat.active.model, name)
                        setattr(self, name, value)
                else:
                    # acdir, acid ....
                    name = field.name
                    value = getattr(mat.active, name)
                    setattr(self, name, value)


m = MAT295()
print(m)
kw = _MaterialHGOMyocardium(1, m)
print(kw)
