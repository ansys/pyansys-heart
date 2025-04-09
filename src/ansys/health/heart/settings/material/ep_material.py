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
from typing import Optional

from ansys.health.heart.settings.defaults import electrophysiology as ep_defaults


@dataclass
class CellModel:
    """Abstract class for different cell models."""

    @dataclass
    class Tentusscher:
        """Hold data for Tentusscher cell model."""

        gas_constant: float = 8314.472
        t: float = 310
        faraday_constant: float = 96485.3415
        cm: float = 0.185
        vc: float = 0.016404
        vsr: float = 0.001094
        vss: float = 0.00005468
        pkna: float = 0.03
        ko: float = 5.4
        nao: float = 140.0
        cao: float = 2.0
        gk1: float = 5.405
        gkr: float = 0.153
        gna: float = 14.838
        gbna: float = 0.0002
        gcal: float = 0.0000398
        gbca: float = 0.000592
        gpca: float = 0.1238
        gpk: float = 0.0146
        pnak: float = 2.724
        km: float = 1.0
        kmna: float = 40.0
        knaca: float = 1000.0
        ksat: float = 0.1
        alpha: float = 2.5
        gamma: float = 0.35
        kmca: float = 1.38
        kmnai: float = 87.5
        kpca: float = 0.0005
        k1: float = 0.15
        k2: float = 0.045
        k3: float = 0.06
        k4: float = 0.005
        ec: float = 1.5
        maxsr: float = 2.5
        minsr: float = 1.0
        vrel: float = 0.102
        vleak: float = 0.00036
        vxfer: float = 0.0038
        vmaxup: float = 0.006375
        kup: float = 0.00025
        bufc: float = 0.2
        kbufc: float = 0.001
        bufsr: float = 10.0
        kbufsf: float = 0.3
        bufss: float = 0.4
        kbufss: float = 0.00025
        # gas_constant=8314.472,
        # faraday_constant=96485.3415,
        gks: float = 0.392
        gto: float = 0.294
        v: float = -85.23
        ki: float = 136.89
        nai: float = 8.604
        cai: float = 0.000126
        cass: float = 0.00036
        casr: float = 3.64
        rpri: float = 0.9073
        xr1: float = 0.00621
        xr2: float = 0.4712
        xs: float = 0.0095
        m: float = 0.00172
        h: float = 0.7444
        j: float = 0.7045
        d: float = 3.373e-5
        f: float = 0.7888
        f2: float = 0.9755
        fcass: float = 0.9953
        s: float = 0.999998
        r: float = 2.42e-8

        def to_dictionary(self):
            """Convert to dictionary."""
            return self.__dict__

        pass

    @dataclass
    class TentusscherEndo(Tentusscher):
        """Hold data for Tentusscher cell model in its endocardium version."""

        gks: float = 0.392
        gto: float = 0.073
        v: float = -86.709
        ki: float = 138.4
        nai: float = 10.355
        cai: float = 0.00013
        cass: float = 0.00036
        casr: float = 3.715
        rpri: float = 0.9068
        xr1: float = 0.00448
        xr2: float = 0.476
        xs: float = 0.0087
        m: float = 0.00155
        h: float = 0.7573
        j: float = 0.7225
        d: float = 3.164e-5
        f: float = 0.8009
        f2: float = 0.9778
        fcass: float = 0.9953
        s: float = 0.3212
        r: float = 2.235e-8

    @dataclass
    class TentusscherEpi(Tentusscher):
        """Hold data for Tentusscher cell model in its epicardium version."""

        gks: float = 0.392
        gto: float = 0.294
        v: float = -85.23
        ki: float = 136.89
        nai: float = 8.604
        cai: float = 0.000126
        cass: float = 0.00036
        casr: float = 3.64
        rpri: float = 0.9073
        xr1: float = 0.00621
        xr2: float = 0.4712
        xs: float = 0.0095
        m: float = 0.00172
        h: float = 0.7444
        j: float = 0.7045
        d: float = 3.373e-5
        f: float = 0.7888
        f2: float = 0.9755
        fcass: float = 0.9953
        s: float = 0.999998
        r: float = 2.42e-8

    @dataclass
    class TentusscherMid(Tentusscher):
        """Hold data for Tentusscher cell model in its mid-myocardium version."""

        gks: float = 0.098
        gto: float = 0.294
        v: float = -85.423
        ki: float = 138.52
        nai: float = 10.132
        cai: float = 0.000153
        cass: float = 0.00042
        casr: float = 4.272
        rpri: float = 0.8978
        xr1: float = 0.0165
        xr2: float = 0.473
        xs: float = 0.0174
        m: float = 0.00165
        h: float = 0.749
        j: float = 0.6788
        d: float = 3.288e-5
        f: float = 0.7026
        f2: float = 0.9526
        fcass: float = 0.9942
        s: float = 0.999998
        r: float = 2.347e-8


@dataclass
class EPMaterialModel:
    """Base class for all EP material models."""

    sigma_fiber: float = None
    sigma_sheet: Optional[float] = None
    sigma_sheet_normal: Optional[float] = None
    beta: Optional[float] = ep_defaults.material["myocardium"]["beta"].m
    cm: Optional[float] = ep_defaults.material["myocardium"]["cm"].m
    lambda_: Optional[float] = None

    def __post_init__(self):
        """Post init method."""
        if self.sigma_sheet is not None and self.sigma_sheet_normal is None:
            self.sigma_sheet_normal = self.sigma_sheet
        if self.sigma_sheet_normal is not None and self.sigma_sheet is None:
            self.sigma_sheet = self.sigma_sheet_normal


@dataclass
class EPMaterial(EPMaterialModel):
    """EP material module."""

    pass

    @dataclass
    class Insulator:
        """Insulator material."""

        sigma_fiber = 0
        cm = 0
        beta = 0
        pass

        def __repr__(self):
            """Print a message."""
            return "Insulator."

    @dataclass
    class Active(EPMaterialModel):
        """Hold data for EP material."""

        solver_type = ep_defaults.analysis["solvertype"]
        if solver_type == "Monodomain":
            sig1 = ep_defaults.material["myocardium"]["sigma_fiber"].m
            sig2 = ep_defaults.material["myocardium"]["sigma_sheet"].m
            sig3 = ep_defaults.material["myocardium"]["sigma_sheet_normal"].m
        elif solver_type == "Eikonal" or solver_type == "ReactionEikonal":
            sig1 = ep_defaults.material["myocardium"]["velocity_fiber"].m
            sig2 = ep_defaults.material["myocardium"]["velocity_sheet"].m
            sig3 = ep_defaults.material["myocardium"]["velocity_sheet_normal"].m

        sigma_fiber: float = sig1
        sigma_sheet: Optional[float] = sig2
        sigma_sheet_normal: Optional[float] = sig3
        cell_model: CellModel = field(default_factory=lambda: CellModel.Tentusscher())

    class ActiveBeam(Active):
        """Hold data for beam active EP material."""

        sigma_fiber = ep_defaults.material["beam"]["sigma"].m
        cell_model = field(default_factory=lambda: CellModel.TentusscherEndo())
        pmjres: float = ep_defaults.material["beam"]["pmjres"].m

    class Passive(EPMaterialModel):
        """Hold data for EP passive material."""

        sigma_fiber: float = ep_defaults.material["myocardium"]["sigma_fiber"].m
        sigma_sheet: Optional[float] = None
        sigma_sheet_normal: Optional[float] = None

    @dataclass
    class DummyMaterial:
        """Just for initialization."""

        pass

        def __repr__(self):
            """Print a message."""
            return "Material is empty."


if __name__ == "__main__":
    m1 = EPMaterial.DummyMaterial()

    tentus1 = CellModel.Tentusscher()
    tentus2 = CellModel.TentusscherEndo()

    mat1 = EPMaterial.Insulator()
    mat2 = EPMaterial.Passive(sigma_fiber=0, sigma_sheet=2)
    mat3 = EPMaterial.Active(sigma_fiber=1, sigma_sheet=2, sigma_sheet_normal=3)
    mat4 = EPMaterial.Active(sigma_fiber=1, sigma_sheet=2, sigma_sheet_normal=3)
    print("done")
