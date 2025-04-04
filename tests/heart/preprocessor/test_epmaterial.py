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

from ansys.heart.core.settings.material.ep_material import CellModel, EPMaterial

midcelldata = dict(
    [
        ("gas_constant", 8314.472),
        ("t", 310),
        ("faraday_constant", 96485.3415),
        ("cm", 0.185),
        ("vc", 0.016404),
        ("vsr", 0.001094),
        ("vss", 5.468e-05),
        ("pkna", 0.03),
        ("ko", 5.4),
        ("nao", 140.0),
        ("cao", 2.0),
        ("gk1", 5.405),
        ("gkr", 0.153),
        ("gna", 14.838),
        ("gbna", 0.0002),
        ("gcal", 3.98e-05),
        ("gbca", 0.000592),
        ("gpca", 0.1238),
        ("gpk", 0.0146),
        ("pnak", 2.724),
        ("km", 1.0),
        ("kmna", 40.0),
        ("knaca", 1000.0),
        ("ksat", 0.1),
        ("alpha", 2.5),
        ("gamma", 0.35),
        ("kmca", 1.38),
        ("kmnai", 87.5),
        ("kpca", 0.0005),
        ("k1", 0.15),
        ("k2", 0.045),
        ("k3", 0.06),
        ("k4", 0.005),
        ("ec", 1.5),
        ("maxsr", 2.5),
        ("minsr", 1.0),
        ("vrel", 0.102),
        ("vleak", 0.00036),
        ("vxfer", 0.0038),
        ("vmaxup", 0.006375),
        ("kup", 0.00025),
        ("bufc", 0.2),
        ("kbufc", 0.001),
        ("bufsr", 10.0),
        ("kbufsf", 0.3),
        ("bufss", 0.4),
        ("kbufss", 0.00025),
        ("gks", 0.098),
        ("gto", 0.294),
        ("v", -85.423),
        ("ki", 138.52),
        ("nai", 10.132),
        ("cai", 0.000153),
        ("cass", 0.00042),
        ("casr", 4.272),
        ("rpri", 0.8978),
        ("xr1", 0.0165),
        ("xr2", 0.473),
        ("xs", 0.0174),
        ("m", 0.00165),
        ("h", 0.749),
        ("j", 0.6788),
        ("d", 3.288e-05),
        ("f", 0.7026),
        ("f2", 0.9526),
        ("fcass", 0.9942),
        ("s", 0.999998),
        ("r", 2.347e-08),
    ]
)

endocelldata = dict(
    [
        ("gas_constant", 8314.472),
        ("t", 310),
        ("faraday_constant", 96485.3415),
        ("cm", 0.185),
        ("vc", 0.016404),
        ("vsr", 0.001094),
        ("vss", 5.468e-05),
        ("pkna", 0.03),
        ("ko", 5.4),
        ("nao", 140.0),
        ("cao", 2.0),
        ("gk1", 5.405),
        ("gkr", 0.153),
        ("gna", 14.838),
        ("gbna", 0.0002),
        ("gcal", 3.98e-05),
        ("gbca", 0.000592),
        ("gpca", 0.1238),
        ("gpk", 0.0146),
        ("pnak", 2.724),
        ("km", 1.0),
        ("kmna", 40.0),
        ("knaca", 1000.0),
        ("ksat", 0.1),
        ("alpha", 2.5),
        ("gamma", 0.35),
        ("kmca", 1.38),
        ("kmnai", 87.5),
        ("kpca", 0.0005),
        ("k1", 0.15),
        ("k2", 0.045),
        ("k3", 0.06),
        ("k4", 0.005),
        ("ec", 1.5),
        ("maxsr", 2.5),
        ("minsr", 1.0),
        ("vrel", 0.102),
        ("vleak", 0.00036),
        ("vxfer", 0.0038),
        ("vmaxup", 0.006375),
        ("kup", 0.00025),
        ("bufc", 0.2),
        ("kbufc", 0.001),
        ("bufsr", 10.0),
        ("kbufsf", 0.3),
        ("bufss", 0.4),
        ("kbufss", 0.00025),
        ("gks", 0.392),
        ("gto", 0.073),
        ("v", -86.709),
        ("ki", 138.4),
        ("nai", 10.355),
        ("cai", 0.00013),
        ("cass", 0.00036),
        ("casr", 3.715),
        ("rpri", 0.9068),
        ("xr1", 0.00448),
        ("xr2", 0.476),
        ("xs", 0.0087),
        ("m", 0.00155),
        ("h", 0.7573),
        ("j", 0.7225),
        ("d", 3.164e-05),
        ("f", 0.8009),
        ("f2", 0.9778),
        ("fcass", 0.9953),
        ("s", 0.3212),
        ("r", 2.235e-08),
    ]
)

epicelldata = dict(
    [
        ("gas_constant", 8314.472),
        ("t", 310),
        ("faraday_constant", 96485.3415),
        ("cm", 0.185),
        ("vc", 0.016404),
        ("vsr", 0.001094),
        ("vss", 5.468e-05),
        ("pkna", 0.03),
        ("ko", 5.4),
        ("nao", 140.0),
        ("cao", 2.0),
        ("gk1", 5.405),
        ("gkr", 0.153),
        ("gna", 14.838),
        ("gbna", 0.0002),
        ("gcal", 3.98e-05),
        ("gbca", 0.000592),
        ("gpca", 0.1238),
        ("gpk", 0.0146),
        ("pnak", 2.724),
        ("km", 1.0),
        ("kmna", 40.0),
        ("knaca", 1000.0),
        ("ksat", 0.1),
        ("alpha", 2.5),
        ("gamma", 0.35),
        ("kmca", 1.38),
        ("kmnai", 87.5),
        ("kpca", 0.0005),
        ("k1", 0.15),
        ("k2", 0.045),
        ("k3", 0.06),
        ("k4", 0.005),
        ("ec", 1.5),
        ("maxsr", 2.5),
        ("minsr", 1.0),
        ("vrel", 0.102),
        ("vleak", 0.00036),
        ("vxfer", 0.0038),
        ("vmaxup", 0.006375),
        ("kup", 0.00025),
        ("bufc", 0.2),
        ("kbufc", 0.001),
        ("bufsr", 10.0),
        ("kbufsf", 0.3),
        ("bufss", 0.4),
        ("kbufss", 0.00025),
        ("gks", 0.392),
        ("gto", 0.294),
        ("v", -85.23),
        ("ki", 136.89),
        ("nai", 8.604),
        ("cai", 0.000126),
        ("cass", 0.00036),
        ("casr", 3.64),
        ("rpri", 0.9073),
        ("xr1", 0.00621),
        ("xr2", 0.4712),
        ("xs", 0.0095),
        ("m", 0.00172),
        ("h", 0.7444),
        ("j", 0.7045),
        ("d", 3.373e-05),
        ("f", 0.7888),
        ("f2", 0.9755),
        ("fcass", 0.9953),
        ("s", 0.999998),
        ("r", 2.42e-08),
    ]
)


def test_cellmodel():
    tentusendo = CellModel.TentusscherEndo()
    tentusmid = CellModel.TentusscherMid()
    tentusepi = CellModel.TentusscherEpi()

    assert tentusendo.to_dictionary().items() == endocelldata.items()
    assert tentusmid.to_dictionary().items() == midcelldata.items()
    assert tentusepi.to_dictionary().items() == epicelldata.items()


def test_active():
    tentusepi = CellModel.TentusscherEpi()
    tentusendo = CellModel.TentusscherEndo()
    active0 = EPMaterial.Active(sigma_fiber=1)
    assert active0.sigma_sheet is not None
    active = EPMaterial.Active(sigma_fiber=1, sigma_sheet=1, sigma_sheet_normal=1)
    assert active.sigma_sheet_normal == 1
    assert active.cell_model.to_dictionary().items() == tentusepi.to_dictionary().items()
    active.cell_model = CellModel.TentusscherEndo()
    assert active.cell_model.to_dictionary().items() == tentusendo.to_dictionary().items()
    active_beam = EPMaterial.ActiveBeam(sigma_fiber=1)
    assert active_beam.pmjres is not None


def test_passive():
    passive = EPMaterial.Passive(sigma_fiber=1, sigma_sheet_normal=4)
    assert not hasattr(passive, "cell_model")
    assert passive.sigma_sheet_normal == 4


def test_insulator():
    insulator = EPMaterial.Insulator()
    assert insulator.sigma_fiber == 0
    assert insulator.beta == 0
    assert insulator.cm == 0
