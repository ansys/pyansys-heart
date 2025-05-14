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

"""
Use PyDYNA keywords module to create commonly used material cards and their default values.

Notes
-----
Examples of material cards include Mat295, Mat077, MatNull.

"""

import dataclasses

# from importlib.resources import files
from importlib.resources import path as resource_path

import numpy as np
import pandas as pd

from ansys.dyna.core.keywords import keywords
from ansys.health.heart.settings.material.material import ActiveModel, Mat295

# import custom keywords in separate namespace
from ansys.health.heart.writer import custom_keywords as custom_keywords


class MaterialCap(keywords.MatNull):
    """Material of the closing cap/valves.

    Parameters
    ----------
    keywords : keywords.MatNull
        Inherits from the Null type material.
    """

    def __init__(self, mid: int = 1) -> None:
        super().__init__(mid=mid, ro=1.04e-6)


class MaterialNeoHook(custom_keywords.Mat077H):
    """Material for the atrium.

    Parameters
    ----------
    mid : int
        Material ID.
    rho : float
        Density of the material.
    c10 : float
        First coefficient of the material.
    nu : float
        Poisson's ratio.
    kappa : float
        Bulk modulus.
    """

    def __init__(
        self,
        mid: int,
        rho: float,
        c10: float,
        nu: float,
        kappa: float,
    ) -> None:
        super().__init__(mid=mid, ro=rho, pr=nu, n=0, c10=c10)
        setattr(self, "user_comment", f"nu deduced from kappa={kappa}")
        return


class MaterialHGOMyocardium(keywords.Mat295):
    """HGO material model, which is derived from Mat295."""

    def __init__(self, id: int, mat: Mat295, ignore_active: bool = False) -> None:
        """Init a keyword of *mat295.

        Parameters
        ----------
        id : int
            Material ID.
        mat : Mat295
            Material data.
        ignore_active : bool, default: False
            Whether to ignore the active module. For example, for stress-free.
        """
        # 1st line
        super().__init__(mid=id)
        setattr(self, "rho", mat.rho)
        setattr(self, "aopt", mat.aopt)
        setattr(self, "user_comment", f"nu deduced from kappa={mat.iso.kappa}")
        # iso
        for field in dataclasses.fields(mat.iso):
            value = getattr(mat.iso, field.name)
            setattr(self, field.name, value)

        # aniso
        if mat.aniso is not None:
            self.atype = mat.aniso.atype
            self.intype = mat.aniso.intype
            self.nf = mat.aniso.nf
            self.ftype = mat.aniso.fibers[0]._ftype  # not used but must be defined

            self.a1 = mat.aniso.vec_a[0]
            self.a2 = mat.aniso.vec_a[1]
            self.a3 = mat.aniso.vec_a[2]

            self.d1 = mat.aniso.vec_d[0]
            self.d2 = mat.aniso.vec_d[1]
            self.d3 = mat.aniso.vec_d[2]

            fiber_sheet = []
            for i in range(len(mat.aniso.fibers)):
                dct = {
                    "theta": mat.aniso.fibers[i]._theta,
                    "a": mat.aniso.fibers[i].a,
                    "b": mat.aniso.fibers[i].b,
                    "ftype": mat.aniso.fibers[i]._ftype,
                    "fcid": mat.aniso.fibers[i]._fcid,
                    "k1": mat.aniso.fibers[i].k1,
                    "k2": mat.aniso.fibers[i].k2,
                }
                fiber_sheet.append(dct)
            self.anisotropic_settings = pd.DataFrame(fiber_sheet)

            if mat.aniso.intype == 1:
                self.coupling_k1 = mat.aniso.k1fs
                self.coupling_k2 = mat.aniso.k2fs

        # active
        if not ignore_active and mat.active is not None:
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


def active_curve(
    curve_type: str = "Strocchi2020",
    endtime: float = 15,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute various (normalized) curves used for the active module.

    Parameters
    ----------
    curve_name : str
        Type of curve to compute.
    """
    # time array
    # T = np.arange( 0, endtime, timestep )
    # NOTE: needs cleaning up
    if curve_type == "Strocchi2020":
        # parameters used in Strocchi in ms
        t_end = 800
        tau_r = 130
        tau_d = 100
        tau_dur = 550
        tau_emd = 0.0  # EM coupling delay
        t_act = 0.0  # activation time from Eikonel model

        # Active tension
        t = np.linspace(0, t_end, 1001)
        active_stress = np.zeros(t.shape)
        ts = t - t_act - tau_emd
        for i, tt in enumerate(ts):
            if 0 < tt < tau_dur:
                active_stress[i] = np.tanh(tt / tau_r) ** 2 * np.tanh((tau_dur - tt) / tau_d) ** 2

        # repeat dataset nCycles times:
        # number of cycles to return
        num_cycles = int(np.ceil(endtime / t_end))

        time_array = t  # time array
        # mock calcium array
        calcium_array0 = 1 / (1 - 0.99 * active_stress) - 1
        calcium_array = np.copy(calcium_array0)
        for ii in range(1, num_cycles):
            time_array = np.append(time_array, t[1:] + ii * t_end)
            calcium_array = np.append(calcium_array, calcium_array0[1:])

    # used for generating multi beats with model actype 1
    elif curve_type == "constant":
        nb_beats = 10
        period = 1.0  # in second

        # define shape pattern
        value = np.array([0, 1, 1])
        time = np.array([0, 0.001 * period, 0.9 * period])

        # repeat for every period
        time_array = time
        for i in range(1, nb_beats):
            time_array = np.append(time_array, time + period * i)
        calcium_array = np.tile(value, nb_beats)

        # append last point
        time_array = np.append(time_array, period * nb_beats)
        calcium_array = np.append(calcium_array, 0.0)

    elif curve_type == "TrueCalcium":
        file_path = resource_path("ansys.health.heart.writer", "calcium_from_EP.txt").__enter__()
        a = np.loadtxt(file_path)
        time_array = a[:, 0] / 1000
        calcium_array = a[:, 1]

    # import matplotlib.pyplot as plt
    # plt.plot(time_array, calcium_array)
    # plt.show()

    return time_array, calcium_array
