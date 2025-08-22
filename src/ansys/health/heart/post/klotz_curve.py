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

"""Klotz curve module."""

import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np


class EDPVR:
    """End diastolic pressure-volume relation.

    Notes
    -----
    Ref: Klotz, et al. Nature protocols 2.9 (2007): 2152-2158.
    """

    # human constant
    An = 27.78  # mmHg
    Bn = 2.76  # mmHg

    def __init__(self, vm: float, pm: float):
        """Initialize Klotz curve with end diastolic volume and pressure.

        Parameters
        ----------
        vm : float
            Volume in mL.
        pm : float
            Pressure in mmHg
        """
        self.vm = vm
        self.pm = pm

        self.v0 = self.vm * (0.6 - 0.006 * self.pm)
        self.v30 = self.v0 + (self.vm - self.v0) / (self.pm / self.An) ** (1 / self.Bn)

        if self.pm <= 22:
            self.Beta = np.log10(self.pm / 30) / np.log10(self.vm / self.v30)
            self.Alpha = 30 / self.v30**self.Beta
        else:
            v15 = 0.8 * (self.v30 - self.v0) + self.v0
            self.Beta = np.log10(self.pm / 15) / np.log10(self.vm / v15)
            self.Alpha = self.pm / self.vm**self.Beta

    def _get_constants(self) -> tuple[float, float]:
        """Get constants."""
        return self.Alpha, self.Beta

    def get_pressure(self, volume: float | np.ndarray) -> float | np.ndarray:
        """Compute pressure from volume.

        Parameters
        ----------
        volume : float | np.ndarray
            Volume in mL.

        Returns
        -------
        float| np.ndarray
            Pressure in mmHg.
        """
        return self.Alpha * volume**self.Beta

    def get_volume(self, pressure: np.ndarray) -> np.ndarray:
        """Compute volume from pressure.

        Parameters
        ----------
        pressure : np.ndarray
            Pressure in mmHg.

        Returns
        -------
        np.ndarray
            Volume in mL.
        """
        if not isinstance(pressure, np.ndarray):
            raise TypeError("Input must be one-dimensional np.array.")
        volume = np.zeros(pressure.shape)
        for i, p in enumerate(pressure):
            volume[i] = (p / self.Alpha) ** (1 / self.Beta)
            # handle singular issue in Klotz curve
            if volume[i] <= self.v0:
                volume[i] = self.v0
        return volume

    def plot_EDPVR(self, simulation_data: list = None) -> matplotlib.figure.Figure:  # noqa: N802
        """Plot Llotz curve  with simulation data if it exists.

        Parameters
        ----------
        simulation_data : list, default: None
            ``[volume, pressure]`` from simulation.

        Returns
        -------
        matplotlib.figure.Figure
            Figure.
        """
        vv = np.linspace(0, 1.1 * self.vm, num=101)
        pp = self.get_pressure(vv)

        fig = plt.figure()

        plt.plot(vv, pp, label="Klotz curve")
        plt.plot(self.v0, 0, "o", label="Klotz V0")
        plt.plot(self.vm, self.pm, "o", label="V_ED, P_ED")

        if simulation_data is not None:
            plt.plot(simulation_data[0], simulation_data[1], "--*", label="Simulation")

        plt.title("EDVPR", fontsize=14)
        plt.xlabel("Volume (mL)", fontsize=14)
        plt.ylabel("Pressure (mmHg)", fontsize=14)
        plt.legend()

        return fig
