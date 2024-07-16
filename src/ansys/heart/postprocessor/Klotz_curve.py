# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
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

    def __init__(self, vm, pm):
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

    def get_constants(self):
        """Get constants."""
        return self.Alpha, self.Beta

    def get_pressure(self, volume):
        """Get pressure."""
        return self.Alpha * volume**self.Beta

    def get_volume(self, pressure):
        """Get volume."""
        volume = np.zeros(pressure.shape)
        for i, p in enumerate(pressure):
            volume[i] = (p / self.Alpha) ** (1 / self.Beta)
            # handle singular issue in Klotz curve
            if volume[i] <= self.v0:
                volume[i] = self.v0
        return volume

    def plot_EDPVR(self, simulation_data=None):
        """
        Plot end-diastolic pressure-volume relation.

        Parameters
        ----------
        simulation_data: optional, plot if exist.
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


if __name__ == "__main__":
    # healthy baseline
    v_ed = 1.752e02  # mL
    p_ed = 15  # mmHg

    # v_ed = 120  # mL
    # p_ed = 7  # mmHg
    edpvr = EDPVR(v_ed, p_ed)

    edpvr.plot_EDPVR()
