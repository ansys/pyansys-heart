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

"""Module for active stress curve."""

from typing import Literal

import matplotlib.pyplot as plt
import numpy as np

from ansys.health.heart import LOG as LOGGER


def strocchi_active(t_end=800, t_act=0) -> tuple[np.ndarray, np.ndarray]:
    """
    Active stress in doi.org/10.1371/journal.pone.0235145.

    T_peak is described in MAT_295

    Parameters
    ----------
    t_end : int, optional
        heart beat period, by default 800
    t_act : int, optional
        start time, by default 0

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        (time, stress) array
    """
    # parameters used in Strocchi in ms
    tau_r = 130 * 800 / t_end
    tau_d = 100 * 800 / t_end
    tau_dur = 550 * 800 / t_end
    tau_emd = 0.0  # EM coupling delay
    # t_act = 0.0  # activation time from Eikonel model

    def _stress():
        # Active tension
        t = np.linspace(0, t_end, 101)
        active_stress = np.zeros(t.shape)
        ts = t - t_act - tau_emd
        for i, tt in enumerate(ts):
            if 0 < tt < tau_dur:
                active_stress[i] = np.tanh(tt / tau_r) ** 2 * np.tanh((tau_dur - tt) / tau_d) ** 2
        return (t, active_stress)

    return _stress()


def kumaraswamy_active(t_end=1000) -> tuple[np.ndarray, np.ndarray]:
    """
    Active stress in  Gaëtan Desrues doi.org/10.1007/978-3-030-78710-3_43.

    T_peak is described in MAT295

    Parameters
    ----------
    t_end : int, optional
        heart beat duration, by default 1000

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        (timen,stress) array
    """
    apd90 = 250 * t_end / 1000  # action potential duration
    time_repolarization = 750 * t_end / 1000  #  repolarization time

    time = np.linspace(0, t_end, 101)
    stress = np.zeros(time.shape)

    def _kumaraswamy(a, b, x):
        return 1 - (1 - x**a) ** b

    for i, t in enumerate(time):
        if t < apd90:
            stress[i] = _kumaraswamy(2, 1.5, t / apd90)
        elif t < time_repolarization:
            stress[i] = -_kumaraswamy(2, 3, (t - apd90) / (time_repolarization - apd90)) + 1
    return (time, stress)


def constant_ca2(tb: float = 800, ca2ionm: float = 4.35) -> tuple[np.ndarray, np.ndarray]:
    """Constant ca2 curve for Active model 1.

    Parameters
    ----------
    tb : float, optional
        heart beat period, by default 800
    ca2ionm : : float, optional
        amplitude which equals ca2ionm in MAT_295

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        (time, stress) array
    """
    t = np.linspace(0, tb, 101)
    v = np.ones((101)) * ca2ionm
    # set to 0 so across threshold at start of beat
    v[0:1] = 0
    # set to below threshold at the last 95% of period
    # eg. active stress will be disabled at 760ms with tb of 800ms
    v[-5:] = 0
    return (t, v)


class ActiveCurve:
    """Active stress or Ca2+ curve."""

    def __init__(
        self,
        func: tuple[np.ndarray, np.ndarray],
        type: Literal["stress", "ca2"] = "ca2",
        threshold: float = 0.5e-6,
        n: int = 5,
    ) -> None:
        """Define a curve for active behavior of MAT295.

        Parameters
        ----------
        func : tuple[np.ndarray, np.ndarray]
            (time, stress or ca2) array for one heart beat
        type : Literal[&quot;stress&quot;, &quot;ca2&quot;], optional
            type of curve, by default "ca2"
        threshold : float, optional
            threshold of des/active active stress, by default 0.5e-6.
        n : int, optional
            No. of heart beat will be written for LS-DYNA, by default 5

        Notes
        -----
        - If type=='stress', threshold is always 0.5e-6 and ca2+ will be shifted up with 1.0e-6
        except t=0. This ensures a continuous activation during simulation.
        """
        self.type = type
        self.n_beat = n

        if type == "stress":
            LOGGER.warning("Threshold will be reset.")
            threshold = 0.5e-6
        self.threshold = threshold

        self.time = func[0]
        self.t_beat = self.time[-1]

        if self.type == "ca2":
            self.ca2 = func[1]
            self.stress = None
        elif self.type == "stress":
            self.stress = func[1]
            self.ca2 = self._stress_to_ca2(func[1])

        self._check_threshold()

    def _check_threshold(self):
        # maybe better to check it cross 1 or 2 times
        if np.max(self.ca2) < self.threshold or np.min(self.ca2) > self.threshold:
            raise ValueError("Threshold must cross ca2+ curve at least once")

    @property
    def dyna_input(self):
        """Return x,y input for k files."""
        return self._repeat((self.time, self.ca2))

    def plot_time_vs_ca2(self):
        """Plot Ca2+ with threshold."""
        fig, ax = plt.subplots(figsize=(8, 4))
        t, v = self._repeat((self.time, self.ca2))
        ax.plot(t, v, label="Ca2+")
        ax.hlines(self.threshold, xmin=t[0], xmax=t[-1], label="threshold", colors="red")
        ax.set_xlabel("time (ms)")
        ax.set_ylabel("Ca2+")
        # ax.set_title('Ca2+')
        ax.legend()
        return fig

    def plot_time_vs_stress(self):
        """Plot stress."""
        if self.stress is None:
            LOGGER.error("Only support stress curve.")
            # self._estimate_stress()
            return None
        t, v = self._repeat((self.time, self.stress))
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(t, v)
        ax.set_xlabel("time (ms)")
        ax.set_ylabel("Normalized active stress")
        # ax.set_title('Ca2+')
        # ax.legend()
        return fig

    def _stress_to_ca2(self, stress):
        if np.min(stress) < 0 or np.max(stress) > 1.0:
            LOGGER.error("Stress curve is not between 0-1.")
            raise ValueError("Stress curve must be between 0-1.")

        # assuming actype=3, eta=0; n=1; Ca2+50=1
        ca2 = 1 / (1 - 0.999 * stress) - 1

        # offset about threshold
        ca2[0] = 0.0
        ca2[1:] += 2 * self.threshold

        return ca2

    def _repeat(self, curve):
        t = np.copy(curve[0])
        v = np.copy(curve[1])

        for ii in range(1, self.n_beat):
            t = np.append(t, curve[0][1:] + ii * self.t_beat)
            v = np.append(v, curve[1][1:])
        return (t, v)

    def _estimate_stress(self):
        # TODO: only with 1
        # TODO: @wenfengye ensure ruff compatibility, see the noqa's
        ca2ionmax = 4.35
        ca2ion = 4.35
        n = 2
        mr = 1048.9
        dtmax = 150
        tr = -1429
        # Range of L 1.78-1.91
        L = 1.85  # noqa N806
        l0 = 1.58
        b = 4.75
        lam = 1
        cf = (np.exp(b * (lam * L - l0)) - 1) ** 0.5
        ca2ion50 = ca2ionmax / cf
        dtr = mr * lam * L + tr
        self.stress = np.zeros(self.ca2.shape)
        for i, t in enumerate(self.time):
            if t < dtmax:
                w = np.pi * t / dtmax
            elif dtmax <= t <= dtmax + dtr:
                w = np.pi * (t - dtmax + dtr) / dtr
            else:
                w = 0
            c = 0.5 * (1 - np.cos(w))
            self.stress[i] = c * ca2ion**n / (ca2ion**n + ca2ion50**n)


if __name__ == "__main__":
    a = ActiveCurve(constant_ca2(), threshold=0.1, type="ca2")
    # a = Ca2Curve(unit_constant_ca2(), type="ca2")
    a.plot_time_vs_ca2()
    a.plot_time_vs_stress()
