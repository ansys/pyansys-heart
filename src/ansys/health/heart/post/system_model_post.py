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

"""Module for postprocessing system model data."""

from dataclasses import dataclass
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.post.dpf_utils import ICVoutReader
from ansys.health.heart.settings.settings import SimulationSettings


@dataclass(init=False)
class Pressure:
    """System state for pressure."""

    cavity: np.ndarray
    artery: np.ndarray
    venous: np.ndarray


@dataclass(init=False)
class Flow:
    """System state for flow."""

    cavity: np.ndarray
    artery: np.ndarray
    venous: np.ndarray
    peripheral: np.ndarray


@dataclass(init=False)
class Volume:
    """System state for volume."""

    cavity: np.ndarray
    artery: np.ndarray
    venous: np.ndarray


@dataclass
class SystemState:
    """
    System state including pressure, flow, and volume.

    Notes
    -----
    Future development.
    """

    pressure: Pressure
    flow: Flow
    volume: Volume


class ZeroDSystem:
    """0D circulation system model (for one cavity)."""

    def __init__(self, csv_path: str, ed_state: list[float, float], name: str = ""):
        """Initialize ZeroDSystem.

        Parameters
        ----------
        csv_path : str
            Path to the CSV file.
        ed_state : list[float,float]
            End of diastole pressure and volume.
        name : str, default: ""
            Cavity name.
        """
        self.name = name
        self.ed = ed_state

        data = pd.read_csv(csv_path)

        self.time = data["time"].to_numpy() / 1000

        self.pressure = Pressure()
        self.pressure.cavity = data["pk"].to_numpy() * 1000
        self.pressure.artery = data["part"].to_numpy() * 1000
        self.pressure.venous = data["pven"].to_numpy() * 1000

        self.flow = Flow()
        self.flow.cavity = data["qk"].to_numpy()
        self.flow.artery = data["qart"].to_numpy()
        self.flow.venous = data["qven"].to_numpy()
        self.flow.peripheral = data["qp"].to_numpy()

        self.volume = Volume()
        self.volume.artery = data["vart"].to_numpy() / 1000
        # integrate volume of cavity
        self.volume.cavity = self._integrate_volume(self.ed[1], self.time, self.flow.cavity)

        pass

    @staticmethod
    def _integrate_volume(v0: float, t: np.ndarray, q: np.ndarray) -> np.ndarray:
        """Integrate the cavity's volume.

        Notes
        -----
        Cavity's volume is not evaluated/saved in the CSV file. This is to ensure
        that volume is consistent with what's in the ICVOUT file.

        This assumes that the implicit solver with ``gamma=0.6`` was used.

        Parameters
        ----------
        v0 : float
            Volume at t0
        t : np.ndarray
            Time array
        q : np.ndarray
            Flow array

        Returns
        -------
        np.ndarray
            Cavity volume.
        """
        gamma = 0.6

        v = np.zeros(len(q))
        v[0] = v0
        for i in range(1, len(t)):
            v[i] = v[i - 1] + (t[i] - t[i - 1]) * ((1 - gamma) * q[i - 1] + gamma * q[i])

        return v


class SystemModelPost:
    """
    Postprocessing system model.

    Notes
    -----
    The units are ms, kPa, and mL.
    """

    def __init__(self, dir: str):
        """Initialize ``SystemModelPost``.

        Parameters
        ----------
        dir : str
            Simulation directory.
        """
        self.dir = dir
        self.model_type = "LV"

        fcsv1 = os.path.join(self.dir, "constant_preload_windkessel_afterload_left.csv")
        fcsv2 = os.path.join(self.dir, "constant_preload_windkessel_afterload_right.csv")
        if os.path.isfile(fcsv2):
            self.model_type = "BV"

        # get EOD pressure
        s = SimulationSettings()
        s.load(os.path.join(self.dir, "simulation_settings.yml"))
        l_ed_pressure = (
            s.mechanics.boundary_conditions.end_diastolic_cavity_pressure.left_ventricle.to(
                "kilopascal"
            ).m
        )
        if self.model_type == "BV":
            r_ed_pressure = (
                s.mechanics.boundary_conditions.end_diastolic_cavity_pressure.right_ventricle.to(
                    "kilopascal"
                ).m
            )

        # get EOD volume
        try:
            icvout = ICVoutReader(os.path.join(self.dir, "binout0000"))
        except FileNotFoundError:
            try:  # from SMP
                icvout = ICVoutReader(os.path.join(self.dir, "binout"))
            except FileNotFoundError as error:
                LOGGER.error(f"Cannot find binout file. {error}")
                raise FileNotFoundError(f"Cannot find binout file. {error}")

        l_ed_volume = icvout.get_volume(1)[0] / 1000
        self.lv_system = ZeroDSystem(fcsv1, [l_ed_pressure, l_ed_volume], name="Left ventricle")

        if self.model_type == "BV":
            r_ed_volume = icvout.get_volume(2)[0] / 1000
            self.rv_system = ZeroDSystem(
                fcsv2, [r_ed_pressure, r_ed_volume], name="Right ventricle"
            )

    def get_ejection_fraction(self, t_start: float = 0, t_end: float = 10e10) -> float:
        """Compute ejection fraction on a given time interval.

        Parameters
        ----------
        t_start : float, default: 0
            Start time.
        t_end : float, default: 10e10
            End time.

        Returns
        -------
        float
            Ejection fraction.
        """
        ef = [None, None]
        start = np.where(self.lv_system.time >= t_start)[0][0]
        end = np.where(self.lv_system.time <= t_end)[0][-1]
        vl = self.lv_system.volume.cavity[start:end]
        try:
            ef[0] = (max(vl) - min(vl)) / max(vl)
        except Exception as e:
            ef[0] = None
            LOGGER.warning(f"Failed to compute ejection fraction. {e}")
        if self.model_type == "BV":
            vr = self.rv_system.volume.cavity[start:end]
            ef[1] = (max(vr) - min(vr)) / max(vr)

        return ef

    def plot_pv_loop(
        self,
        t_start: float = 0,
        t_end: float = 10e10,
        show_ed: bool = True,
        ef: list[float, float] = [None, None],
    ) -> plt.Figure:
        """Plot PV loop.

        Parameters
        ----------
        t_start : float, default: 0
            Start time to plot.
        t_end : float, default: 10e10
            End time to plot.
        show_ed : bool, default: True
            Whether to show the end of the diastole state in zeropressure.
        ef : list[float, float], default: [None, None]
            Ejection fraction to show in the legend.

        Returns
        -------
        plt.Figure
            Figure handle.
        """
        start = np.where(self.lv_system.time >= t_start)[0][0]
        end = np.where(self.lv_system.time <= t_end)[0][-1]

        fig, axis = plt.subplots()
        fig.suptitle("Pressure Volume Loop")
        if self.model_type == "LV":
            axis.set_xlim(
                [
                    0.95 * np.min(self.lv_system.volume.cavity),
                    1.05 * np.max(self.lv_system.volume.cavity),
                ]
            )
            axis.set_ylim(
                [
                    0.8 * np.min(self.lv_system.pressure.cavity),
                    1.2 * np.max(self.lv_system.pressure.cavity),
                ]
            )
        else:
            axis.set_xlim(
                [
                    0.95 * np.min(self.lv_system.volume.cavity),
                    1.05 * np.max(self.rv_system.volume.cavity),
                ]
            )
            axis.set_ylim(
                [
                    0.8 * np.min(self.rv_system.pressure.cavity),
                    1.2 * np.max(self.lv_system.pressure.cavity),
                ]
            )

        def add_pv(cavity, color, ef=None):
            v = cavity.volume.cavity[start:end]
            p = cavity.pressure.cavity[start:end]

            # label
            label = "{0}".format(cavity.name)
            if ef is not None:
                label = "{0},EF={1:.1f}%".format(label, ef * 100)

            # plot
            axis.plot(v, p, color, label=label)

            if show_ed:
                axis.scatter(cavity.ed[1], cavity.ed[0], facecolor=color, label=cavity.name + "@ED")
            else:  # highlight last point
                if len(v) > 0:  # safety
                    axis.scatter(v[-1], p[-1], facecolor=color)
            return

        add_pv(self.lv_system, "blue", ef=ef[0])
        if self.model_type == "BV":
            add_pv(self.rv_system, "red", ef=ef[1])

        axis.set_xlabel("Volume (mL)")
        axis.set_ylabel("Pressure (kPa)")

        ax2 = axis.twinx()
        mn, mx = axis.get_ylim()
        ax2.set_ylim(mn * 7.50062, mx * 7.50062)  # kPa --> mmHg
        ax2.set_ylabel("(mmHg)")

        axis.legend()

        return fig
