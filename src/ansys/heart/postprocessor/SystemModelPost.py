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

"""Module for postprocessing system model data."""

from dataclasses import dataclass
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from ansys.heart.core import LOG as LOGGER
from ansys.heart.postprocessor.dpf_utils import ICVoutReader
from ansys.heart.simulator.settings.settings import SimulationSettings


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
    System state including pressure, flow, volume.

    Notes
    -----
    future development
    """

    pressure: Pressure
    flow: Flow
    volume: Volume


class ZeroDSystem:
    """0D circulation system model (for one cavity)."""

    def __init__(self, csv_path, ed_state, name=""):
        """
        Initialize ZeroDSystem.

        Notes
        -----
        # from "ms, MPa, mm^3" to "s, kPa, mL"

        Parameters
        ----------
        csv_path: str, csv file path
        ed_state: List, End of Diastole pressure and volume
        name: str, system name, like 'left_ventricle'

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
        self.volume.cavity = self.integrate_volume(self.ed[1], self.time, self.flow.cavity)

        pass

    @staticmethod
    def integrate_volume(v0, t, q):
        """
        Integrate cavity's volume.

        Notes
        -----
        Cavity's volume is not evaluated/saved in csv file.
        Use implicit with gamma=0.6

        Parameters
        ----------
        v0: float, volume at time of 0
        t: time array
        q: flow array

        Returns
        -------
            volume array
        """
        gamma = 0.6

        v = np.zeros(len(q))
        v[0] = v0
        for i in range(1, len(t)):
            v[i] = v[i - 1] + (t[i] - t[i - 1]) * ((1 - gamma) * q[i - 1] + gamma * q[i])

        return v


class SystemModelPost:
    """
    Class for post-processing system model.

    Notes
    -----
    unit: ms, kPa, mL
    """

    def __init__(self, dir):
        """
        Initialize SystemModelPost.

        Parameters
        ----------
        dir: simulation directory
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
            except FileNotFoundError:
                LOGGER.error("Cannot find binout file.")
                exit()
        l_ed_volume = icvout.get_volume(1)[0] / 1000
        self.lv_system = ZeroDSystem(fcsv1, [l_ed_pressure, l_ed_volume], name="Left ventricle")

        if self.model_type == "BV":
            r_ed_volume = icvout.get_volume(2)[0] / 1000
            self.rv_system = ZeroDSystem(
                fcsv2, [r_ed_pressure, r_ed_volume], name="Right ventricle"
            )

    def get_ejection_fraction(self, t_start=0, t_end=10e10):
        """
        Compute ejection fraction at a given time interval.

        Parameters
        ----------
        t_start: start time
        t_end: end time

        Returns
        -------
        Ejection fraction
        """
        ef = [None, None]
        start = np.where(self.lv_system.time >= t_start)[0][0]
        end = np.where(self.lv_system.time <= t_end)[0][-1]
        vl = self.lv_system.volume.cavity[start:end]
        try:
            ef[0] = (max(vl) - min(vl)) / max(vl)
        except:
            ef[0] = None
            LOGGER.warning("Failed to compute ejection fraction.")
        if self.model_type == "BV":
            vr = self.rv_system.volume.cavity[start:end]
            ef[1] = (max(vr) - min(vr)) / max(vr)

        return ef

    def plot_pv_loop(self, t_start=0, t_end=10e10, show_ed=True, ef=[None, None]):
        """
        Plot PV loop.

        Parameters
        ----------
        ef: Default None, else plot ejection fraction in legend.
        show_ed: Default False, else plot ED state
        t_start: start time
        t_end: end time
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

    @staticmethod
    def plot_pressure_flow_volume(cavity: ZeroDSystem, t_start: float = 0, t_end: float = 10e5):
        """Plot pressure/flow/volume curves.

        Parameters
        ----------
        cavity: ZeroDSystem,
        t_end: start time
        t_start: end time

        """
        fig, axis = plt.subplots(3, figsize=(8, 4), sharex="all")
        fig.suptitle(f"{cavity.name}: Pressure & Flow & Volume")

        # define plot x range
        start = np.where(cavity.time >= t_start)[0][0]
        end = np.where(cavity.time <= t_end)[0][-1]
        axis[0].set_xlim([cavity.time[start], cavity.time[end]])

        # find where both valves are closed: iso-volume
        iso_vol = (cavity.pressure.cavity < cavity.pressure.artery) & (
            cavity.pressure.cavity > cavity.pressure.venous
        )
        # start and end time when this state change
        tt = cavity.time[np.where(iso_vol[:-1] != iso_vol[1:])[0]]
        # plot iso-volume phase in grey zone
        try:
            tt = tt.reshape(-1, 2)
        except ValueError:
            tt = np.append(tt, 10e5)
            tt = tt.reshape(-1, 2)

        for i in range(len(tt)):
            for j in range(3):
                axis[j].axvspan(tt[i, 0], tt[i, 1], facecolor="grey", alpha=0.3)

        # do plot
        axis[0].plot(cavity.time, cavity.pressure.cavity, label="cavity")
        axis[0].plot(cavity.time, cavity.pressure.artery, label="artery")
        axis[0].plot(cavity.time, cavity.pressure.venous, label="venous")
        axis[0].set_ylabel("Pressure (kPa)")
        axis[0].legend()

        axis[1].plot(cavity.time, cavity.flow.cavity, label="cavity")
        axis[1].plot(cavity.time, cavity.flow.artery, label="artery")
        axis[1].plot(cavity.time, cavity.flow.venous, label="venous")
        axis[1].plot(cavity.time, cavity.flow.peripheral, label="peripheral")
        axis[1].set_ylabel("Flow (mL/s)")
        axis[1].legend()

        axis[2].plot(cavity.time, cavity.volume.cavity, label="cavity")
        # axis[2].plot(cavity.time, cavity.volume.artery, label="artery")
        axis[2].set_ylabel("Volume (mL)")
        axis[2].set_xlabel("Time (s)")
        axis[1].legend()

        return fig

    def _check_output(self):
        """
        Check if system states == FEM states.

        Notes
        -----
        Only for debug
        Require qd to read binout
        It's normal that the cavity volume is slight different.
        """
        from ansys.heart.postprocessor.deprecated_binout_helper import IcvOut

        try:
            self.bin = IcvOut(os.path.join(self.dir, "binout"))
        except FileExistsError:
            self.bin = IcvOut(os.path.join(self.dir, "binout0000"))

        self.bin.time /= 1000  # ms ->s
        self.bin.pressure *= 1000  # MPa -> kPa
        self.bin.volume /= 1000  # mm^3 -> mL
        self.bin.flow /= 1  # mm^3/ms ->mL/s

        fig, axis = plt.subplots(3, figsize=(8, 4))

        axis[0].plot(self.bin.time, self.bin.pressure[:, 0], label="pressure_binout")
        axis[0].plot(
            self.lv_system.time, self.lv_system.pressure.cavity, "--", label="pressure_csv"
        )
        axis[0].legend()

        axis[1].plot(self.bin.time, self.bin.flow[:, 0], label="flow_binout")
        axis[1].plot(self.lv_system.time, -self.lv_system.flow.cavity, "--", label="flow_csv")
        axis[1].legend()

        axis[2].plot(self.bin.time, self.bin.volume[:, 0], label="volume_binout")
        axis[2].plot(self.lv_system.time, self.lv_system.volume.cavity, "--", label="volume_csv")
        axis[2].legend()
        axis[2].set_xlabel("Time (s)")

        return fig

    def _check_total_volume(self, plot_all=False):
        """Check if total volume is constant for a closed loop.

        Parameters
        ----------
        plot_all

        """
        # if not self.closed_loop:
        if True:
            LOGGER.error("future development.")
            return

        fig, axis = plt.subplots(figsize=(8, 4))

        if self.model_type == "LV":
            vlv = interp1d(self.bin.time, self.bin.volume[:, 0])(self.lv_system["time"])
            v_total = vlv + self.lv_system["vart"] + self.lv_system["vven"]
            if plot_all:
                axis.plot(self.lv_system["time"], self.lv_system["vart"], label="vart")
                axis.plot(self.lv_system["time"], self.lv_system["vven"], label="vven")
                axis.plot(self.lv_system["time"], vlv, label="vlv")

        elif self.model_type == "BV":
            vlv = interp1d(self.bin.time, self.bin.volume[:, 0])(self.lv_system["time"])
            vrv = interp1d(self.bin.time, self.bin.volume[:, 1])(self.rv_system["time"])
            v_total = (
                vlv
                + self.lv_system["vart"]
                + self.lv_system["vven"]
                + vrv
                + self.lv_system["v_part"]
                + self.lv_system["v_pven"]
            )
            if plot_all:
                axis.plot(self.lv_system["time"], self.lv_system["vart"], label="vart")
                axis.plot(self.lv_system["time"], self.lv_system["v_part"], label="v_part")
                axis.plot(self.lv_system["time"], self.lv_system["vven"], label="vven")
                axis.plot(self.lv_system["time"], self.lv_system["v_pven"], label="v_pven")
                axis.plot(self.lv_system["time"], vlv, label="vlv")
                axis.plot(self.lv_system["time"], vrv, label="vrv")
        variation = abs(max(v_total) - min(v_total)) / max(v_total)

        axis.plot(
            self.lv_system["time"][1:],
            v_total[1:],
            label="V_total, var={0:.1f}%".format(variation * 100),
        )
        axis.set_ylabel("Volume (mL)")
        axis.set_xlabel("Time (s)")
        axis.legend()

        return fig
