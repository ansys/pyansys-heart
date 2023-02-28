"""Module for postprocessing system model data."""
from dataclasses import dataclass
import json
import os


from ansys.heart.simulator.settings.settings import SimulationSettings
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


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

    Notes: future use.
    """

    pressure: Pressure
    flow: Flow
    volume: Volume


class ZeroDSystem:
    """0D circulation system model (for one cavity)."""

    def __init__(self, csv_path, ed_state, name=""):
        """
        Init.

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

        Note
        ----
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

    Note
    ----
    unit: ms, kPa, mL
    """

    def __init__(self, dir):
        """
        Init.

        Parameters
        ----------
        dir: simulation directory
        """
        self.dir = dir

        # get EOD pressure
        s = SimulationSettings()
        s.load(os.path.join(self.dir, "simulation_settings.yml"))
        lp = s.mechanics.boundary_conditions.end_diastolic_cavity_pressure.left_ventricle.to(
            "kilopascal"
        ).m
        rp = s.mechanics.boundary_conditions.end_diastolic_cavity_pressure.right_ventricle.to(
            "kilopascal"
        ).m

        # get EOD volume
        # todo: get this information from binout:icvout
        try:
            # load simulated EOD volume
            with open(os.path.join(self.dir, "Post_report.json")) as f:
                dct = json.load(f)
                lv = dct["Simulation Left ventricle volume (mm3)"][-1] / 1000
                rv = dct["Simulation Right ventricle volume"][-1] / 1000
        except FileExistsError:
            # load Input EOD volume
            with open(os.path.join(self.dir, "volumes.json")) as f:
                dct = json.load(f)
                lv = dct["Left ventricle"] / 1000
                rv = dct["Right ventricle"] / 1000

        f = os.path.join(self.dir, "constant_preload_windkessel_afterload_left.csv")
        self.lv = ZeroDSystem(f, [lp, lv], name="Left ventricle")

        f = os.path.join(self.dir, "constant_preload_windkessel_afterload_right.csv")
        self.rv = ZeroDSystem(f, [rp, rv], name="Right ventricle")

    def plot_pv_loop(self, t_start=0, t_end=10e10):
        """
        Plot PV loop.

        Parameters
        ----------
        t_start: start time
        t_end: end time
        """
        start = np.where(self.lv.time >= t_start)[0][0]
        end = np.where(self.lv.time <= t_end)[0][-1]

        fig, axis = plt.subplots()
        fig.suptitle("Pressure Volume Loop")

        def add_pv(cavity, color):
            v = cavity.volume.cavity[start:end]
            ef = (max(v) - min(v)) / max(v)
            p = cavity.pressure.cavity[start:end]
            axis.plot(v, p, color, label="{0},EF={1:.1f}%".format(cavity.name, ef * 100))
            axis.scatter(cavity.ed[1], cavity.ed[0], facecolor=color, label=cavity.name + "@ED")
            return

        add_pv(self.lv, "blue")
        try:
            add_pv(self.rv, "red")
        except:
            pass

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
        try:
            from ansys.heart.postprocessor.binout_helper import IcvOut
            self.bin = IcvOut(os.path.join(self.dir, "binout"))
        except FileExistsError:
            from ansys.heart.postprocessor.binout_helper import IcvOut
            self.bin = IcvOut(os.path.join(self.dir, "binout0000"))

        self.bin.time /= 1000  # ms ->s
        self.bin.pressure *= 1000  # MPa -> kPa
        self.bin.volume /= 1000  # mm^3 -> mL
        self.bin.flow /= 1  # mm^3/ms ->mL/s

        fig, axis = plt.subplots(3, figsize=(8, 4))

        axis[0].plot(self.bin.time, self.bin.pressure[:, 0], label="pressure_binout")
        axis[0].plot(self.lv.time, self.lv.pressure.cavity, "--", label="pressure_csv")
        axis[0].legend()

        axis[1].plot(self.bin.time, self.bin.flow[:, 0], label="flow_binout")
        axis[1].plot(self.lv.time, -self.lv.flow.cavity, "--", label="flow_csv")
        axis[1].legend()

        axis[2].plot(self.bin.time, self.bin.volume[:, 0], label="volume_binout")
        axis[2].plot(self.lv.time, self.lv.volume.cavity, "--", label="volume_csv")
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
            print("future development.")
            return

        fig, axis = plt.subplots(figsize=(8, 4))

        if self.type == "LV":
            vlv = interp1d(self.bin.time, self.bin.volume[:, 0])(self.lv["time"])
            v_total = vlv + self.lv["vart"] + self.lv["vven"]
            if plot_all:
                axis.plot(self.lv["time"], self.lv["vart"], label="vart")
                axis.plot(self.lv["time"], self.lv["vven"], label="vven")
                axis.plot(self.lv["time"], vlv, label="vlv")

        elif self.type == "BV":
            vlv = interp1d(self.bin.time, self.bin.volume[:, 0])(self.lv["time"])
            vrv = interp1d(self.bin.time, self.bin.volume[:, 1])(self.rv["time"])
            v_total = (
                vlv
                + self.lv["vart"]
                + self.lv["vven"]
                + vrv
                + self.lv["v_part"]
                + self.lv["v_pven"]
            )
            if plot_all:
                axis.plot(self.lv["time"], self.lv["vart"], label="vart")
                axis.plot(self.lv["time"], self.lv["v_part"], label="v_part")
                axis.plot(self.lv["time"], self.lv["vven"], label="vven")
                axis.plot(self.lv["time"], self.lv["v_pven"], label="v_pven")
                axis.plot(self.lv["time"], vlv, label="vlv")
                axis.plot(self.lv["time"], vrv, label="vrv")
        variation = abs(max(v_total) - min(v_total)) / max(v_total)

        axis.plot(
            self.lv["time"][1:],
            v_total[1:],
            label="V_total, var={0:.1f}%".format(variation * 100),
        )
        axis.set_ylabel("Volume (mL)")
        axis.set_xlabel("Time (s)")
        axis.legend()

        return fig
