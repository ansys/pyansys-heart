"""Methods and classes for postprocessing system model data."""
import os

from ansys.heart.postprocessor.binout_helper import IcvOut
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


class SystemModelPost:
    """
    Class for post-processing system model.

    Note
    ----
    unit: ms, kPa, mL
    """

    def __init__(self, dir, closed_loop=False):
        """
        Init.

        Parameters
        ----------
        dir: simulation directory
        closed_loop: if it's a closed loop, default False
        """
        self.dir = dir
        self.closed_loop = closed_loop

        self.bin = IcvOut(os.path.join(self.dir, "binout0000"))
        self.bin.pressure *= 1000  # kPa
        self.bin.volume /= 1000  # mL
        self.bin.flow /= 1000  # mL/time

        # Get the initial cavity volume
        # Note: some reasonable difference between LSDYNA volume and vtk volume
        self.cavity0 = self.bin.volume[0]

        if self.bin.volume.shape[1] == 1:
            self.type = "LV"
            print("LV system model")
        elif self.bin.volume.shape[1] == 2:
            self.type = "BV"
            print("BV system model")
        else:
            print("Unknown system model")

        self.p_ed = np.array([2, 0.5333])  # kPa
        self.v_ed = self.cavity0

        self.compute_ejection_ratio()
        # todo: check if last loop converge

        self._load_csv()

    def _load_csv(self):
        """Load system states written from define function."""
        try:
            self.lv = pd.read_csv(
                os.path.join(self.dir, "constant_preload_windkessel_afterload_left.csv")
            )
        except IOError:
            print("Only accept default file name constant_preload_windkessel_afterload_*.csv")

        #  special fix for Constant pre/after load case
        #  part/pven are not written in Csv files
        if "part" not in self.lv.columns:
            self.lv["part"] = 8 * np.ones(len(self.lv["time"]))
            self.lv["pven"] = 2 * np.ones(len(self.lv["time"]))

        # convert mm^3 to mL
        for name in self.lv.columns:
            if name[0] == "v" or name[0] == "q":
                self.lv[name] /= 1000
            if name[0] == "p":
                self.lv[name] *= 1000

        if self.type == "BV":
            self.rv = pd.read_csv(
                os.path.join(self.dir, "constant_preload_windkessel_afterload_right.csv")
            )

            if "part" not in self.rv.columns:
                self.rv["part"] = 2 * np.ones(len(self.lv["time"]))
                self.rv["pven"] = 0.5333 * np.ones(len(self.lv["time"]))

            for name in self.rv.columns:
                if name[0] == "v" or name[0] == "q":
                    self.rv[name] /= 1000
                if name[0] == "p":
                    self.rv[name] *= 1000

    def compute_ejection_ratio(self, cycle_duration=1000):
        """Compute ejection ratio of last loop."""
        # get PV of last loop
        _, volume = self.get_PV(self.bin.time[-1] - cycle_duration)

        self.lv_ef = (max(volume[:, 0]) - min(volume[:, 0])) / max(volume[:, 0])
        if self.type == "BV":
            self.rv_ef = (max(volume[:, 1]) - min(volume[:, 1])) / max(volume[:, 1])
        return

    def get_PV(self, t_start=0, t_end=10e10):
        """Get Pressure & volume.

        Parameters
        ----------
        t_end
        t_start

        """
        i_start = np.where(self.bin.time >= t_start)[0][0]
        i_end = np.where(self.bin.time <= t_end)[0][-1]

        volume = self.bin.volume[i_start:i_end, :]
        pressure = self.bin.pressure[i_start:i_end, :]
        return pressure, volume

    def plot_pv_loop(self, t_start=0, t_end=10e10):
        """
        Plot PV loop.

        Parameters
        ----------
        t_start: start time
        t_end: end time
        """
        pressure, volume = self.get_PV(t_start, t_end)

        fig, axis = plt.subplots()
        fig.suptitle("Pressure Volume Loop")
        axis.plot(
            volume[:, 0],
            pressure[:, 0],
            "b",
            label="LV,EF={0:.1f}%".format(self.lv_ef * 100),
        )
        axis.scatter(self.v_ed[0], self.p_ed[0], facecolor="blue", label="LV@ED")
        if self.type == "BV":
            axis.plot(
                volume[:, 1],
                pressure[:, 1],
                "r",
                label="RV, EF={0:.1f}%".format(self.rv_ef * 100),
            )
            axis.scatter(self.v_ed[1], self.p_ed[1], facecolor="red", label="RV@ED")
        axis.set_xlabel("Volume (mL)")
        axis.set_ylabel("Pressure (kPa)")
        axis.legend()

        return fig

    def plot_pressure_flow_volume(self, cavity, t_start=0, t_end=10e10):
        """Plot curves.

        Parameters
        ----------
        cavity: 'lv' or 'rv'
        t_end: start time
        t_start: end time

        """
        fig, axis = plt.subplots(3, figsize=(8, 4), sharex=True)

        if cavity == "lv":
            system = self.lv
            id = 0
            fig.suptitle("Left ventricle: Pressure & Flow & Volume")

        elif cavity == "rv":
            system = self.rv
            id = 1
            fig.suptitle("Right ventricle: Pressure & Flow & Volume")

        p_names = ["pk", "pven", "part"]
        q_names = ["qk", "qven", "qart"]
        if "qp" in system.columns:  # add peripheral flow
            q_names.append("qp")

        if self.type == "BV" and self.closed_loop:  # special case for BV closed loop
            if cavity == "lv":
                p_names = ["pvv", "p_pven", "part"]
                q_names = ["qvv", "q_pven", "qart", "qp"]
            elif cavity == "rv":
                p_names = ["pvv", "pven", "p_part"]
                q_names = ["qvv", "qven", "q_part", "q_pp"]

        time = system["time"].values
        # prepare pressure and flow
        pressure = []
        flow = []
        for p_name in p_names:
            pressure.append(system[p_name].values)
        for q_name in q_names:
            flow.append(system[q_name].values)

        # define plot x range
        i_start = np.where(time >= t_start)[0][0]
        i_end = np.where(time <= t_end)[0][-1]
        axis[0].set_xlim([time[i_start], time[i_end]])

        # find where both valves are closed: iso-volume
        iso_vol = (pressure[0] > pressure[1]) & (pressure[0] < pressure[2])
        # start and end time when this state change
        tt = time[np.where(iso_vol[:-1] != iso_vol[1:])[0]]
        # plot iso-volume phase in grey zone
        try:  # because reshape can fail
            tt = tt.reshape(-1, 2)
            for i in range(len(tt)):
                for j in range(3):
                    axis[j].axvspan(tt[i, 0], tt[i, 1], facecolor="grey", alpha=0.3)
        except:
            print("Cannot find iso-volume stage automatically")

        # do plot
        for pname, pp in zip(p_names, pressure):
            axis[0].plot(time, pp, label=pname)
        axis[0].set_ylabel("Pressure (kPa)")
        axis[0].legend()

        for qname, qq in zip(q_names, flow):
            axis[1].plot(time, qq, label=qname)
        axis[1].set_ylabel("Flow (mL/s)")
        axis[1].legend()

        axis[2].plot(self.bin.time, self.bin.volume[:, id], label="Volume")
        axis[2].set_ylabel("Volume (mL)")
        axis[2].set_xlabel("Time (s)")

        return fig

    def __check_prefilling(self, cavity, offset=0.0, prefill_duration=1000):
        """
        Check prefilling process.

        Obsolete
        """
        fig, axis = plt.subplots(3, figsize=(8, 4), sharex=True)

        if cavity == "lv":
            id = 0
            if self.type == "BV" and self.closed_loop:

                pven0 = self.lv["p_pven"]
            else:
                pven0 = self.lv["pven"]

            fig.suptitle("Left ventricle: Prefilling Check ")
        elif cavity == "rv":
            id = 1
            pven0 = self.rv["pven"]
            fig.suptitle("Right ventricle: Prefilling Check")

        # define the last step
        i_end = np.where(self.bin.time <= prefill_duration + offset)[0][-1] + 1
        i_end2 = np.where(self.lv["time"] <= prefill_duration + offset)[0][-1] + 1

        axis[0].plot(self.bin.time[0:i_end], self.bin.pressure[0:i_end, id])
        axis[0].plot(self.lv["time"][0:i_end2], pven0[0:i_end2], "--", color="orange")
        axis[0].set_ylabel("Pressure (kPa)")

        axis[1].plot(self.bin.time[0:i_end], self.bin.volume[0:i_end, id])
        axis[1].hlines(self.v_ed[id], 0, prefill_duration, linestyles="--", color="orange")
        axis[1].set_ylabel("Volume (mL)")

        axis[2].plot(self.bin.time[0:i_end], self.bin.flow[0:i_end, id])
        axis[2].hlines(0, 0, prefill_duration, linestyles="--", color="orange")
        axis[2].set_ylabel("Flow (mL/s)")
        axis[2].set_xlabel("Time (s)")

        return fig

    def _check_output(self, cavity="lv"):
        """Check system states == FEM states."""
        fig, axis = plt.subplots(2, figsize=(8, 4))

        if cavity == "lv":
            system = self.lv
            id = 0
            fig.suptitle("Left ventricle: Output Check ")

        elif cavity == "rv":
            system = self.rv
            id = 1
            fig.suptitle("Right ventricle: Output Check")
        p_name = "pk"
        q_name = "qk"
        if self.type == "BV" and self.closed_loop:
            p_name = "pvv"
            q_name = "qvv"
        axis[0].plot(self.bin.time, self.bin.pressure[:, id], label="P_binout")
        axis[0].plot(system["time"], system[p_name], "--", label="P_sys")
        axis[0].legend()

        axis[1].plot(self.bin.time, -self.bin.flow[:, id], label="Q_binout")
        axis[1].plot(system["time"], system[q_name], "--", label="Q_sys")

        axis[1].set_xlabel("Time (s)")
        axis[1].legend()

        # a special case where cavity volume is integrated in System
        if "vlv" in system.columns:
            fig2, axis = plt.subplots(figsize=(8, 4))
            fig2.suptitle("Cavity Volume Check ")
            axis.plot(self.bin.time, self.bin.volume[:, id], label="V_binout")
            axis.plot(system["time"], system["vlv"], "--", label="V_sys")
            axis.legend()

        return fig

    def _check_total_volume(self, plot_all=False):
        """Check if total volume is constant for a closed loop.

        Parameters
        ----------
        plot_all

        """
        if not self.closed_loop:
            print("Only closed loop can check volume")
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
