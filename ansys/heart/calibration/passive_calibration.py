"""Calibration passive material parameter by Klotz curve."""
import os
import subprocess

from ansys.heart.postprocessor.Klotz_curve import EDPVR
from ansys.heart.postprocessor.binout_helper import NodOut
from ansys.heart.postprocessor.compute_volume import get_cavity_volume2
import matplotlib.pyplot as plt
import numpy as np


class PassiveCalibration:
    """Passive calibration."""

    def __init__(self, work_directory):
        """
        Initialize Passive Calibration class.

        Parameters
        ----------
        work_directory
        """
        self.work_directory = work_directory

        nodes_file = os.path.join(work_directory, "nodes.k")
        sgm_file = os.path.join(work_directory, "left_ventricle.segment")

        # load nodes.k: EOD geometry
        data = []
        with open(nodes_file) as f:
            for line in f.readlines():
                if line[0] != "*" and line[0] != "$":
                    data.append(line)
        x_ed = np.genfromtxt(data, delimiter=[8, 16, 16, 16])
        x_ed = x_ed[x_ed[:, 0].argsort()][:, 1:]

        # load left cavity segment
        self.lv_cavity = np.loadtxt(sgm_file, delimiter=",", dtype=int)

        # compute Klotz curve
        self.v_ed = get_cavity_volume2(x_ed, self.lv_cavity)
        self.p_ed = 2 * 7.5  # 2kPa to  mmHg
        self.klotz = EDPVR(self.v_ed, self.p_ed)

    def load_results(self):
        """Load zerop simulation results."""
        # load inflation simulation
        # todo: filename iter3 is hard coded
        nodout = NodOut(os.path.join(self.work_directory, "iter3.binout0000"))
        time = nodout.time
        coords = nodout.get_coordinates()

        # compute volume at different simulation time
        self.v_sim = np.zeros(time.shape)
        for i, coord in enumerate(coords):
            lv_volume = get_cavity_volume2(coord, self.lv_cavity)
            self.v_sim[i] = lv_volume

    def compute_error(self):
        """
        Compute rsm error between Klotz curve and simulation.

        Returns
        -------
        float
        error
        """
        self.pressure = np.linspace(0, self.p_ed, num=len(self.v_sim))
        self.v_kz = self.klotz.get_volume(self.pressure)

        # RSM error between analytics and simulation
        self.rsm = np.linalg.norm(self.v_sim - self.v_kz)

        return self.rsm

    def plot(self):
        """
        Plot Klotz and simulation EDPVR.

        Returns
        -------
        plot curve
        """
        vv = np.linspace(0, 1.1 * self.v_ed, num=101)
        pp = self.klotz.get_pressure(vv)
        fig = plt.figure()
        plt.plot(vv, pp, label="Klotz")
        plt.plot(self.v_kz, self.pressure, "o", label="Klotz")
        plt.plot(self.v_sim, self.pressure, "--*", label="FEM")
        plt.legend()
        plt.title("RSM={0:10.5e}".format(self.rsm))
        return fig

    def run_lsdyna(self):
        """
        Run lsdyna in wsl.

        Returns
        -------
        str
        wsl output
        """
        # TODO: cooperate with simulator.
        command = ["wsl", "-e", "bash", "-lic", "./run_lsdyna.sh"]
        run_command_display = " ".join([str(s) for s in command])
        process = subprocess.run(
            ["powershell", "-Command", run_command_display], capture_output=True
        )
        return process.stdout.decode()


if __name__ == "__main__":
    pass
