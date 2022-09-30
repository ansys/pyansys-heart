"""Calibration passive material parameter by Klotz curve."""
import os
import pathlib

from ansys.heart.postprocessor.Klotz_curve import EDPVR
from ansys.heart.postprocessor.compute_volume import get_cavity_volume2
from ansys.heart.postprocessor.read_binout import read_nodout
import matplotlib.pyplot as plt
import numpy as np


class PassiveCalibration:
    """Passive calibration."""

    def __init__(self, folder):
        """
        Initialize Passive Calibration class.

        Parameters
        ----------
        folder
        """
        nodes_file = os.path.join(folder, "nodes.k")
        sgm_file = os.path.join(folder, "left_ventricle.segment")
        # todo: l2a
        bin_file = os.path.join(folder, "nodout")

        # load nodes.k: EOD geometry
        data = []
        with open(nodes_file) as f:
            for line in f.readlines():
                if line[0] != "*" and line[0] != "$":
                    data.append(line)
        x_ed = np.genfromtxt(data, delimiter=[8, 16, 16, 16])
        x_ed = x_ed[x_ed[:, 0].argsort()][:, 1:]

        # load left cavity segment
        lv_cavity = np.loadtxt(sgm_file, delimiter=",", dtype=int)

        # load inflation simulation
        time, coords = read_nodout(bin_file)

        # compute EOD volume and pressure
        self.v_ed = get_cavity_volume2(x_ed, lv_cavity)
        self.p_ed = 2 * 7.5

        self.klotz = EDPVR(self.v_ed, self.p_ed)

        # compute volume at different simulation time
        self.v_sim = np.zeros(time.shape)
        for i, coord in enumerate(coords):
            lv_volume = get_cavity_volume2(coord, lv_cavity)
            self.v_sim[i] = lv_volume

        self.pressure = np.linspace(0, self.p_ed, num=len(self.v_sim))
        self.v_kz = self.klotz.get_volume(self.pressure)

        # RSM error between analytics and simulation
        self.rsm = np.linalg.norm(self.v_sim - self.v_kz)
        return

    def plot(self):
        """
        Plot Klotz and simulation EDPVR.

        Returns
        -------
        plot curve
        """
        plt.plot(self.v_kz, self.pressure, "--o", label="Klotz")
        plt.plot(self.v_sim, self.pressure, "--*", label="FEM")
        plt.legend("RSM={0:10.5e}".format(self.rsm))
        plt.show()
        return


if __name__ == "__main__":
    "test"
    path_to_case = os.path.join(pathlib.Path(__file__).parents[3], "tests", "heart", "calibration")
    case = PassiveCalibration(path_to_case)
    print(case.rsm)
    case.plot()
