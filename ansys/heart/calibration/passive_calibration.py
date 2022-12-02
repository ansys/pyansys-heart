"""Calibration passive material parameter by Klotz curve."""
import os
import shutil
import sys

from ansys.heart.general import run_lsdyna
from ansys.heart.postprocessor.Klotz_curve import EDPVR
from ansys.heart.postprocessor.binout_helper import NodOut
from ansys.heart.preprocessor.mesh.objects import Cavity, SurfaceMesh
import matplotlib.pyplot as plt
import numpy as np
import pkg_resources


def create_calibration_folder(target_dir, python_exe: str = ""):
    """
    Create necessary files for calibration.

    Parameters
    ----------
    target_dir: target folder
    python_exe: venv path

    """
    # todo: test if it's a legitimate folder
    here = os.path.dirname(os.path.realpath(__file__))
    file_path = pkg_resources.resource_filename(
        "ansys.heart.calibration", "PassiveCalibration.lsopt"
    )
    shutil.copy(
        file_path,
        os.path.join(target_dir, "PassiveCalibration.lsopt"),
    )

    file_path = pkg_resources.resource_filename("ansys.heart.calibration", "material.k")
    shutil.copy(
        file_path,
        os.path.join(target_dir, "material.k"),
    )

    if python_exe == "":
        python_exe = f"{sys.prefix}\\Scripts\\python.exe"

    with open(os.path.join(target_dir, "run.bat"), "w") as f:
        f.write(f"{python_exe} run.py")

    file_path = pkg_resources.resource_filename("ansys.heart.calibration", "run.template")
    shutil.copy(
        file_path,
        os.path.join(target_dir, "run.py"),
    )


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
        # change ID index start by 0
        lv_cavity_faces = np.loadtxt(sgm_file, delimiter=",", dtype=int) - 1

        # create cavity object by using surface mesh
        self.lv_cavity = Cavity(SurfaceMesh(name="lv_cavity", faces=lv_cavity_faces, nodes=x_ed))
        self.lv_cavity.compute_volume()

        # compute Klotz curve
        self.v_ed = self.lv_cavity.volume
        self.p_ed = 2 * 7.5  # 2kPa to  mmHg
        self.klotz = EDPVR(self.v_ed, self.p_ed)

        self.volume_sim = None
        self.pressure_sim = None

    def load_results(self):
        """Load zerop simulation results."""
        # load inflation simulation
        # todo: filename iter3 is hard coded

        try:
            nodout = NodOut(os.path.join(self.work_directory, "iter3.binout"))
        except IOError:
            nodout = NodOut(os.path.join(self.work_directory, "iter3.binout0000"))
        finally:
            Exception("Cannot load binout file")

        # time need to be normalized
        time = nodout.time / nodout.time[-1]
        coords = nodout.get_coordinates()

        # compute volume at different simulation time
        self.pressure_sim = time * self.p_ed
        self.volume_sim = np.zeros(time.shape)
        for i, coord in enumerate(coords):
            self.lv_cavity.surface.nodes = coord
            self.lv_cavity.compute_volume()
            self.volume_sim[i] = self.lv_cavity.volume

    def compute_error(self):
        """
        Compute rsm error between Klotz curve and simulation.

        Returns
        -------
        float
        error
        """
        v_kz = self.klotz.get_volume(self.pressure_sim)

        # RSM error between analytics and simulation
        rsm = np.linalg.norm(self.volume_sim - v_kz)

        return rsm

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
        plt.plot(self.klotz.get_volume(self.pressure_sim), self.pressure_sim, "o", label="Klotz")
        plt.plot(self.volume_sim, self.pressure_sim, "--*", label="FEM")
        plt.legend()
        plt.title("RSM={0:10.5e}".format(self.compute_error()))
        return fig

    def run_one_step_calibration(self):
        """Run zerop simulation and compare with Klotz curve."""
        run_lsdyna("main.k", options="case")
        self.load_results()
        error = self.compute_error()
        with open("result", "a") as f:
            f.write("{0:10.5e}".format(error))
        fig = self.plot()
        fig.savefig("vs.png")


if __name__ == "__main__":
    pass
