"""Calibration active material parameters."""
import os
import shutil
import sys

from ansys.heart.general import run_lsdyna
from ansys.heart.postprocessor.SystemModelPost import SystemModelPost
import pkg_resources


class ActiveCalibration:
    """ActiveCalibration."""

    def __init__(self, work_directory):
        """
        Init.

        Parameters
        ----------
        work_directory
        """
        self.work_directory = work_directory
        self.system = None

    def create_calibration_folder(self, python_exe: str = ""):
        """
        Create necessary files for calibration.

        Parameters
        ----------
        python_exe: venv path

        """
        # todo: test if it's a legitimate folder
        here = os.path.dirname(os.path.realpath(__file__))
        shutil.copy(
            pkg_resources.resource_filename("ansys.heart.calibration", "ActiveCalibration.lsopt"),
            os.path.join(self.work_directory, "ActiveCalibration.lsopt"),
        )

        shutil.copy(
            pkg_resources.resource_filename("ansys.heart.calibration", "material_active.k"),
            os.path.join(self.work_directory, "material.k"),
        )

        if python_exe == "":
            python_exe = f"{sys.prefix}\\Scripts\\python.exe"

        with open(os.path.join(self.work_directory, "run.bat"), "w") as f:
            f.write(f"{python_exe} run.py")

        shutil.copy(
            pkg_resources.resource_filename("ansys.heart.calibration", "run_active.template"),
            os.path.join(self.work_directory, "run.py"),
        )

    def load_results(self):
        """Load simulation results."""
        self.system = SystemModelPost(self.work_directory)
        fig = self.system.plot_pv_loop()
        fig.savefig("pv_loop.png")

    def compute_objective_function(self):
        """Define objective function."""
        return abs(self.system.lv_ef - 0.6)

    def run_one_step_calibration(self):
        """Run one step simulation."""
        run_lsdyna("main.k")
        self.load_results()
        error = self.compute_objective_function()
        with open("result", "a") as f:
            f.write("{0:10.5e}".format(error))
