"""Calibration active material parameters."""
# from importlib.resources import files
from importlib.resources import path as resource_path
import os
import pathlib
import shutil
import sys
import textwrap

import numpy as np

from ansys.heart import LOG as LOGGER
from ansys.heart.calibration.ivc import IVCSimulator
from ansys.heart.postprocessor.SystemModelPost import SystemModelPost
from ansys.heart.preprocessor.models.v0_1.models import HeartModel
from ansys.heart.simulator.settings import settings


class ActiveCalibration:
    """ActiveCalibration."""

    def __init__(self, work_directory, model_path, settings_file=""):
        """
        Initialize Calibration class.

        Parameters
        ----------
        work_directory : str
            work directory.
        model_path : str
            HeartModel path
        settings_file : str, optional
            yaml or json settings file, None will load package default settings

        Notes
        -----
            Fiber information must be saved in model_path.

            See :func:`heart.simulator.simulator.BaseSimulator.compute_fibers`

        """
        self.work_directory = work_directory
        self.model = HeartModel.load_model(model_path)
        if settings_file == "":
            self.settings_file = None
        else:
            self.settings_file = settings_file

        zerop_folder = os.path.join(pathlib.Path(model_path).parent, "zeropressure")
        if not os.path.isdir(zerop_folder):
            LOGGER.error("A zeropressure must be present.")
            exit()
        else:
            self.zerop_path = zerop_folder

    @staticmethod
    def create_calibration_project(
        work_directory: str,
        lsdyna_path: str,
        ncpu: int,
        dynatype: str,
        model_path: str,
        settings: str = "",
        python_exe: str = f"{sys.prefix}\\Scripts\\python.exe",
    ):
        """
        Create necessary files for calibration.

        Parameters
        ----------
        settings
        work_directory
        lsdyna_path
        ncpu
        dynatype
        model_path
        python_exe

        """
        os.makedirs(work_directory)
        # LSOPT project file

        shutil.copy(
            resource_path("ansys.heart.calibration", "ActiveCalibration.lsopt").__enter__(),
            os.path.join(work_directory, "ActiveCalibration.lsopt"),
        )

        # LSOPT job file
        with open(os.path.join(work_directory, "run.bat"), "w") as f:
            f.write(f"{python_exe} run.py")

        # LSOPT job is to run this run.py
        with open(os.path.join(work_directory, "run.py"), "w") as f:
            f.write(
                textwrap.dedent(
                    f"""
            import os
            import pathlib
            from ansys.heart.calibration.active_calibration import ActiveCalibration
            \n
            if __name__ == "__main__":
                path_to_working_directory = pathlib.Path(__file__).absolute().parents[0]
                os.chdir(path_to_working_directory)
                case = ActiveCalibration(path_to_working_directory,
                                         r"{model_path}",
                                         r"{settings}")
                case.run_one_step_calibration("{lsdyna_path}",{ncpu},"{dynatype}")
            """
                )
            )

        # LSOPT input file
        with open(os.path.join(work_directory, "parameters.k"), "w") as f:
            f.write("ps1,<<ps1>>\n")
            f.write("ps2,<<ps2>>\n")

    def run_one_step_calibration(self, lsdyna_path, ncpu, dynatype):
        """
        Run zerop pressure simulation.

        Parameters
        ----------
        lsdyna_path
        ncpu
        dynatype

        """
        simulator = IVCSimulator(
            self.model,
            lsdyna_path,
            dynatype,
            num_cpus=ncpu,
            simulation_directory=os.path.join(self.work_directory),
        )

        if self.settings_file is None:
            simulator.load_default_settings()
        else:
            simulator.settings.load(self.settings_file)

        self.apply_input_parameter(simulator.settings)
        # run
        simulator.simulate(folder_name="ivc", zerop_folder=self.zerop_path, auto_post=False)

        self.define_objective()

    def apply_input_parameter(self, setting: settings):
        """Apply parameters."""
        with open(os.path.join(self.work_directory, "parameters.k")) as f:
            l = f.readlines()
            p1 = float(l[0].split(",")[1])
            p2 = float(l[1].split(",")[1])
        from pint import Quantity

        # dummy input parameters
        setting.mechanics.material.myocardium["active"]["taumax"] = Quantity(p1, "MPa")
        setting.mechanics.material.myocardium["active"]["ss"] = p2

    def define_objective(self):
        """Define objective function."""
        s = SystemModelPost(os.path.join(self.work_directory, "ivc"))
        t = s.lv_system.time
        p = s.lv_system.pressure.cavity

        # dummy objective function
        error = (np.max(p) - 20) ** 2

        with open(os.path.join(self.work_directory, "result"), "w") as f:
            f.write("{0:10.5e}".format(error))
