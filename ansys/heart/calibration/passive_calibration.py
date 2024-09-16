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

"""Calibration passive material parameter by Klotz curve."""

# from importlib.resources import files
from importlib.resources import path as resource_path
import json
import os
import shutil
import sys
import textwrap

from ansys.heart.postprocessor.Klotz_curve import EDPVR
from ansys.heart.preprocessor.models import HeartModel
import numpy as np

from ansys.heart.simulator.settings import settings
from ansys.heart.simulator.simulator import MechanicsSimulator


class PassiveCalibration:
    """Passive calibration."""

    def __init__(self, work_directory, model_path, settings_file=""):
        """
        Initialize Passive Calibration class.

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
        settings : str
            settings file
        work_directory : str
            work directory
        lsdyna_path : str
            path to lsdyna
        ncpu : int
            number of cpu
        dynatype : str
            dynatype
        model_path : str
            model path
        python_exe : str
            python executable path

        """
        os.makedirs(work_directory)
        # LSOPT project file
        shutil.copy(
            resource_path("ansys.heart.calibration", "PassiveCalibration.lsopt").__enter__(),
            os.path.join(work_directory, "PassiveCalibration.lsopt"),
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
            from ansys.heart.calibration.passive_calibration import PassiveCalibration
            \n
            if __name__ == "__main__":
                path_to_working_directory = pathlib.Path(__file__).absolute().parents[0]
                os.chdir(path_to_working_directory)
                case = PassiveCalibration(path_to_working_directory,
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
        lsdyna_path: str
            path to lsdyna executable
        ncpu: int
            number of cpus
        dynatype: str
            type of lsdyna
        """
        simulator = MechanicsSimulator(
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
        simulator.compute_stress_free_configuration()

        self.define_objective()

    def apply_input_parameter(self, setting: settings):
        """Apply parameters."""
        with open(os.path.join(self.work_directory, "parameters.k")) as f:
            l = f.readlines()
            p1 = float(l[0].split(",")[1])
            p2 = float(l[1].split(",")[1])

        setting.mechanics.material.myocardium["isotropic"]["k1"] *= p1
        setting.mechanics.material.myocardium["isotropic"]["k2"] *= p2
        setting.mechanics.material.myocardium["anisotropic"]["k1f"] *= p1
        setting.mechanics.material.myocardium["anisotropic"]["k2f"] *= p2

        try:
            setting.mechanics.material.myocardium["anisotropic"]["k1s"] *= p1
            setting.mechanics.material.myocardium["anisotropic"]["k2s"] *= p2
            setting.mechanics.material.myocardium["anisotropic"]["k1fs"] *= p1
            setting.mechanics.material.myocardium["anisotropic"]["k2fs"] *= p2
        except:
            pass

    def define_objective(self):
        """Define objective function."""
        with open(
            os.path.join(self.work_directory, "zeropressure", "post", "Post_report.json")
        ) as fjson:
            dct = json.load(fjson)

        v_ed = dct["True left ventricle volume (mm3)"]
        p_ed = dct["Left ventricle EOD pressure (mmHg)"]
        klotz = EDPVR(v_ed / 1000, p_ed)

        v_sim = np.array(dct["Simulation Left ventricle volume (mm3)"]) / 1000
        p_sim = np.array(dct["Simulation output time (ms)"]) * p_ed / 1000

        v_klotz = klotz.get_volume(p_sim)

        with open(os.path.join(self.work_directory, "result"), "w") as f:
            error = np.linalg.norm(v_sim - v_klotz)
            f.write("{0:10.5e}".format(error))
