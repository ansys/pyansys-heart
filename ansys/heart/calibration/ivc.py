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

"""Overload mechanical writer and simulator for iso-volumic contraction(IVC) simulation."""

import copy
import os
import pathlib as Path
from typing import Literal

from ansys.heart import LOG as LOGGER
from ansys.heart.preprocessor.models.v0_1.models import HeartModel
from ansys.heart.simulator.settings.settings import SimulationSettings
from ansys.heart.simulator.simulator import MechanicsSimulator
from ansys.heart.writer.dynawriter import MechanicsDynaWriter
from pint import Quantity


class IVCWriter(MechanicsDynaWriter):
    """Overload MechanicsDynaWriter for IVC."""

    def __init__(
        self,
        model: HeartModel,
        settings: SimulationSettings = None,
    ) -> None:
        """Overload MechanicsDynaWriter for IVC."""
        super().__init__(model=model, settings=settings)

        self.settings.mechanics.analysis.end_time = Quantity(1000, "millisecond")

    def _update_system_model(self):
        """Create a system model so valves are closed anyway."""
        system_settings = copy.deepcopy(self.settings.mechanics.system)
        system_settings._remove_units()

        from ansys.heart.writer.system_models import define_function_windkessel

        if self.system_model_name != system_settings.name:
            LOGGER.error("Circulation system parameters cannot be rad from Json")

        for cavity in self.model.cavities:
            if "Left ventricle" in cavity.name:
                define_function_wk = define_function_windkessel(
                    function_id=10,
                    function_name="constant_preload_windkessel_afterload_left",
                    implicit=True,
                    constants=dict(system_settings.left_ventricle["constants"]),
                    initialvalues=system_settings.left_ventricle["initial_value"]["part"],
                    ivc=True,
                )
                self.kw_database.control_volume.append(define_function_wk)

            elif "Right ventricle" in cavity.name:
                define_function_wk = define_function_windkessel(
                    function_id=11,
                    function_name="constant_preload_windkessel_afterload_right",
                    implicit=True,
                    constants=dict(system_settings.right_ventricle["constants"]),
                    initialvalues=system_settings.right_ventricle["initial_value"]["part"],
                    ivc=True,
                )
                self.kw_database.control_volume.append(define_function_wk)


class IVCSimulator(MechanicsSimulator):
    """Overload MechanicsSimulator for IVC."""

    def __init__(
        self,
        model: HeartModel,
        lsdynapath: Path,
        dynatype: Literal["smp", "intelmpi", "platformmpi"],
        num_cpus: int = 1,
        simulation_directory: Path = "",
        initial_stress: bool = True,
    ) -> None:
        """Overload MechanicsSimulator for IVC."""
        super().__init__(
            model, lsdynapath, dynatype, num_cpus, simulation_directory, initial_stress
        )

    def _write_main_simulation_files(self, folder_name):
        """Overload to call IVC writer."""
        export_directory = os.path.join(self.root_directory, folder_name)
        self.directories["main-mechanics"] = export_directory

        dyna_writer = IVCWriter(
            self.model,
            self.settings,
        )
        dyna_writer.update(with_dynain=self.initial_stress)
        dyna_writer.export(export_directory)

        return export_directory
