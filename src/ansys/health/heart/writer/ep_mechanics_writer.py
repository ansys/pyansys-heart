# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
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
"""Module containing classes for writing LS-DYNA keyword files for ep-mechanics simulations."""

from typing import Callable, Optional

from ansys.dyna.core.keywords import keywords
from ansys.health.heart.models import FourChamber, HeartModel
from ansys.health.heart.settings.settings import SimulationSettings
from ansys.health.heart.writer.base_writer import BaseDynaWriter
from ansys.health.heart.writer.ep_writer import ElectrophysiologyDynaWriter
from ansys.health.heart.writer.heart_decks import ElectroMechanicsDecks
from ansys.health.heart.writer.mechanics_writer import MechanicsDynaWriter


class ElectroMechanicsDynaWriter(MechanicsDynaWriter, ElectrophysiologyDynaWriter):
    """Class for preparing the input for LS-DYNA electromechanical simulation."""

    def __init__(
        self,
        model: HeartModel,
        settings: Optional[SimulationSettings] = None,
    ) -> None:
        if isinstance(model, FourChamber):
            model._create_atrioventricular_isolation()

        BaseDynaWriter.__init__(self, model=model, settings=settings)

        self.kw_database = ElectroMechanicsDecks()
        """Collection of keyword decks relevant for mechanics."""

        self.set_flow_area = True
        """from MechanicWriter."""

    def update(
        self, dynain_name: Optional[str] = None, robin_bcs: Optional[list[Callable]] = None
    ) -> None:
        """Update the keyword database.

        Parameters
        ----------
        dynain_name : str, default: None
            Dynain file from stress-free configuration computation.
        robin_bcs : list[Callable], default: None
            List of lambda functions to apply Robin-type boundary conditions.

        Notes
        -----
        You do not need to write mesh files if a Dynain file is given.
        """
        if isinstance(self.model, FourChamber):
            self.model.left_atrium.fiber = True
            self.model.left_atrium.active = True
            self.model.right_atrium.fiber = True
            self.model.right_atrium.active = True

        MechanicsDynaWriter.update(self, dynain_name=dynain_name, robin_bcs=robin_bcs)

        if self.model.conduction_mesh.number_of_cells != 0:
            # Coupling enabled, EP beam nodes follow the motion of surfaces
            self.kw_database.ep_settings.append(keywords.EmControlCoupling(thcoupl=1, smcoupl=0))
            beam_pid = self._update_use_Purkinje()
            self.include_to_main("beam_networks.k")
        else:
            beam_pid = None

        self._update_parts_cellmodels()
        self.include_to_main("cell_models.k")

        self._update_ep_settings(beam_pid)
        self._update_stimulation()

        # coupling parameters
        coupling_str = (
            "*EM_CONTROL_COUPLING\n$    THCPL     SMCPL    THLCID    SMLCID\n         1         0\n"
        )
        self.kw_database.ep_settings.append("$ EM-MECA coupling control")
        self.kw_database.ep_settings.append(coupling_str)
        self.include_to_main("ep_settings.k")

        return

    def _update_material_db(self, add_active: bool = True) -> None:
        """Update the database of material keywords."""
        MechanicsDynaWriter._update_material_db(self, add_active=add_active, em_couple=True)
        ElectrophysiologyDynaWriter._update_ep_material_db(self)
        return
