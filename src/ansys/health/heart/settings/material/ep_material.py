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

"""EP material module."""

from enum import Enum
from typing import Optional

from pydantic import BaseModel, Field, model_validator

from ansys.health.heart.settings.material.cell_models import Tentusscher


class EPSolverType(Enum):
    """Enumeration of EP solver types."""

    MONODOMAIN = "Monodomain"
    EIKONAL = "Eikonal"
    REACTION_EIKONAL = "ReactionEikonal"


class EPMaterialModel(BaseModel):
    """Base class for all EP material models."""

    sigma_fiber: Optional[float] = None
    sigma_sheet: Optional[float] = None
    sigma_sheet_normal: Optional[float] = None
    beta: Optional[float] = None
    cm: Optional[float] = None
    lambda_: Optional[float] = None

    @model_validator(mode="after")
    def check_inputs(self):
        """Post init method."""
        if self.sigma_sheet is not None and self.sigma_sheet_normal is None:
            self.sigma_sheet_normal = self.sigma_sheet
        if self.sigma_sheet_normal is not None and self.sigma_sheet is None:
            self.sigma_sheet = self.sigma_sheet_normal

        return self


class Insulator(EPMaterialModel):
    """Insulator material."""

    sigma_fiber: float = 0.0
    sigma_sheet_normal: float = 0.0
    sigma_sheet: float = 0.0
    cm: float = 0.0
    beta: float = 0.0


class Passive(EPMaterialModel):
    """Hold data for a passive EP material."""


class Active(EPMaterialModel):
    """Hold data for an active EP material."""

    # NOTE: an active EP material has a cell model associated with it.

    cell_model: Tentusscher = Field(default_factory=lambda: Tentusscher())


class ActiveBeam(Active):
    """Hold data for beam active EP material."""

    # TODO: replace by TentusscherEndo
    cell_model: Tentusscher = Tentusscher()

    @model_validator(mode="after")
    def check_inputs(self):
        """Post init method."""
        # ActiveBeam is by definition isotropic, so remove sheet conductivities if set
        self.sigma_sheet is None
        self.sigma_sheet_normal is None

        return self
