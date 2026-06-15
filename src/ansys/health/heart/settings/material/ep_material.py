# Copyright (C) 2023 - 2026 ANSYS, Inc. and/or its affiliates.
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

from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.settings.material.cell_models import Tentusscher, TentusscherEndo


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
    """
    Insulator material.

    Will use ``*EM_MAT_001`` with mtype=1.
    """

    sigma_fiber: float = 0.0
    sigma_sheet_normal: float = 0.0
    sigma_sheet: float = 0.0
    cm: float = 0.0
    beta: float = 0.0


class Passive(EPMaterialModel):
    """
    Hold data for a passive EP material.

    Will use ``*EM_MAT_001`` or ``*EM_MAT_003`` but no cell model associated with it.
    """


class Active(EPMaterialModel):
    """
    Hold data for an active EP material.

    Use ``*EM_MAT_003`` with mtype=2, associated with a cell model.
    """

    cm: Optional[float] = 0.01
    beta: Optional[float] = 140.0

    cell_model: Tentusscher = Field(default_factory=lambda: Tentusscher())


class ActiveNew(Active):
    """
    Same with :class:`Active` except use ``*EM_MAT_005``.

    In Eikonal-reaction model, it uses anisotropic conductivity tensor
    for ECG computation.
    Other advantages in Bi-domain model but not supported in current version.
    """

    cond_sigma_fiber: Optional[float] = None
    cond_sigma_sheet: Optional[float] = None
    cond_sigma_sheet_normal: Optional[float] = None

    @model_validator(mode="after")
    def warn_on_init(self):
        """Warn that ActiveNew is experimental."""
        LOGGER.warning(
            "ActiveNew is initialized. This material uses ``*EM_MAT_005`` which "
            "requires LS-DYNA R16.2 or later."
        )
        return self


class ActiveBeam(Active):
    """
    Hold data for conduction beam's active EP material.

    Will use ``*EM_MAT_001`` with mtype=2, associated with a cell model same to the endocardium.
    """

    cell_model: Tentusscher = TentusscherEndo()

    @model_validator(mode="after")
    def check_inputs(self):
        """Post init method."""
        # ActiveBeam is isotropic, so remove sheet conductivity if set
        self.sigma_sheet = None
        self.sigma_sheet_normal = None

        return self
