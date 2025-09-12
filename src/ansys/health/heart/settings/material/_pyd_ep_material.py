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

"""Material module."""

from typing import Literal, Optional

from pydantic import BaseModel, Field, model_validator

from ansys.health.heart.settings.defaults import electrophysiology as ep_defaults
from ansys.health.heart.settings.material.cell_models import Tentusscher, TentusscherEndo


class EPMaterialModel(BaseModel):
    """Base class for all EP material models."""

    sigma_fiber: Optional[float] = None
    sigma_sheet: Optional[float] = None
    sigma_sheet_normal: Optional[float] = None
    beta: Optional[float] = ep_defaults.material["myocardium"]["beta"].m
    cm: Optional[float] = ep_defaults.material["myocardium"]["cm"].m
    lambda_: Optional[float] = None

    @model_validator(mode="after")
    def check_inputs(self):
        """Post init method."""
        if self.sigma_sheet is not None and self.sigma_sheet_normal is None:
            self.sigma_sheet_normal = self.sigma_sheet
        if self.sigma_sheet_normal is not None and self.sigma_sheet is None:
            self.sigma_sheet = self.sigma_sheet_normal

        return self


class Insulator(BaseModel):
    """Insulator material."""

    sigma_fiber: float = 0.0
    cm: float = 0.0
    beta: float = 0.0


class Active(EPMaterialModel):
    """Hold data for EP material."""

    solver_type: Literal["Monodomain", "Eikonal", "Reaction-Eikonal"] = ep_defaults.analysis[
        "solvertype"
    ]

    sigma_fiber: Optional[float] = None
    sigma_sheet: Optional[float] = None
    sigma_sheet_normal: Optional[float] = None

    cell_model: Tentusscher = Field(default_factory=lambda: Tentusscher())

    # NOTE: complicated logic and conditional default values. Should split into different classes
    @model_validator(mode="after")
    def check_sigmas(self):
        if self.solver_type == "Monodomain":
            if self.sigma_fiber is None:
                self.sigma_fiber = ep_defaults.material["myocardium"]["sigma_fiber"].m
            if self.sigma_sheet is None:
                self.sigma_sheet = ep_defaults.material["myocardium"]["sigma_sheet"].m
            if self.sigma_sheet_normal is None:
                self.sigma_sheet_normal = ep_defaults.material["myocardium"]["sigma_sheet_normal"].m
        elif self.solver_type == "Eikonal" or self.solver_type == "ReactionEikonal":
            if self.sigma_fiber is None:
                self.sigma_fiber = ep_defaults.material["myocardium"]["velocity_fiber"].m
            if self.sigma_sheet is None:
                self.sigma_sheet = ep_defaults.material["myocardium"]["velocity_sheet"].m
            if self.sigma_sheet_normal is None:
                self.sigma_sheet_normal = ep_defaults.material["myocardium"][
                    "velocity_sheet_normal"
                ].m

        return self


class ActiveBeam(Active):
    """Hold data for beam active EP material."""

    sigma_fiber: float = ep_defaults.material["beam"]["sigma"].m
    cell_model: TentusscherEndo = Field(default_factory=lambda: TentusscherEndo())
    pmjres: float = ep_defaults.material["beam"]["pmjres"].m


class Passive(EPMaterialModel):
    """Hold data for EP passive material."""

    sigma_fiber: float = ep_defaults.material["myocardium"]["sigma_fiber"].m
    sigma_sheet: Optional[float] = None
    sigma_sheet_normal: Optional[float] = None
