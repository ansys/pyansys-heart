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
"""Module to system model."""

from dataclasses import dataclass

from ansys.health.heart.models import BiVentricle, FourChamber, HeartModel, LeftVentricle
from ansys.health.heart.parts import Chamber
from ansys.health.heart.writer.define_function_templates import _define_function_0d_system


@dataclass
class CVInteraction:
    """Template to define control volume interaction."""

    id: int
    cvid1: int
    cvid2: int
    lcid: int
    flow_name: str
    parameters: dict

    def _define_function_keyword(self):
        if self.flow_name == "closed-loop":
            return ""
        else:
            return _define_function_0d_system(self.lcid, self.flow_name, self.parameters)


@dataclass
class ControlVolume:
    """Template to define control volume."""

    part: Chamber
    id: int
    Interactions: list[CVInteraction]


def _create_open_loop(id_offset: int, model: HeartModel, settings) -> list[ControlVolume]:
    """Create open loop system model.

    Parameters
    ----------
    id_offset : int
        ID of the first control volume
    model : HeartModel
        Heart model
    settings : _type_
        parameters for the control volume

    Returns
    -------
    list[ControlVolume]
        list of control volumes
    """
    if isinstance(model, LeftVentricle):
        system_map = [
            ControlVolume(
                part=model.left_ventricle,
                id=1,
                Interactions=[
                    CVInteraction(
                        id=1,
                        cvid1=1,
                        cvid2=0,
                        lcid=id_offset,
                        flow_name="constant_preload_windkessel_afterload_left",
                        parameters=settings.left_ventricle,
                    )
                ],
            )
        ]
    elif isinstance(model, BiVentricle):
        system_map = [
            ControlVolume(
                part=model.left_ventricle,
                id=1,
                Interactions=[
                    CVInteraction(
                        id=1,
                        cvid1=1,
                        cvid2=0,
                        lcid=id_offset,
                        flow_name="constant_preload_windkessel_afterload_left",
                        parameters=settings.left_ventricle,
                    )
                ],
            ),
            ControlVolume(
                part=model.right_ventricle,
                id=2,
                Interactions=[
                    CVInteraction(
                        id=2,
                        cvid1=2,
                        cvid2=0,
                        lcid=id_offset + 1,
                        flow_name="constant_preload_windkessel_afterload_right",
                        parameters=settings.right_ventricle,
                    )
                ],
            ),
        ]
    elif isinstance(model, FourChamber):
        system_map = [
            ControlVolume(
                part=model.left_ventricle,
                id=1,
                Interactions=[
                    CVInteraction(
                        id=1,
                        cvid1=1,
                        cvid2=0,
                        lcid=id_offset,
                        flow_name="afterload_windkessel_left",
                        parameters=settings.left_ventricle,
                    ),
                ],
            ),
            ControlVolume(
                part=model.right_ventricle,
                id=2,
                Interactions=[
                    CVInteraction(
                        id=2,
                        cvid1=2,
                        cvid2=0,
                        lcid=id_offset + 1,
                        flow_name="afterload_windkessel_right",
                        parameters=settings.right_ventricle,
                    ),
                ],
            ),
            ControlVolume(
                part=model.left_atrium,
                id=3,
                Interactions=[
                    CVInteraction(
                        id=3,
                        cvid1=3,
                        cvid2=0,
                        lcid=id_offset + 2,
                        flow_name="constant_flow_left_atrium",
                        parameters={"flow": -83.0},  # ~5 L/min
                    ),
                    CVInteraction(
                        id=4,
                        cvid1=3,
                        cvid2=1,
                        lcid=id_offset + 3,
                        flow_name="valve_mitral",
                        parameters={"Rv": 1e-6},
                    ),
                ],
            ),
            ControlVolume(
                part=model.right_atrium,
                id=4,
                Interactions=[
                    CVInteraction(
                        id=5,
                        cvid1=4,
                        cvid2=0,
                        lcid=id_offset + 4,
                        flow_name="constant_flow_right_atrium",
                        parameters={"flow": -83.0},  # ~5 L/min
                    ),
                    CVInteraction(
                        id=6,
                        cvid1=4,
                        cvid2=2,
                        lcid=id_offset + 5,
                        flow_name="valve_tricuspid",
                        parameters={"Rv": 1e-6},
                    ),
                ],
            ),
        ]

    return system_map


def _create_closed_loop(model: HeartModel) -> list[ControlVolume]:
    """Create close loop system model.

    Parameters
    ----------
    model : HeartModel
        Heart model

    Returns
    -------
    list[ControlVolume]
        list of control volumes
    """
    interaction_id = [-1, -2, -3, -4]

    if isinstance(model, LeftVentricle):
        control_volumes = [model.left_ventricle]
    elif isinstance(model, BiVentricle):
        control_volumes = [model.left_ventricle, model.right_ventricle]
    elif isinstance(model, FourChamber):
        control_volumes = [
            model.left_ventricle,
            model.right_ventricle,
            model.left_atrium,
            model.right_atrium,
        ]

    system_map = []
    for i, part in enumerate(control_volumes):
        system_map.append(
            ControlVolume(
                part=part,
                id=i + 1,
                Interactions=[
                    CVInteraction(
                        id=i + 1,
                        cvid1=i + 1,
                        cvid2=0,
                        lcid=interaction_id[i],
                        flow_name="closed-loop",
                        parameters={},
                    )
                ],
            )
        )
    return system_map
