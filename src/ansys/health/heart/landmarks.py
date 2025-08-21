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

"""Module containing heart anatomical landmarks."""

from dataclasses import dataclass, field

from ansys.health.heart.objects import Point


@dataclass
class LandMarks:
    """Heart anatomical points."""

    sa_node: Point = field(default_factory=lambda: Point("SA_node", xyz=None, node_id=None))
    """Sinoatrial node."""
    av_node: Point = field(default_factory=lambda: Point("AV_node", xyz=None, node_id=None))
    """Atrioventricular node."""
    his_bif_node: Point = field(
        default_factory=lambda: Point("His_bifurcation", xyz=None, node_id=None)
    )
    """His bundle bifurcation node."""
    his_left_end_node: Point = field(
        default_factory=lambda: Point("His_left_end", xyz=None, node_id=None)
    )
    """His bundle left end node."""
    his_right_end_node: Point = field(
        default_factory=lambda: Point("His_right_end", xyz=None, node_id=None)
    )
    """His bundle right end node."""
    bachmann_end_node: Point = field(
        default_factory=lambda: Point("Bachmann_end", xyz=None, node_id=None)
    )
    """Bachmann bundle end node."""
    left_fascicle_end_node: Point = field(
        default_factory=lambda: Point("Left_fasicle_end", xyz=None, node_id=None)
    )
    """Left fascicle end node."""
    left_apex: Point = field(default_factory=lambda: Point("Left_apex", xyz=None, node_id=None))
    """Left ventricle apex (endocardium)."""
    right_apex: Point = field(default_factory=lambda: Point("Right_apex", xyz=None, node_id=None))
    """Right ventricle apex (endocardium)."""
