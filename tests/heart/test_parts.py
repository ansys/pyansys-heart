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

import pytest
import pyvista as pv

from ansys.health.heart.objects import Cap, Cavity, SurfaceMesh
from ansys.health.heart.parts import (
    Artery,
    Atrium,
    Myocardium,
    Part,
    Septum,
    Ventricle,
    _PartType,
)


@pytest.mark.parametrize(
    "cls, name, expected_keys",
    [
        (Part, "Part", ["part-id", "part-type", "surfaces"]),
        (Septum, "Septum", ["part-id", "part-type", "surfaces"]),
        (Artery, "Artery", ["part-id", "part-type", "surfaces"]),
        (Myocardium, "Myocardium", ["part-id", "part-type", "surfaces"]),
        (
            Ventricle,
            "Ventricle",
            ["part-id", "part-type", "surfaces", "caps", "cavity"],
        ),
        (Atrium, "Atrium", ["part-id", "part-type", "surfaces", "caps", "cavity"]),
    ],
)
def test_part_get_info(cls, name, expected_keys):
    """Test getting part info for each of the Part classes."""
    part: Part = cls(name)

    info = part._get_info()

    assert list(info[name].keys()) == expected_keys


def test_part_get_info_with_data():
    """Test getting part info when data is present.."""
    # Prepare a mock part with some data.
    part = Ventricle("Part1")
    part.pid = 1
    part.endocardium = SurfaceMesh(pv.Tube(), name="tube1", id=10)
    part.epicardium = SurfaceMesh(pv.Tube(radius=1.2), name="tube2", id=11)

    # Create two mock caps.
    cap1 = Cap("cap1")
    cap1._mesh = SurfaceMesh(pv.Circle(0.2), id=100, name="cap1")
    cap2 = Cap("cap2")
    cap2._mesh = SurfaceMesh(pv.Circle(0.1), id=101, name="cap2")
    part.caps.extend([cap1, cap2])

    # Add a mock cavity.
    part.cavity = Cavity(
        surface=SurfaceMesh(
            pv.merge([part.endocardium, cap1._mesh, cap2._mesh]), name="cavity1", id=1000
        )
    )

    info = part._get_info()

    assert info["Part1"]["part-id"] == 1
    assert info["Part1"]["part-type"] == _PartType.VENTRICLE.value
    assert info["Part1"]["surfaces"] == {"tube1": 10, "tube2": 11}
    assert info["Part1"]["caps"] == {"cap1": 100, "cap2": 101}
    assert info["Part1"]["cavity"] == {"cavity1": 1000}
