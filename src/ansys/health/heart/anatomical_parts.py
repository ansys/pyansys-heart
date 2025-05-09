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

"""Defines anatomical heart part classes such as Ventricle, Atrium, Septum, Artery, and Myocardium.

These classes extend the base Part class and provide specialized attributes for different
heart structures.
"""

from ansys.health.heart.objects import Cap, Cavity, Part, PartType, Point, SurfaceMesh


class Septum(Part):
    """Septum part."""

    def __init__(self, name: str = None):
        super().__init__(name=name, part_type=PartType.SEPTUM)


class _Chamber(Part):
    """Intermediate class for heart chambers with endocardium and epicardium."""

    def __init__(self, name: str = None, part_type: PartType = None):
        super().__init__(name=name, part_type=part_type)
        self.endocardium = SurfaceMesh(name=f"{self.name} endocardium")
        """Endocardial surface."""
        self.epicardium = SurfaceMesh(name=f"{self.name} epicardium")
        """Epicardial surface."""
        self.septum = SurfaceMesh(name="{0} endocardium septum".format(self.name))
        """Septal part."""

        self.myocardium = Part(name="myocardium", part_type=PartType.MYOCARDIUM)
        """Myocardial part."""

        self.caps: list[Cap] = []
        """List of caps belonging to the part."""
        self.cavity: Cavity = None
        """Cavity belonging to the part."""

        self.active: bool = True
        """Flag indicating if active stress should be included."""
        self.fiber: bool = True
        """Flag indicating if fiber/sheet data should be included."""


class Ventricle(_Chamber):
    """Ventricle part."""

    def __init__(self, name: str = None):
        super().__init__(name=name, part_type=PartType.VENTRICLE)

        self.apex_points: list[Point] = []
        """List of apex points."""


class Atrium(_Chamber):
    """Atrium part."""

    def __init__(self, name: str = None):
        super().__init__(name=name, part_type=PartType.ATRIUM)


class Artery(Part):
    """Artery part."""

    def __init__(self, name: str = None):
        super().__init__(name=name, part_type=PartType.ARTERY)

        self.wall = SurfaceMesh(name="{0} wall".format(self.name))


class Myocardium(Part):
    """Myocardium part."""

    def __init__(self, name: str = None):
        super().__init__(name=name, part_type=PartType.MYOCARDIUM)
