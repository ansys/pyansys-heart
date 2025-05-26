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

"""Define classes for anatomical parts.

Each class extends the base Part class and provides specialized attributes for different
heart structures.
"""

from enum import Enum

import numpy as np
import yaml

from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.objects import Cap, Cavity, Point, SurfaceMesh
from ansys.health.heart.settings.material.ep_material import EPMaterial
from ansys.health.heart.settings.material.material import MechanicalMaterialModel


class _PartType(Enum):
    """Stores valid part types."""

    VENTRICLE = "ventricle"
    ATRIUM = "atrium"
    SEPTUM = "septum"
    ARTERY = "artery"
    MYOCARDIUM = "myocardium"
    UNDEFINED = "undefined"


class Part:
    """Base part class."""

    @property
    def surfaces(self) -> list[SurfaceMesh]:
        """List of surfaces belonging to the part."""
        surfaces = []
        for key, value in self.__dict__.items():
            if isinstance(value, SurfaceMesh):
                surfaces.append(value)
        return surfaces

    @property
    def surface_names(self) -> list[str]:
        """List of surface names belonging to the part."""
        surface_names = []
        for surface in self.surfaces:
            surface_names.append(surface.name)
        return surface_names

    def get_point(self, pointname: str) -> Point | None:
        """Get a point from the part."""
        for point in self.points:
            if point.name == pointname:
                return point
        LOGGER.error("Cannot find point {0:s}.".format(pointname))
        return None

    def __init__(self, name: str = None, part_type: _PartType = _PartType.UNDEFINED) -> None:
        self.name: str = name
        """Part name."""
        self.pid: int | None = None
        """Part ID."""
        self.mid: int | None = None
        """Material ID associated with the part."""
        self._part_type: _PartType = part_type
        """Type of the part."""
        self.element_ids: np.ndarray = np.empty((0, 4), dtype=int)
        """Array holding element IDs that make up the part."""
        self.points: list[Point] = []
        """Points of interest belonging to the part."""

        self.fiber: bool = False
        """Flag indicating if the part has fiber/sheet data."""
        self.active: bool = False
        """Flag indicating if active stress is established."""

        self.meca_material: MechanicalMaterialModel = MechanicalMaterialModel.DummyMaterial()
        """Material model to assign in the simulator."""

        self.ep_material: EPMaterial = EPMaterial.DummyMaterial()
        """EP material model to assign in the simulator."""

    def __str__(self) -> str:
        """Return a string representation of the part."""
        return yaml.dump(self._to_dict(), indent=4)

    def _to_dict(self) -> dict:
        """Get part information to reconstruct from a mesh file."""
        info = {
            self.name: {
                "part-id": self.pid,
                "part-type": self._part_type.value,
                "surfaces": {},
            }
        }

        info2 = {}
        info2["surfaces"] = {}

        for surface in self.surfaces:
            if isinstance(surface, SurfaceMesh):
                if surface.id:
                    info2["surfaces"][surface.name] = surface.id

        if hasattr(self, "caps"):
            info2["caps"] = {}
            for cap in self.caps:
                info2["caps"][cap.name] = cap._mesh.id

        if hasattr(self, "cavity"):
            info2["cavity"] = {}
            if self.cavity is not None:
                info2["cavity"][self.cavity.surface.name] = self.cavity.surface.id

        info[self.name].update(info2)

        return info


class Septum(Part):
    """Septum part."""

    def __init__(self, name: str = None) -> None:
        super().__init__(name=name, part_type=_PartType.SEPTUM)


class Chamber(Part):
    """Intermediate class for heart chambers with endocardium and epicardium."""

    def __init__(self, name: str = None, part_type: _PartType = None) -> None:
        super().__init__(name=name, part_type=part_type)
        self.endocardium: SurfaceMesh = SurfaceMesh(name=f"{self.name} endocardium")
        """Endocardial surface."""
        self.epicardium: SurfaceMesh = SurfaceMesh(name=f"{self.name} epicardium")
        """Epicardial surface."""

        self.myocardium: Myocardium = Myocardium(name="myocardium")
        """Myocardial part."""

        self.caps: list[Cap] = []
        """List of caps belonging to the part."""
        self.cavity: Cavity | None = None
        """Cavity belonging to the part."""

        self.active: bool = True
        """Flag indicating if active stress should be included."""
        self.fiber: bool = True
        """Flag indicating if fiber/sheet data should be included."""


class Ventricle(Chamber):
    """Ventricle part."""

    def __init__(self, name: str = None) -> None:
        super().__init__(name=name, part_type=_PartType.VENTRICLE)

        self.septum: SurfaceMesh = SurfaceMesh(name="{0} endocardium septum".format(self.name))
        """Septal surface."""

        self.apex_points: list[Point] = []
        """List of apex points."""


class Atrium(Chamber):
    """Atrium part."""

    def __init__(self, name: str = None) -> None:
        super().__init__(name=name, part_type=_PartType.ATRIUM)


class Artery(Part):
    """Artery part."""

    def __init__(self, name: str = None) -> None:
        super().__init__(name=name, part_type=_PartType.ARTERY)

        self.wall: SurfaceMesh = SurfaceMesh(name="{0} wall".format(self.name))


class Myocardium(Part):
    """Myocardium part."""

    def __init__(self, name: str = None) -> None:
        super().__init__(name=name, part_type=_PartType.MYOCARDIUM)
