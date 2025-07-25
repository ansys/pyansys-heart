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

from __future__ import annotations

from enum import Enum

from deprecated import deprecated
import numpy as np
import yaml

from ansys.health.heart import LOG as LOGGER, __version__
from ansys.health.heart.objects import Cap, CapType, Cavity, Mesh, Point, SurfaceMesh
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

    @property
    def _attribute_name(self):
        """Equivalent attribute name of the part."""
        return self.name.lower().replace(" ", "_").replace("-", "_")

    def get_point(self, pointname: str) -> Point | None:
        """Get a point from the part."""
        for point in self.points:
            if point.name == pointname:
                return point
        LOGGER.error("Cannot find point {0:s}.".format(pointname))
        return None

    def get_element_ids(self, mesh: Mesh = None) -> np.ndarray:
        """Get element IDs that make up the part.

        Parameters
        ----------
        mesh : Mesh, default: None
            The mesh object where to get the element IDs from.

        Returns
        -------
        np.ndarray
            Array of element IDs that make up the part.
        """
        if mesh is None:
            LOGGER.error("Mesh is not provided to get element IDs.")
            return np.empty((0,), dtype=int)

        if self.pid is None:
            LOGGER.error("Part ID is not set. Cannot get element IDs.")
            return np.empty((0,), dtype=int)

        if "_volume-id" not in mesh.cell_data.keys():
            LOGGER.error("Mesh does not contain '_volume-id' cell data.")
            return np.empty((0,), dtype=int)

        return np.argwhere(mesh.cell_data["_volume-id"] == self.pid).flatten()

    @property
    @deprecated(
        """`element_ids` as an attribute is deprecated. Use `part.get_element_ids(mesh)` instead.
        To modify element IDs of a part use the `_volume-id` cell data of the mesh object.""",
    )
    def element_ids(self) -> np.ndarray:
        """Get element IDs that make up the part."""
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
        data = {
            self.name: {
                "part-id": self.pid,
                "part-type": self._part_type.value,
                "fiber": self.fiber,
                "active": self.active,
                "_version": __version__,
                "surfaces": {},
            }
        }

        data2 = {}
        data2["surfaces"] = {}

        for surface in self.surfaces:
            if isinstance(surface, SurfaceMesh):
                if surface.id:
                    data2["surfaces"][surface.name] = surface.id

        if hasattr(self, "caps"):
            data2["caps"] = {}
            for cap in self.caps:
                data2["caps"][cap.name] = cap._mesh.id

        if hasattr(self, "cavity"):
            data2["cavity"] = {}
            if self.cavity is not None:
                data2["cavity"][self.cavity.surface.name] = self.cavity.surface.id

        data[self.name].update(data2)

        return data

    def _get_predefined_surfaces(self) -> list[str]:
        """Get a list of allowed surfaces for the part."""
        return [
            value.name for key, value in self.__dict__.items() if isinstance(value, SurfaceMesh)
        ]

    @staticmethod
    def _set_from_dict(
        data: dict, mesh: Mesh | None = None
    ) -> Part | Septum | Chamber | Ventricle | Atrium | Artery | Myocardium:
        """Reconstruct a part from a dictionary.

        Parameters
        ----------
        data : dict
            Dictionary that describes the part in JSON format.
        mesh : Mesh, optional
            Mesh object to use for reconstructing surfaces and cavities. If not provided,
            the part is created without surfaces and cavities.

        Returns
        -------
        Part | Septum | Chamber | Ventricle | Atrium | Artery | Myocardium
            Reconstructed part.
        """
        if not isinstance(data, dict):
            LOGGER.error("Data must be a dictionary.")
            return None

        name = next(iter(data), None)

        part_data: dict = data[name]

        if not part_data.get("_version", None):
            import json

            json_str = json.dumps(Part()._to_dict(), indent=4)
            LOGGER.error(
                f"""Part data does not contain version information.
                Consider regenerating in the following format:\n{json_str}."""
            )

        try:
            _part_type: str = _PartType(part_data.get("part-type", _PartType.UNDEFINED.value))
        except ValueError as e:
            LOGGER.error(
                f"""Invalid part type: {part_data.get("part-type", "undefined")}.
                Defaulting to UNDEFINED. {e}"""
            )
            _part_type = _PartType.UNDEFINED

        _part_type_to_class_map = {
            _PartType.SEPTUM: Septum,
            _PartType.VENTRICLE: Ventricle,
            _PartType.ATRIUM: Atrium,
            _PartType.ARTERY: Artery,
            _PartType.MYOCARDIUM: Myocardium,
            _PartType.UNDEFINED: Part,
        }

        part_cls = _part_type_to_class_map[_part_type]
        part: Part = part_cls(name=name)

        # assign part id, active, fiber, and surfaces
        part.pid = part_data.get("part-id", None)
        part.active = part_data.get("active", False)
        part.fiber = part_data.get("fiber", False)

        if mesh:
            # try to set the surfaces and cavities with mesh data.
            for surface_name, surface_id in part_data.get("surfaces", {}).items():
                if surface_name not in part.surface_names:
                    LOGGER.error(
                        f"Surface {surface_name} is not a standard surface for part {name}."
                    )
                    continue
                surface = mesh.get_surface(surface_id)
                attribute_name = surface_name.replace(part.name + " ", "")
                setattr(part, attribute_name, surface)

            # try to initialize cavity object.
            if isinstance(part, Chamber):
                if part_data.get("cavity", {}) != {}:
                    cavity_name, cavity_id = next(iter(part_data.get("cavity").items()))
                    part.cavity = Cavity(surface=mesh.get_surface(cavity_id), name=cavity_name)

                for cap_name, cap_id in part_data.get("caps", {}).items():
                    #! note that we assume cap name equals cap type here.
                    cap = Cap(cap_name, cap_type=CapType(cap_name))
                    cap._mesh = mesh.get_surface(cap_id)
                    part.caps.append(cap)

        return part


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

        self.septum: SurfaceMesh = SurfaceMesh(name="{0} septum".format(self.name))
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

        self.ep_material = EPMaterial.Insulator()
        """EP material model for the artery part."""


class Myocardium(Part):
    """Myocardium part."""

    def __init__(self, name: str = None) -> None:
        super().__init__(name=name, part_type=_PartType.MYOCARDIUM)
