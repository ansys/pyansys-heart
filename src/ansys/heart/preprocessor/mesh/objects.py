# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
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

"""
Module that contains classes relevant for the mesh.

Such as a Mesh object, Part object, Features, etc.

"""

from enum import Enum
import pathlib
from typing import List, Union

import numpy as np

from ansys.heart.core import LOG as LOGGER
import ansys.heart.preprocessor.mesh.vtkmethods as vtkmethods
from ansys.heart.simulator.settings.material.ep_material import EPMaterialModel
from ansys.heart.simulator.settings.material.material import MechanicalMaterialModel

try:
    import pyvista as pv
except ImportError:
    LOGGER.warning("Importing pyvista failed. Install with: pip install pyvista")


class Feature:
    """Feature class."""

    def __init__(self, name: str = None) -> None:
        LOGGER.error(DeprecationWarning("Deprecated"))
        self.name = name
        """Name of feature."""
        self.type = None
        """Type of feature."""
        self.nsid: int = None
        """Node set id associated with feature."""
        self.pid: int = None
        """Part id associated with the feature."""

        pass


class SurfaceMesh(pv.PolyData, Feature):
    """Surface class."""

    @property
    def nodes(self):
        """Node coordinates."""
        return np.array(self.points)

    @nodes.setter
    def nodes(self, array: np.ndarray):
        if isinstance(array, type(None)):
            return
        try:
            num_extra_points = array.shape[0] - self.points.shape[0]
            self.points = array
            if num_extra_points > 0:
                for key in self.point_data.keys():
                    shape = self.point_data[key].shape
                    dtype = self.point_data[key].dtype

                    # vectors
                    if len(shape) > 1:
                        append_shape = (num_extra_points, shape[1])
                        self.point_data[key] = np.vstack(
                            [self.point_data[key], np.empty(append_shape, dtype) * np.nan]
                        )
                    # scalars
                    else:
                        append_shape = (num_extra_points,)
                        self.point_data[key] = np.append(
                            self.point_data[key], np.empty(append_shape, dtype) * np.nan
                        )

            elif num_extra_points < 0:
                raise NotImplementedError(
                    "Assigning less nodes than the original not implemented yet."
                )

        except:
            LOGGER.warning("Failed to set nodes.")
            return

    @property
    def triangles(self):
        """Triangular faces of the surface num_faces x 3."""
        faces = np.reshape(self.faces, (self.n_cells, 3 + 1))[:, 1:]
        return faces

    @triangles.setter
    def triangles(self, value: np.ndarray):
        # sets faces of PolyData
        try:
            num_faces = value.shape[0]
            faces = np.hstack([np.full((num_faces, 1), 3, dtype=np.int8), value])
            self.faces = faces
        except:
            return

    @property
    def boundary_edges(self):
        """Get boundary edges of self."""
        boundary_edges = vtkmethods.get_boundary_edge_loops(self, remove_open_edge_loops=False)
        boundary_edges = np.vstack(list(boundary_edges.values()))
        return boundary_edges

    def __init__(
        self,
        name: str = None,
        triangles: np.ndarray = None,
        nodes: np.ndarray = None,
        id: int = None,
    ) -> None:
        super().__init__(self)
        Feature.__init__(self, name)

        self.type = "surface"
        """Surface type."""
        self.id: int = id
        """ID of surface."""
        self.nsid: int = None
        """ID of corresponding set of nodes."""

        self.triangles = triangles
        """Triangular faces of the surface num_faces x 3."""
        self.nodes = nodes
        """Node coordinates."""

    @property
    def node_ids(self) -> np.ndarray:
        """Global node ids - sorted by earliest occurrence."""
        _, idx = np.unique(self.triangles.flatten(), return_index=True)
        node_ids = self.triangles.flatten()[np.sort(idx)]
        return node_ids

    @property
    def _boundary_nodes(self) -> np.ndarray:
        """Global node ids of nodes on the boundary of the mesh (if any)."""
        _, idx = np.unique(self.boundary_edges.flatten(), return_index=True)
        node_ids = self.boundary_edges.flatten()[np.sort(idx)]
        return node_ids

    def write_to_stl(self, filename: pathlib.Path = None) -> None:
        """Write the surface to a vtk file."""
        if not filename:
            filename = "_".join(self.name.lower().split()) + ".stl"
        if filename[-4:] != ".stl":
            filename = filename + ".stl"

        # NOTE: The below should yield the same stls, but somehow fluent meshing
        # produces a slightly different mesh. Should still be valid though
        # cleaned = self.clean()
        # cleaned.save(filename)
        # vtkmethods.add_solid_name_to_stl(filename, self.name, file_type="binary")

        vtk_surface = vtkmethods.create_vtk_surface_triangles(self.nodes, self.triangles)
        vtkmethods.vtk_surface_to_stl(vtk_surface, filename, self.name)
        return


class BeamMesh(pv.UnstructuredGrid, Feature):
    """Beam class."""

    all_beam_nodes = []
    # beam nodes array

    @property
    def nodes(self):
        """Node coordinates."""
        return np.array(self.points)

    @nodes.setter
    def nodes(self, array: np.ndarray):
        if isinstance(array, type(None)):
            return
        try:
            self.points = array
        except:
            LOGGER.warning("Failed to set nodes.")
            return

    @property
    def edges(self):
        """Tetrahedrons num_tetra x 4."""
        return self.cells_dict[pv.CellType.LINE]

    @edges.setter
    def edges(self, value: np.ndarray):
        # sets lines of UnstructuredGrid
        try:
            points = self.points
            celltypes = np.full(value.shape[0], pv.CellType.LINE, dtype=np.int8)
            lines = np.hstack([np.full(len(celltypes), 2)[:, None], value])
            super().__init__(lines, celltypes, points)
        except:
            LOGGER.warning("Failed to set lines.")
            return

    def __init__(
        self,
        name: str = None,
        edges: np.ndarray = None,
        nodes: np.ndarray = None,
        beam_nodes_mask: np.ndarray = None,
        pid: int = None,
        nsid: int = -1,
    ) -> None:
        super().__init__(self)
        Feature.__init__(self, name)

        self.edges = edges
        """Beams edges."""

        self.nodes = nodes
        """Node coordinates."""

        self.beam_nodes_mask = beam_nodes_mask
        """True for beam nodes, False for solid nodes."""

        self.pid = pid
        """Part id associated with the network."""

        self.nsid: int = nsid
        """Surface id associated with the network."""

        self._all_beam_nodes: np.ndarray = np.empty((0, 3))
        """Temporary attribute to save all previously created beam nodes."""

        self.ep_material: EPMaterialModel = EPMaterialModel.DummyMaterial()
        """Initialize dummy ep material model"""


class Cavity(Feature):
    """Cavity class."""

    def __init__(self, surface: SurfaceMesh = None, centroid: np.ndarray = None, name=None) -> None:
        super().__init__(name)

        self.surface: SurfaceMesh = surface
        """Surface mesh making up the cavity."""
        self.centroid: np.ndarray = centroid
        """Centroid of the cavity."""

    @property
    def volume(self):
        """Volume of the cavity."""
        return self.surface.volume

    def compute_centroid(self):
        """Compute the centroid of the cavity."""
        # self.centroid = np.mean(self.surface.nodes[np.unique(self.surface.triangles), :], axis=0)
        self.centroid = self.surface.center
        return self.centroid


class Cap(Feature):
    """Cap class."""

    def __init__(self, name: str = None, node_ids: Union[List[int], np.ndarray] = []) -> None:
        super().__init__(name)
        self.node_ids = node_ids
        """(Global) node ids of the cap."""
        self.triangles = None
        """Triangulation of cap."""
        self.normal = None
        """Normal of cap."""
        self.centroid = None
        """Centroid of cap."""
        self.centroid_id = None
        """Centroid of cap ID (in case centroid node is created)."""
        return


class Point(Feature):
    """Point class. Can be used to collect relevant points in the mesh."""

    def __init__(self, name: str = None, xyz: np.ndarray = None, node_id: int = None) -> None:
        super().__init__(name)

        self.xyz: np.ndarray = xyz
        """XYZ Coordinates of point."""
        self.node_id: int = node_id
        """Global node id of point."""


class Mesh(pv.UnstructuredGrid):
    """Mesh class: inherits from pyvista UnstructuredGrid.

    Notes
    -----
    Only tetrahedrons are supported.
    Additional attributes are added on top of the pyvista UnstructuredGrid class
    """

    @property
    def nodes(self):
        """Node coordinates."""
        return np.array(self.points)

    @nodes.setter
    def nodes(self, array: np.ndarray):
        if isinstance(array, type(None)):
            return
        try:
            num_extra_points = array.shape[0] - self.points.shape[0]
            self.points = array
            if num_extra_points > 0:
                for key in self.point_data.keys():
                    shape = self.point_data[key].shape
                    dtype = self.point_data[key].dtype

                    # vectors
                    if len(shape) > 1:
                        append_shape = (num_extra_points, shape[1])
                        self.point_data[key] = np.vstack(
                            [self.point_data[key], np.empty(append_shape, dtype) * np.nan]
                        )
                    # scalars
                    else:
                        append_shape = (num_extra_points,)
                        self.point_data[key] = np.append(
                            self.point_data[key], np.empty(append_shape, dtype) * np.nan
                        )

            elif num_extra_points < 0:
                raise NotImplementedError(
                    "Assigning less nodes than the original not implemented yet."
                )

        except:
            LOGGER.warning("Failed to set nodes.")
            return

    @property
    def tetrahedrons(self):
        """Tetrahedrons num_tetra x 4."""
        return self.cells_dict[pv.CellType.TETRA]

    @tetrahedrons.setter
    def tetrahedrons(self, value: np.ndarray):
        # sets tetrahedrons of UnstructuredGrid
        try:
            points = self.points
            celltypes = np.full(value.shape[0], pv.CellType.TETRA, dtype=np.int8)
            tetra = np.hstack([np.full(len(celltypes), 4)[:, None], value])
            super().__init__(tetra, celltypes, points)
        except:
            LOGGER.warning("Failed to set tetrahedrons.")
            return

    def __init__(self, *args):
        super().__init__(*args)

        self.boundaries: List[SurfaceMesh] = []
        """List of boundary surface meshes within the part."""
        pass

    @property
    def part_ids(self) -> np.ndarray:
        """Array of part ids indicating to which part the tetrahedron belongs.

        Notes
        -----
        This is derived from the "part-id" field in cell data
        """
        # NOTE "tags" should be removed.
        try:
            value = self.cell_data["tags"].astype(int)
            return value
        except (KeyError, NameError):
            LOGGER.warning("'tags' field not found in self.cell_data")
            value = None
        try:
            value = self.cell_data["part-id"].astype(int)
            return value
        except (KeyError, NameError):
            LOGGER.warning("'part-id' field not found in self.cell_data")
            value = None
        return value

    @property
    def boundary_names(self) -> List[str]:
        """Iterate over boundaries and returns their names."""
        return [b.name for b in self.boundaries]

    def _sync_nodes_of_surfaces(self):
        """Synchronize the node array of each associated surface.

        Notes
        -----
        Temporary until this module is refactored.
        """
        for b in self.boundaries:
            b.nodes = self.nodes

        return

    def get_surface_from_name(self, name: str = None):
        """Return a list of surfaces that match the given list of names.

        Notes
        -----
        Returns single surface. When multiple matches are found returns list of surfaces
        """
        surfaces_search = self.boundaries
        surfaces = [s for s in surfaces_search if s.name == name]
        if len(surfaces) == 0:
            return None
        if len(surfaces) == 1:
            return surfaces[0]
        else:
            return surfaces


class PartType(Enum):
    """Stores valid part types."""

    VENTRICLE = "ventricle"
    ATRIUM = "atrium"
    SEPTUM = "septum"
    ARTERY = "artery"
    MYOCARDIUM = "myocardium"
    UNDEFINED = "undefined"


class Part:
    """Part class."""

    @property
    def surfaces(self) -> List[SurfaceMesh]:
        """List of surfaces belonging to part."""
        surfaces = []
        for key, value in self.__dict__.items():
            if isinstance(value, SurfaceMesh):
                surfaces.append(value)
        return surfaces

    @property
    def surface_names(self) -> List[str]:
        """List of surface names belonging to part."""
        surface_names = []
        for key, value in self.__dict__.items():
            if isinstance(value, SurfaceMesh):
                surface_names.append(value.name)
        return surface_names

    def get_point(self, pointname: str) -> Point:
        """Get point from part."""
        for point in self.points:
            if point.name == pointname:
                return point
        LOGGER.error("Cannot find point {0:s}.".format(pointname))
        return None

    def __init__(self, name: str = None, part_type: PartType = PartType.UNDEFINED) -> None:
        self.name = name
        """Name of the part."""
        self.pid = None
        """Part ID."""
        self.mid = None
        """Material id associated with part."""
        self.part_type: PartType = part_type
        """Type of the part."""
        self.element_ids: np.ndarray = np.empty((0, 4), dtype=int)
        """Array holding element ids that make up this part."""
        self.points: List[Point] = []
        """Points of interest belonging to the part."""
        self.caps: List[Cap] = []
        """List of caps belonging to the part."""
        self.cavity: Cavity = None

        self.fiber: bool = False
        """If this part has fiber/sheet data."""
        self.active: bool = False
        """If active stress will be established."""

        self.meca_material: MechanicalMaterialModel = MechanicalMaterialModel.DummyMaterial()
        """Material model will be assiggned in Simulator."""

        self.ep_material: EPMaterialModel = EPMaterialModel.DummyMaterial()
        """EP Material model will be assiggned in Simulator."""

        """Cavity belonging to the part."""
        if self.part_type in [PartType.VENTRICLE]:
            self.apex_points: List[Point] = []
            """Points on apex."""

        self._add_surfaces()

    def _add_surfaces(self):
        """Add surfaces to the part."""
        if self.part_type in [PartType.VENTRICLE, PartType.ATRIUM]:
            self.endocardium = SurfaceMesh("{0} endocardium".format(self.name))
            """Endocardium."""
            self.epicardium = SurfaceMesh("{0} epicardium".format(self.name))
            """Epicardium."""
            if self.part_type == PartType.VENTRICLE:
                self.septum = SurfaceMesh("{0} septum".format(self.name))
                """Septum surface."""
        elif self.part_type in [PartType.ARTERY]:
            self.wall = SurfaceMesh("{0} wall".format(self.name))
            """Wall."""
        return

    def _add_myocardium_part(self):
        self.myocardium = Part(name="myocardium", part_type=PartType.MYOCARDIUM)
        return

    def _add_septum_part(self):
        self.septum = Part(name="septum", part_type=PartType.SEPTUM)
        return
