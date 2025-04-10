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

"""Conduction system class."""

from enum import Enum
from typing import Literal

import numpy as np
import pyvista as pv

from ansys.health.heart.objects import Mesh
from ansys.health.heart.settings.material.ep_material import EPMaterial


class ConductionBeamType(Enum):
    """Conduction beam types."""

    LEFT_PURKINJE = "Left-purkinje"
    """Left Purkinje network."""
    RIGHT_PURKINJE = "Right-purkinje"
    """Right Purkinje network."""
    SAN_AVN = "SAN_to_AVN"
    """Sino-atrial node to atrio-ventricular node."""
    LEFT_BUNDLE_BRANCH = "Left bundle branch"
    """Left bundle branch."""
    RIGHT_BUNDLE_BRANCH = "Right bundle branch"
    """Right bundle branch."""
    HIS_TOP = "His_top"
    """His bundle top part."""
    HIS_LEFT = "His_left"
    HIS_RIGHT = "His_right"
    BACHMANN_BUNDLE = "Bachman bundle"
    """Bachmann bundle."""


class ConductionBeams:
    """Conduction beams class."""

    def __init__(
        self,
        name: ConductionBeamType,
        mesh: Mesh,
        id: int,
        is_connceted: np.ndarray,
        relying_surface: pv.PolyData,
        material: EPMaterial = EPMaterial.DummyMaterial(),
    ):
        """Create a conduction beam.

        Parameters
        ----------
        name : ConductionBeamType
            name of the conduction beam.
        mesh : Mesh
            Beam mesh.
        id : int
            id of the conduction beam.
        is_connceted : np.ndarray
            mask array of points connected to solid mesh.
        relying_surface : pv.PolyData
            surface mesh where the conduction beam is relying on.
        material : EPMaterial, optional
            EP material property, by default EPMaterial.DummyMaterial()
        """
        self.name = name
        self.mesh = mesh.copy()
        self.id = id
        self.is_connected = is_connceted

        # TODO: check if mesh are on relying_surface
        self.relying_surface = relying_surface
        self.ep_material = material

        self._assign_data()

    def _assign_data(self):
        # tempo script to support old structure
        self.mesh.point_data["_is-connected"] = self.is_connected
        self.mesh.cell_data["_line-id"] = self.id * np.ones(self.mesh.n_cells)

    def plot(self):
        """Plot the conduction beam."""
        plotter = pv.Plotter()
        plotter.add_mesh(self.relying_surface, color="w", opacity=0.5)
        # self.conduction_system.set_active_scalars("_line-id")
        # beams = self.conduction_system
        plotter.add_mesh(self.mesh, line_width=2)
        plotter.show()

    @property
    def length(self):
        """Length of the conduction beam."""
        return self.mesh.length

    @staticmethod
    def create_from_keypoints(
        name: ConductionBeamType,
        keypoints: list[np.ndarray],
        id: int,
        base_mesh: pv.PolyData | pv.UnstructuredGrid,
        connection: Literal["none", "first", "last", "all"] = "none",
        refine_length: float = 1.5,
    ):
        """Create a conduction beam by providing keypoints.

        Parameters
        ----------
        name : ConductionBeamType
             name of the conduction beam.
        keypoints : list[np.ndarray]
            keypoints are used to create a path on the base mesh.
        id : int
            id of the conduction beam.
        base_mesh : pv.PolyData | pv.UnstructuredGrid
            base mesh where the conduction beam is created.
            If PolyData, results are geodesic line on the surface.
            If UnstructuredGrid, results are lines in the solid.
        connection : Literal[&quot;none&quot;, &quot;first&quot;, &quot;last&quot;, &quot;all&quot;]
        , optional
            describe how beam is connected to solid mesh, by default "none"
        refine_length : float, optional
            beam length, by default 1.5

        Returns
        -------
        _type_
            _description_
        """
        if isinstance(base_mesh, pv.PolyData):
            under_surface = base_mesh
            mesh = _create_path_by_geodesic(keypoints, under_surface, refine_length)
        else:
            NotImplementedError

        is_connceted = np.zeros(mesh.n_points)
        if connection == "first":
            is_connceted[0] = 1
        elif connection == "last":
            is_connceted[-1] = 1
        elif connection == "all":
            is_connceted[:] = 1

        return ConductionBeams(name, mesh, id, is_connceted, under_surface)


def _create_line(point_start: np.array, point_end: np.array, beam_length: float) -> np.ndarray:
    """Create points in a line defined by a start point and an end point.

    Parameters
    ----------
    point_start : np.array
        Start Point.
    point_end : np.array
        End point.
    beam_length : float
        Beam length.

    Returns
    -------
    np.ndarray:
        List of created points.
    """
    line_vector = point_end - point_start
    line_length = np.linalg.norm(line_vector)
    n_points = int(np.round(line_length / beam_length)) + 1
    points = np.zeros([n_points, 3])
    # beams = np.zeros([n_points - 1, 2])
    points = np.linspace(point_start, point_end, n_points)
    # beams[:, 0] = np.linspace(0, n_points - 2, n_points - 1, dtype=int)
    # beams[:, 1] = np.linspace(0, n_points - 2, n_points - 1, dtype=int) + 1
    return points


def _refine_line(nodes: np.array, beam_length: float) -> np.ndarray:
    new_nodes = [nodes[0, :]]
    for beam_id in range(len(nodes) - 1):
        point_start = nodes[beam_id, :]
        point_end = nodes[beam_id + 1, :]
        points = _create_line(point_start, point_end, beam_length=beam_length)
        new_nodes = np.vstack((new_nodes, points[1:, :]))
    return new_nodes


def _create_path_by_geodesic(
    key_points: list[np.ndarray], surface: pv.PolyData, refine_length: float
) -> pv.PolyData:
    """Create a path by geodesic line between key points.

    Parameters
    ----------
    key_points : list[np.ndarray]
        points to be connected by geodesic lines.
    surface : pv.PolyData
        surface where the path is created.
    refine_length : float
        refine beam length.

    Returns
    -------
    pv.PolyData
        lines created by geodesic lines.
    """
    path_points = []
    for i in range(len(key_points) - 1):
        p1 = key_points[i]
        p2 = key_points[i + 1]

        path = surface.geodesic(surface.find_closest_point(p1), surface.find_closest_point(p2))
        for point in path.points:
            path_points.append(point)

    path_points = _refine_line(np.array(path_points), beam_length=refine_length)

    return pv.lines_from_points(path_points)
