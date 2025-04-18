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

from __future__ import annotations

from enum import Enum
from typing import Literal

import networkx as nx
import numpy as np
import pyvista as pv

from ansys.health.heart.objects import Mesh, SurfaceMesh
from ansys.health.heart.settings.material.ep_material import EPMaterial


class ConductionPathType(Enum):
    """Conduction beam types."""

    LEFT_PURKINJE = "Left-purkinje"
    """Left Purkinje network."""
    RIGHT_PURKINJE = "Right-purkinje"
    """Right Purkinje network."""
    SAN_AVN = "SAN_to_AVN"
    """Sino-atrial node to atrio-ventricular node."""
    MID_SAN_AVN = "MID_SAN_to_AVN"
    """Sino-atrial node to atrio-ventricular node."""
    POST_SAN_AVN = "POST_SAN_to_AVN"
    """Sino-atrial node to atrio-ventricular node."""

    LEFT_BUNDLE_BRANCH = "Left bundle branch"
    """Left bundle branch."""
    RIGHT_BUNDLE_BRANCH = "Right bundle branch"
    """Right bundle branch."""
    HIS_TOP = "His_top"
    """Top part of the His bundle."""
    HIS_LEFT = "His_left"
    """Left part of the His bundle."""
    HIS_RIGHT = "His_right"
    """Right part of the His bundle."""
    BACHMANN_BUNDLE = "Bachmann bundle"
    """Bachmann bundle."""


connections = {
    ConductionPathType.LEFT_PURKINJE: [ConductionPathType.HIS_LEFT, None],
    ConductionPathType.RIGHT_PURKINJE: [ConductionPathType.HIS_RIGHT, None],
    ConductionPathType.SAN_AVN: [None, ConductionPathType.HIS_TOP],
    ConductionPathType.MID_SAN_AVN: [ConductionPathType.SAN_AVN, ConductionPathType.HIS_TOP],
    ConductionPathType.POST_SAN_AVN: [ConductionPathType.SAN_AVN, ConductionPathType.HIS_TOP],
    ConductionPathType.HIS_TOP: [
        ConductionPathType.SAN_AVN,
        None,  # [ConductionBeamType.HIS_LEFT, ConductionBeamType.HIS_RIGHT],
    ],
    ConductionPathType.HIS_LEFT: [
        ConductionPathType.HIS_TOP,
        ConductionPathType.LEFT_BUNDLE_BRANCH,
    ],
    ConductionPathType.HIS_RIGHT: [
        ConductionPathType.HIS_TOP,
        ConductionPathType.RIGHT_BUNDLE_BRANCH,
    ],
    ConductionPathType.LEFT_BUNDLE_BRANCH: [
        ConductionPathType.HIS_LEFT,
        ConductionPathType.LEFT_PURKINJE,
    ],
    ConductionPathType.RIGHT_BUNDLE_BRANCH: [
        ConductionPathType.HIS_RIGHT,
        ConductionPathType.RIGHT_PURKINJE,
    ],
    ConductionPathType.BACHMANN_BUNDLE: [ConductionPathType.SAN_AVN, None],
}


class ConductionPath:
    """ConductionPath class."""

    def __init__(
        self,
        name: ConductionPathType,
        mesh: Mesh,
        id: int,
        is_connected: np.ndarray,
        relying_surface: pv.PolyData,
        material: EPMaterial = EPMaterial.DummyMaterial(),
    ):
        """Create a conduction path.

        Parameters
        ----------
        name : ConductionPathType
            Name of the conduction path.
        mesh : Mesh
            Line mesh pf the path.
        id : int
            ID of the conduction path.
        is_connected : np.ndarray
            Mask array of points connected to solid mesh.
        relying_surface : pv.PolyData
            Surface mesh where the conduction path is relying on.
        material : EPMaterial, default: EPMaterial.DummyMaterial()
            EP Material property.
        """
        self.name = name
        self.mesh = mesh.copy()
        self.id = id
        self.is_connected = is_connected

        # TODO: check if the mesh lays on the relying_surface
        self.relying_surface = relying_surface
        self.ep_material = material

        self._assign_data()
        self._assign_connection()

    def _assign_connection(self):
        self._up = connections[self.name][0]
        self._down = connections[self.name][1]

    def _assign_data(self):
        # save data into mesh
        self.mesh.point_data["_is-connected"] = self.is_connected
        self.mesh.cell_data["_line-id"] = self.id * np.ones(self.mesh.n_cells)

    def plot(self):
        """Plot the conduction path with underlying surface."""
        plotter = pv.Plotter()
        plotter.add_mesh(self.relying_surface, color="w", opacity=0.5)
        plotter.add_mesh(self.mesh, line_width=2)
        plotter.show()

    @property
    def length(self):
        """Length of the conduction path."""
        return self.mesh.length

    @staticmethod
    def create_from_keypoints(
        name: ConductionPathType,
        keypoints: list[np.ndarray],
        id: int,
        base_mesh: pv.PolyData | pv.UnstructuredGrid,
        connection: Literal["none", "first", "last", "all"] = "none",
        line_length: float | None = 1.5,
    ) -> ConductionPath:
        """Create a conduction path on a base mesh through a set of keypoints.

        Parameters
        ----------
        name : ConductionPathType
            Name of the conduction path.
        keypoints : list[np.ndarray]
            Keypoints used to construct the path on the base mesh.
        id : int
            ID of the conduction path.
        base_mesh : pv.PolyData | pv.UnstructuredGrid
            Base mesh where the conductionn path is created. If PolyData, then the
            result is a geodesic path on the surface. If an pv.UnstructuredGrid, then the
            result the shortest path in the solid.
        connection : Literal[&quot;none&quot;, &quot;first&quot;, &quot;last&quot;, &quot;all&quot;]
        , default: "none"
            Describes how the path is connected to the solid mesh.
        line_length : float | None, default: 1.5
            Length of line element in case of refinement.

        Returns
        -------
        ConductionPath
            ConductionPath.
        """
        if isinstance(base_mesh, pv.PolyData):
            under_surface = base_mesh
            path_mesh = _create_path_on_surface(keypoints, under_surface, line_length)
        else:
            path_mesh, under_surface = _create_path_in_solid(keypoints, base_mesh, line_length)

        is_connceted = np.zeros(path_mesh.n_points)
        if connection == "first":
            is_connceted[0] = 1
        elif connection == "last":
            is_connceted[-1] = 1
        elif connection == "all":
            # only connect nodes located on the basemesh (if there is an refinement)
            is_connceted[path_mesh.point_data["base_mesh_nodes"]] = 1

        return ConductionPath(name, path_mesh, id, is_connceted, under_surface)

    @staticmethod
    def create_from_k_file(
        name: ConductionPathType, k_file: str, id: int, base_mesh: pv.PolyData, model
    ) -> ConductionPath:
        """Build conduction path from LS-DYNA k-file.

        Parameters
        ----------
        name : ConductionPathType
            Conduction path name.
        k_file : str
            Path to LS-DYNA k-file.
        id : int
            ID of the conduction path.
        base_mesh : pv.PolyData
            Surface mesh where the conduction path is relying on.
        model : HeartModel
            HeartModel object.

        Returns
        -------
        ConductionPath
            Conduction path object.
        """
        beam_nodes, edges, mask, _ = _read_purkinje_kfile(k_file)

        # build tree: beam_nodes and solid_points
        original_points_order = np.unique(edges[np.invert(mask)])
        """TODO: remove model from input, read points from k files."""
        solid_points = model.mesh.points[original_points_order]
        connectivity = np.empty_like(edges)
        np.copyto(connectivity, edges)

        # create ids of solid points and fill connectivity
        _, _, inverse_indices = np.unique(
            connectivity[np.logical_not(mask)], return_index=True, return_inverse=True
        )
        connectivity[np.logical_not(mask)] = inverse_indices + max(connectivity[mask]) + 1
        celltypes = np.full((connectivity.shape[0], 1), 2)
        connectivity = np.hstack((celltypes, connectivity))
        beam_points = np.vstack([beam_nodes, solid_points])
        # NOTE LS-DYNA create a new apex node as origin
        # TODO: merge it to apex of solid mesh
        is_connected = np.concatenate(
            [np.zeros(len(beam_nodes)), np.ones(len(solid_points))]
        ).astype(np.int64)

        beam_net = pv.PolyData(beam_points, lines=connectivity)

        return ConductionPath(name, beam_net, id, is_connected, base_mesh)


def _fill_points(point_start: np.array, point_end: np.array, length: float) -> np.ndarray:
    """Create points in a line defined by a start and an end point.

    Parameters
    ----------
    point_start : np.array
        Start point.
    point_end : np.array
        End point.
    length : float
        Length.

    Returns
    -------
    np.ndarray
        List of created points.
    """
    line_vector = point_end - point_start
    line_length = np.linalg.norm(line_vector)
    n_points = int(np.round(line_length / length)) + 1
    points = np.zeros([n_points, 3])
    points = np.linspace(point_start, point_end, n_points)
    return points


def _refine_points(nodes: np.array, length: float = None) -> tuple[np.ndarray, np.ndarray]:
    """Add new points between two points.

    Parameters
    ----------
    nodes : np.array
        Nodes to be refined.
    length : float, default None
        Length of the line element.
        If None, no refinement is done.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Refined nodes and mask of original nodes.
    """
    if length is None:  # No refinement
        return nodes, np.ones(len(nodes), dtype=bool)

    org_node_id = []
    refined_nodes = [nodes[0, :]]
    org_node_id.append(0)

    for i_cell in range(len(nodes) - 1):
        point_start = nodes[i_cell, :]
        point_end = nodes[i_cell + 1, :]
        points = _fill_points(point_start, point_end, length=length)

        refined_nodes = np.vstack((refined_nodes, points[1:, :]))
        org_node_id.append(len(refined_nodes) - 1)

    # set to True if it's an original node
    mask = np.zeros(len(refined_nodes), dtype=bool)
    mask[org_node_id] = True
    return refined_nodes, mask


def _create_path_on_surface(
    key_points: list[np.ndarray], surface: pv.PolyData, line_length: float
) -> pv.PolyData:
    """Create a geodesic path between key points.

    Parameters
    ----------
    key_points : list[np.ndarray]
        Points to be connected by the geodesic path.
    surface : pv.PolyData
        Surface on which the path is created.
    refine_length : float
        Length of the line element.

    Returns
    -------
    pv.PolyData
        Lines created by the geodesic path.
    """
    path_points = []
    for i in range(len(key_points) - 1):
        p1 = key_points[i]
        p2 = key_points[i + 1]

        path = surface.geodesic(surface.find_closest_point(p1), surface.find_closest_point(p2))
        for point in path.points:
            path_points.append(point)

    path_points, mask = _refine_points(np.array(path_points), length=line_length)

    path = pv.lines_from_points(path_points)
    path.point_data["base_mesh_nodes"] = mask
    return path


def _create_path_in_solid(
    key_points: list[np.ndarray], volume: pv.UnstructuredGrid, line_length: float
) -> tuple[pv.PolyData, pv.PolyData]:
    """Create a path in the solid mesh.

    Parameters
    ----------
    key_points : list[np.ndarray]
        Key points to be connected by the path, 2 points are required.
    volume : pv.UnstructuredGrid
        Solid mesh where the path is created.
    line_length : float
        Length of the line element.

    Returns
    -------
    tuple[pv.PolyData, pv.PolyData]
        Path mesh and surface mesh where the path is created.
    """
    if len(key_points) != 2:
        TypeError("Can only define 2 keypoints.")
        return

    # keep only tetra cells
    mesh = volume.extract_cells_by_type(pv.CellType.TETRA)

    # do the search in a small region for efficiency
    start = key_points[0]
    end = key_points[1]
    center = 0.5 * (start + end)
    radius = 3 * np.linalg.norm(start - center)
    sphere = pv.Sphere(center=center, radius=radius)

    # extract region
    cell_center = mesh.cell_centers()
    ids = np.where(cell_center.select_enclosed_points(sphere)["SelectedPoints"])[0]
    sub_mesh = mesh.extract_cells(ids)

    # search shortes path across cells
    source_id = sub_mesh.find_closest_point(start)
    target_id = sub_mesh.find_closest_point(end)
    graph = _mesh_to_nx_graph(sub_mesh)

    # ids are in submesh
    ids = nx.shortest_path(graph, source=source_id, target=target_id)
    coords = sub_mesh.points[ids]

    #
    path_points, mask = _refine_points(coords, length=line_length)
    path = pv.lines_from_points(path_points)
    path.point_data["base_mesh_nodes"] = mask

    # seg
    # TODO: split function
    tetras = sub_mesh.cells.reshape(-1, 5)[:, 1:]
    triangles = np.vstack(
        (
            tetras[:, [0, 1, 2]],
            tetras[:, [0, 1, 3]],
            tetras[:, [0, 2, 3]],
            tetras[:, [1, 2, 3]],
        )
    )  # TODO: replace by pv extract_surface()
    segment = []
    for i, j in zip(ids[0:-1], ids[1:]):
        for tri in triangles:
            if i in tri and j in tri:
                segment.append(tri)
                break
    segment = np.array(segment)

    # segment2 = sub_mesh["_global-point-ids"][segment]

    surf = SurfaceMesh(
        name="his_bundle_segment",  # NOTE
        triangles=segment,
        nodes=sub_mesh.points,
    )
    return path, surf


def _mesh_to_nx_graph(mesh: pv.UnstructuredGrid) -> nx.Graph:
    """Convert tetra mesh to graph."""
    graph = nx.Graph()
    # Add nodes
    for i, point in enumerate(mesh.points):
        graph.add_node(i, pos=tuple(point))

    # Assume all cells are tetra
    cells = np.array(mesh.cells).reshape(-1, 5)[:, 1:]
    # Add edges
    for cell in cells:
        graph.add_edge(cell[0], cell[1])
        graph.add_edge(cell[1], cell[2])
        graph.add_edge(cell[2], cell[0])
        graph.add_edge(cell[0], cell[3])
        graph.add_edge(cell[1], cell[3])
        graph.add_edge(cell[2], cell[3])

    return graph


def _read_purkinje_kfile(filename: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Read purkinje k file.

    It contains new created nodes to create Purkinje network
    and all the beam elements of Purkinje network.

    Parameters
    ----------
    filename : str
        Filename of the LS-DYNA keyword file that contains the Purkinje network.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        Coordinates of new created nodes.

        Connectivity of the beam elements.
            If mask is True, ID is new created node.
            If mask is False, ID is original node.

        Mask of connectivity.
            True for new created nodes and False for original nodes.

        Part ID of the beam elements.
    """
    # Open file and import beams and created nodes
    with open(filename, "r") as file:
        start_nodes = 0
        lines = file.readlines()

    # find line ids delimiting node data and edge data
    start_nodes = np.array(np.where(["*NODE" in line for line in lines]))[0][0]
    end_nodes = np.array(np.where(["*" in line for line in lines]))
    end_nodes = end_nodes[end_nodes > start_nodes][0]
    start_beams = np.array(np.where(["*ELEMENT_BEAM" in line for line in lines]))[0][0]
    end_beams = np.array(np.where(["*" in line for line in lines]))
    end_beams = end_beams[end_beams > start_beams][0]

    # load node data
    node_data = np.loadtxt(filename, skiprows=start_nodes + 1, max_rows=end_nodes - start_nodes - 1)
    node_ids = node_data[:, 0].astype(int) - 1  # 0 based
    coords = node_data[:, 1:4]

    # load beam data
    beam_data = np.loadtxt(
        filename, skiprows=start_beams + 1, max_rows=end_beams - start_beams - 1, dtype=int
    )
    edges = beam_data[:, 2:4] - 1  # 0 based
    pid = beam_data[:, 1]

    edges_mask = np.isin(edges, node_ids)  # True for new created nodes
    edges[edges_mask] -= node_ids[0]  # beam nodes id start from 0

    return coords, edges, edges_mask, pid
