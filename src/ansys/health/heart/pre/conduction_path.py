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

from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.objects import Mesh, SurfaceMesh
from ansys.health.heart.settings.material.ep_material import EPMaterial


class ConductionPathType(Enum):
    """Conduction Path types."""

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
    LEFT_ANTERIOR_FASCILE = "Left anterior fascicle"
    """Left anterior fascicle."""
    LEFT_POSTERIOR_FASCICLE = "Left posterior fascicle"
    """Left posterior fascicle."""
    USER_PAHT_1 = "User path 1"
    """User path 1."""
    USER_PAHT_2 = "User path 2"
    """User path 2."""
    USER_PAHT_3 = "User path 3"
    """User path 3."""


class ConductionPath:
    """Conduction path class."""

    def __init__(
        self,
        name: ConductionPathType,
        mesh: Mesh,
        id: int,
        is_connected: np.ndarray,
        relying_surface: pv.PolyData,
        material: EPMaterial = EPMaterial.DummyMaterial(),
        up_path: ConductionPath | None = None,
        down_path: ConductionPath | None = None,
    ):
        """
        Initialize a conduction path.

        Parameters
        ----------
        name : ConductionPathType
            Name of the conduction path.
        mesh : Mesh
            Line mesh of the path.
        id : int
            ID of the conduction path.
        is_connected : np.ndarray
            Mask array of points connected to the solid mesh.
        relying_surface : pv.PolyData
            Surface mesh that the conduction path relies on.
        material : EPMaterial, default: EPMaterial.DummyMaterial()
            EP Material property.
        up_path : ConductionPath | None, default: None
            Upstream conduction path. Its closest point connects to the first point of this path.
        down_path : ConductionPath | None, default: None
            Downstream conduction path. Its closest point connects to the last point of this path.

        Notes
        -----
        up_path and down_path can be parallel paths, such as the 3 SA-AV paths.
        """
        self.name = name
        self.mesh = mesh.copy()
        self.id = id
        self.is_connected = is_connected
        self.relying_surface = relying_surface

        # Check that the line length is greater than 0
        if not self.mesh.compute_cell_sizes()["Length"].min() > 0:
            raise ValueError(f"{self.name} contains cells with length <= 0. ")

        # check if the mesh lays on the relying_surface
        dst = self.mesh.compute_implicit_distance(self.relying_surface)["implicit_distance"]
        LOGGER.info(
            f"Maximal distance of {self.name} to its relying surface is: {np.max(abs(dst))}."
        )

        self.ep_material = material

        self._assign_data()
        self.up_path = up_path
        self.down_path = down_path

    @property
    def up_path(self) -> ConductionPath | None:
        """Get upstream conduction path."""
        return self._up_path

    @property
    def down_path(self) -> ConductionPath | None:
        """Get downstream conduction path."""
        return self._down_path

    @up_path.setter
    def up_path(self, value: ConductionPath | None):
        """Set upstream conduction path.

        Parameters
        ----------
        value : ConductionPath | None
            Upstream conduction path, its closest point will be connected to
            the first point of this path.
        """
        if value is not None:
            origin = self.mesh.points[0]
            target_id = value.mesh.find_closest_point(origin)
            target = value.mesh.points[target_id]
            dst = np.linalg.norm(origin - target)
            LOGGER.info(f"Distance between {self.name} and {value.name} is: {dst}.")

        self._up_path = value

    @down_path.setter
    def down_path(self, value: ConductionPath | None):
        """Set downstream conduction path.

        Parameters
        ----------
        value : ConductionPath | None
            Downstream conduction path, its closest point will be connected to
            the last point of this path.
        """
        if value is not None:
            origin = self.mesh.points[-1]
            target_id = value.mesh.find_closest_point(origin)
            target = value.mesh.points[target_id]
            dst = np.linalg.norm(origin - target)
            LOGGER.info(f"Distance between {self.name} and {value.name} is: {dst}.")
        self._down_path = value

    def _assign_data(self):
        # save data into mesh
        self.mesh.point_data["_is-connected"] = self.is_connected
        self.mesh.cell_data["_line-id"] = self.id * np.ones(self.mesh.n_cells)

    def get_terminal_nodes(self) -> np.ndarray:
        """Get the terminal nodes of the conduction path.

        Notes
        -----
        The terminal nodes are the points that are referenced only once
        in the line segments.

        Returns
        -------
        np.ndarray
            Array of terminal node indices.
        """
        p1 = np.unique(self.mesh.lines.reshape(-1, 3)[:, 1])
        p2 = np.unique(self.mesh.lines.reshape(-1, 3)[:, 2])
        return np.setdiff1d(p2, p1)

    def get_terminal_coordinates(self) -> np.ndarray:
        """Get the nodal coordinates of the terminal nodes.

        Returns
        -------
        np.ndarray
            Nx3 array with the coordinates of the terminal nodes.
        """
        return self.mesh.points[self.get_terminal_nodes(), :]

    def plot(self, show_plotter: bool = True) -> pv.Plotter | None:
        """
        Plot the conduction path with its underlying surface.

        This method creates a PyVista plotter, adds the relying surface (in semi-transparent white)
        and the conduction path (as a line), and either shows the plot or returns the plotter
        for further customization.

        Parameters
        ----------
        show_plotter : bool, default: True
            If True, immediately show the plot window. If False, return the plotter
            object for further modification (e.g., adding more meshes).

        Returns
        -------
        plotter : pyvista.Plotter or None
            The PyVista plotter object if ``show_plotter`` is False, otherwise None.

        Examples
        --------
        >>> plotter = conduction_path.plot(show_plotter=False)
        >>> plotter.add_mesh(other_mesh, color="red")
        >>> plotter.show()
        """
        plotter = pv.Plotter()
        plotter.add_mesh(self.relying_surface, color="w", opacity=0.5)
        plotter.add_mesh(self.mesh, line_width=2)
        plotter.add_text(self.name.value, font_size=12, color="black", position="upper_edge")
        if show_plotter:
            plotter.show()
            return
        return plotter

    @property
    def length(self):
        """Length of the conduction path."""
        return self.mesh.length

    def add_pmj_path(
        self, pmj_list: list[int], merge_with: Literal["node", "cell"] = "cell"
    ) -> ConductionPath:
        """
        Add Purkinje-Myocardial Junction branches to the current conduction path.

        Parameters
        ----------
        pmj_list : list[int]
            Indices of points to create a Purkinje-Myocardial Junction.
        merge_with : Literal['node','cell'], default: 'cell'
            Whether to merge with a neighbor node (1 split) or cell (3 splits).

        Returns
        -------
        ConductionPath
            The updated conduction path.

        Notes
        -----
        PMJ resistance is controlled by pmjres in *EM_EP_PURKINJE_NETWORK2.
        """
        # TODO: make sure we do not create path with length of 0
        if merge_with not in ("node", "cell"):
            raise ValueError(f"merge_with must be 'node' or 'cell', got '{merge_with}'")

        new_points = []
        new_lines = []
        new_is_connected = []

        for ii in pmj_list:
            p0 = self.mesh.points[ii]
            cell_id = self.relying_surface.find_containing_cell(p0)
            neighbour_ids = self.relying_surface.get_cell(cell_id).point_ids

            for i in range(len(neighbour_ids)):
                neigh_coord = self.relying_surface.points[neighbour_ids[i]]
                new_points.append(neigh_coord)
                # The new point index will be len(points) + len(new_points) - 1
                new_idx = self.mesh.n_points + len(new_points) - 1
                new_lines.append([2, ii, new_idx])
                new_is_connected.append(1)

                if merge_with == "node":
                    break  # stop with first node

        pmj_mesh, is_connected = _path_merge(
            self.mesh, new_points, new_lines, self.is_connected, new_is_connected
        )

        # Update self in-place
        self.mesh = pmj_mesh
        self.is_connected = is_connected
        self._assign_data()

        return self

    @staticmethod
    def create_from_keypoints(
        name: ConductionPathType,
        keypoints: list[np.ndarray],
        id: int,
        base_mesh: pv.PolyData | pv.UnstructuredGrid,
        connection: Literal["first"] | None = None,
        line_length: float | None = 1.5,
        center: bool = False,
    ) -> ConductionPath:
        """
        Create a conduction path on a base mesh through a set of keypoints.

        Notes
        -----
        To add PMJ (Purkinje-Myocardial Junction) points, use the :meth:`add_pmj_path` method
        after creating the path.

        Parameters
        ----------
        name : ConductionPathType
            Name of the conduction path.
        keypoints : list[np.ndarray]
            Keypoints used to construct the path on the base mesh.
        id : int
            ID of the conduction path.
        base_mesh : pv.PolyData | pv.UnstructuredGrid
            Base mesh where the conduction path is created. If ``PolyData``, the
            result is a geodesic path on the surface. If ``pv.UnstructuredGrid``, the
            result is the shortest path in the solid.
        connection : Literal["first"] | None, default: None
            If "first", the first point of the path is marked as connected to the solid mesh.
            If None, no points are marked as connected.
        line_length : float | None, default: 1.5
            Length of the line element in case of refinement.
        center : bool, default: False
            Whether to use a geodesic path through the centers of the surface cells.
            Length of the line element in case of refinement.
        center : bool, default: False
            If True, the geodesic path passes through the centers of the surface cells.

        Returns
        -------
        ConductionPath
            The created conduction path.
        """
        # Check element types
        if isinstance(base_mesh, pv.PolyData):
            cell_types = np.unique(base_mesh.faces.reshape(-1, base_mesh.faces[0] + 1)[:, 0])
            if not np.all(cell_types == 3):  # 3 = triangle
                LOGGER.error(
                    "Base mesh contains non-triangle elements. Only triangles are supported."
                )
                return
        else:
            cell_types = np.unique(base_mesh.celltypes)
            if not np.all(cell_types == pv.CellType.TETRA):
                LOGGER.error(
                    "Base mesh contains non-tetrahedral elements. Only tetras are supported."
                )
                return

        if isinstance(base_mesh, pv.PolyData):
            under_surface = base_mesh
            if center:
                path_mesh = _create_path_on_surface_center(keypoints, under_surface, line_length)
            else:
                path_mesh = _create_path_on_surface(keypoints, under_surface, line_length)
        else:
            path_mesh, under_surface = _create_path_in_solid(keypoints, base_mesh, line_length)

        is_connected = np.zeros(path_mesh.n_points)
        if connection == "first":
            is_connected[0] = 1

        return ConductionPath(name, path_mesh, id, is_connected, under_surface)

    @staticmethod
    def create_from_k_file(
        name: ConductionPathType,
        k_file: str,
        id: int,
        base_mesh: pv.PolyData,
        model,
        merge_apex: bool = True,
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
            Surface mesh that the conduction path is relying on.
        model : HeartModel
            HeartModel object.
        merge_apex : bool, default: True
            Whether to merge the apex node with the solid mesh.

        Returns
        -------
        ConductionPath
            Conduction path.
        """
        # The method is now unnecessarily complex to build polydata of path,
        # we can just read solid + beam nodes and beam elements, then build it
        # following with clean() to remove unused nodes.

        beam_nodes, edges, mask, _ = _read_purkinje_kfile(k_file)

        # get solid points which are not in k_file
        # alternatively, can be get from reading nodes.k
        solid_points_ids = np.unique(edges[np.invert(mask)])
        solid_points = model.mesh.points[solid_points_ids]

        # create connectivity
        connectivity = np.empty_like(edges)
        np.copyto(connectivity, edges)
        _, _, inverse_indices = np.unique(
            connectivity[np.logical_not(mask)], return_index=True, return_inverse=True
        )
        connectivity[np.logical_not(mask)] = inverse_indices + max(connectivity[mask]) + 1

        # build polydata
        points = np.vstack([beam_nodes, solid_points])
        celltypes = np.full((connectivity.shape[0], 1), 2)
        connectivity = np.hstack((celltypes, connectivity))
        path = pv.PolyData(points, lines=connectivity)

        # LS-DYNA creates a new node at apex as origin of Purkinje network
        is_connected = np.concatenate(
            [np.zeros(len(beam_nodes)), np.ones(len(solid_points))]
        ).astype(np.int64)

        if merge_apex:
            is_connected[0] = 1
        return ConductionPath(name, path, id, is_connected, base_mesh)


def _path_merge(
    path: pv.PolyData,
    new_points: np.ndarray,
    new_lines: np.ndarray,
    is_connected: np.ndarray,
    new_is_connected: np.ndarray,
) -> tuple[pv.PolyData, np.ndarray]:
    """
    Merge new lines into the path and keep the last line at the end.

    Parameters
    ----------
    path : pv.PolyData
        The original path mesh.
    new_points : np.ndarray
        Array of new points to add.
    new_lines : np.ndarray
        Array of new lines to add (shape: [n_lines, 3], where each row is [2, start_idx, end_idx]).
    is_connected : np.ndarray
        Array indicating which points in the original path are connected.
    new_is_connected : np.ndarray
        Array indicating which new points are connected.

    Returns
    -------
    tuple[pv.PolyData, np.ndarray]
        The merged path mesh and the updated is_connected array.
    """
    points = path.points.copy()
    lines = path.lines.copy().reshape(-1, 3)
    orig_last_idx = len(points) - 1

    # Add new points and lines
    points = np.vstack([points, np.array(new_points)])
    # Ensure the last cell remains at the end, as required by the merging method.
    lines = np.vstack([lines[0:-1], np.array(new_lines), lines[-1]])
    is_connected = np.concatenate([is_connected, np.array(new_is_connected)])

    # Swap last point in original points with last point in merged array
    merged_last_idx = len(points) - 1
    # Swap points
    points[[orig_last_idx, merged_last_idx]] = points[[merged_last_idx, orig_last_idx]]

    # Swap is_connected
    is_connected[[orig_last_idx, merged_last_idx]] = is_connected[[merged_last_idx, orig_last_idx]]

    # Swap lines
    arr = lines[:, 1:].copy()
    arr[arr == orig_last_idx] = -1  # temp
    arr[arr == merged_last_idx] = orig_last_idx
    arr[arr == -1] = merged_last_idx

    # Rebuild PolyData
    new_path = pv.PolyData(points, lines=np.insert(arr, 0, 2, axis=1))

    return new_path, is_connected


def _fill_points(point_start: np.array, point_end: np.array, length: float) -> np.ndarray:
    """
    Create additional points in a line defined by a start and an end point.

    Parameters
    ----------
    point_start : np.array
        Start point.
    point_end : np.array
        End point.
    length : float
        Length of each segment.

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
    """
    Add new points between two points.

    Parameters
    ----------
    nodes : np.array
        Nodes to be refined.
    length : float, default None
        Length of the line element. If None, no refinement is done.

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
    key_points: list[np.ndarray], surface: pv.PolyData, line_length: float | None = None
) -> pv.PolyData:
    """
    Create a geodesic path between key points on a surface.

    Parameters
    ----------
    key_points : list[np.ndarray]
        Points to be connected by the geodesic path.
    surface : pv.PolyData
        Surface on which the path is created.
    line_length : float | None, default: None
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


def _create_path_on_surface_center(
    key_points: list[np.ndarray], surface: pv.PolyData, line_length: float | None = None
) -> pv.PolyData:
    """
    Create a geodesic path between key points through the center of surface cells.

    Parameters
    ----------
    key_points : list[np.ndarray]
        Points to be connected by the geodesic path.
    surface : pv.PolyData
        Surface on which the path is created.
    line_length : float | None, default: None
        Length of the line element.

    Returns
    -------
    pv.PolyData
        Lines created by the geodesic path.
    """
    cell_centers = surface.cell_centers().points

    graph = nx.Graph()

    for cell_id in range(surface.n_cells):
        center = cell_centers[cell_id]
        graph.add_node(cell_id, pos=center)

        neighbors = surface.cell_neighbors(cell_id)

        for neighbor_id in neighbors:
            # distance is the weight
            dist = np.linalg.norm(cell_centers[cell_id] - cell_centers[neighbor_id])
            graph.add_edge(cell_id, neighbor_id, weight=dist)

    path_points = []

    for i in range(len(key_points) - 1):
        start_cell = np.argmin(np.linalg.norm(cell_centers - key_points[i], axis=1))
        end_cell = np.argmin(np.linalg.norm(cell_centers - key_points[i + 1], axis=1))

        # Compute shortest path
        path_cells = nx.shortest_path(graph, source=start_cell, target=end_cell, weight="weight")

        # Get corresponding points
        for point in cell_centers[path_cells[:-1]]:
            path_points.append(point)

    path_points.append(cell_centers[path_cells[-1]])
    path_points = np.array(path_points)

    # replace first and end points are on surface nodes
    path_points[0] = surface.points[surface.find_closest_point(key_points[0])]
    path_points[-1] = surface.points[surface.find_closest_point(key_points[-1])]

    path_points, mask = _refine_points(path_points, length=line_length)

    path = pv.lines_from_points(path_points)
    path.point_data["base_mesh_nodes"] = mask

    return path


def _create_path_in_solid(
    key_points: list[np.ndarray], volume: pv.UnstructuredGrid, line_length: float | None = None
) -> tuple[pv.PolyData, pv.PolyData]:
    """
    Create a path in the solid mesh.

    Parameters
    ----------
    key_points : list[np.ndarray]
        Key points to be connected by the path. Requires 2 points.
    volume : pv.UnstructuredGrid
        Solid mesh where the path is created.
    line_length : float | None, default: None
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
    radius = 10 * np.linalg.norm(start - center)
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

    surf = SurfaceMesh(
        name="his_bundle_segment",  # NOTE
        triangles=segment,
        nodes=sub_mesh.points,
    )
    return path, surf


def _mesh_to_nx_graph(mesh: pv.UnstructuredGrid) -> nx.Graph:
    """
    Convert a tetrahedral mesh to a NetworkX graph.

    Parameters
    ----------
    mesh : pv.UnstructuredGrid
        The tetrahedral mesh.

    Returns
    -------
    nx.Graph
        The corresponding graph.
    """
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
    """
    Read a Purkinje k-file.

    The file contains newly created nodes for the Purkinje network
    and all the beam elements of the Purkinje network.

    Parameters
    ----------
    filename : str
        Filename of the LS-DYNA keyword file that contains the Purkinje network.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        Coordinates of new created nodes.
        Connectivity of the beam elements.
            If mask is True, ID is a new created node.
            If mask is False, ID is an original node.
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
