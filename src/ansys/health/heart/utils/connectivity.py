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

"""Module containing methods for mesh connectivity."""

import copy
from typing import Tuple

import numpy as np

from ansys.health.heart import LOG as LOGGER


def get_faces_tetra(tetra: np.ndarray) -> np.ndarray:
    """Get faces that make up the tetrahedrons."""
    num_tetra = tetra.shape[0]
    faces = np.zeros((num_tetra, 3, 4), dtype=int)
    masks = np.array(
        [
            [True, True, True, False],
            [True, True, False, True],
            [True, False, True, True],
            [False, True, True, True],
        ]
    )
    for ii, mask in enumerate(masks):
        faces[:, :, ii] = tetra[:, mask]

    return faces


def face_tetra_connectivity(tetra: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute the tetra-face connectivity tables."""
    import time as time

    LOGGER.debug("Establishing tetra-face connectivity...")
    t0 = time.time()

    faces = get_faces_tetra(tetra)
    # NOTE: reshape faces such that shape is (NumTetra*4, 3) with following structure:
    # i n1 n2 n3
    # 0 tet1_face1
    # 1 tet1_face2
    # 2 tet1_face3
    # 3 tet1_face4
    # 4 tet2_face1
    # 5 tet2_face2
    # 6 tet2_face3
    # 7 tet2_face4
    # ...
    num_tetra = tetra.shape[0]
    faces_1 = np.reshape(faces.transpose(0, 2, 1), (4 * num_tetra, 3))

    # sort faces to find duplicates
    faces_sorted = np.sort(faces_1, axis=1)
    np.sort(faces_sorted, axis=0)

    # find duplicate rows
    faces_unique, index, inverse, counts = np.unique(
        faces_sorted, return_index=True, return_inverse=True, return_counts=True, axis=0
    )

    # find duplicate rows in reversed order
    faces_sorted_flip = np.flipud(faces_sorted)
    faces_unique_r, index_r, inverse_r, counts_r = np.unique(
        faces_sorted_flip, return_index=True, return_inverse=True, return_counts=True, axis=0
    )

    tetra_ids = np.repeat(np.arange(0, num_tetra, 1), 4)
    tetra_ids_flip = np.flipud(tetra_ids)

    # get connected tetra ID for each face (two for interior and one for boundary face)
    c0 = tetra_ids[index][inverse]
    c1 = np.flip(tetra_ids_flip[index_r][inverse_r])

    # remove any duplicate faces
    mapper = np.sort(index)
    faces_1 = faces_1[mapper, :]
    c0 = c0[mapper]
    c1 = c1[mapper]

    t1 = time.time()
    LOGGER.debug("Time elapsed: {:.1f} s".format(t1 - t0))

    return faces_1, c0, c1


def get_face_type(faces: np.ndarray, face_cell_connectivity: np.ndarray) -> np.ndarray:
    """Establish the face type, which indicates whether it is a boundary or an interior face.

    Parameters
    ----------
    faces : np.ndarray
        Array with face definitions.
    face_cell_connectivity : np.ndarray
        Array describing the cells that each of the faces is connected to.
        For example, ``np.array([c0, c1])``.

    Returns
    -------
    np.ndarray
        Type of face, which is either interior ``(face_type = 1)``
        or boundary ``(face_type = 2)``.
    """
    interior_face_ids = face_cell_connectivity[:, 0] != face_cell_connectivity[:, 1]
    boundary_face_ids = face_cell_connectivity[:, 0] == face_cell_connectivity[:, 1]
    face_types = np.zeros((faces.shape[0]), dtype=int)
    face_types[interior_face_ids] = 1
    face_types[boundary_face_ids] = 2
    num_assigned = np.sum(boundary_face_ids) + np.sum(interior_face_ids)
    if num_assigned != faces.shape[0]:
        un_assigned = faces.shape[0] - num_assigned
        raise Exception(
            f"Not all faces are assigned a type. {un_assigned}/{faces.shape[0]} are not assigned."
        )
    return face_types


def get_edges_from_triangles(triangles: np.ndarray) -> np.ndarray:
    """Generate an array of edges from an array of triangles."""
    num_triangles = triangles.shape[0]
    num_edges = num_triangles * 3
    edges = np.repeat(triangles, 3, axis=0)
    mask = np.tile(
        np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]], dtype=bool),
        (num_triangles, 1),
    )
    edges = np.reshape(edges[mask], (num_edges, 2))

    return edges


def get_free_edges(
    triangles: np.ndarray, return_free_triangles: bool = False
) -> np.ndarray | Tuple[np.ndarray, np.ndarray]:
    """Get the boundary edges that are only referenced once.

    Parameters
    ----------
    triangles : np.ndarray
        Array of triangles.
    return_free_triangles : bool, default: False
        Whether to return the free triangles.

    Returns
    -------
    free_edges : np.ndarray
        NumPy array with the free edges.
    free_triangles: np.ndarray, optional
        Numpy array with the triangles that use these free edges
    """
    edges = get_edges_from_triangles(triangles)

    edges_sort = np.sort(edges, axis=1)

    unique_edges, idx, counts = np.unique(edges_sort, axis=0, return_counts=True, return_index=True)
    free_edges = edges[idx, :][counts == 1, :]

    if not return_free_triangles:
        return free_edges

    elif return_free_triangles:
        # get free triangles
        free_triangles = triangles[
            np.argwhere(np.sum(np.isin(triangles, free_edges), axis=1) == 2).flatten(), :
        ]

        return free_edges, free_triangles


#! NOTE: method not used anymore, but may still be useful.
def edge_connectivity(
    edges: np.ndarray, return_type: bool = False, sort_closed: bool = False
) -> np.ndarray:
    """Group edges by connectivity.

    Parameters
    ----------
    edges : np.array
        NumEdges x 2 NumPy arrays with edge definitions.
    return_type : bool, default: False
        Whether to return the edge group type. If ``True``, the function
        returns a list of strings with these types:

        - ``open``: Edge group is open-ended.
        - ``closed``: Edge group forms a closed edge loop.

    sort_closed : bool, default: False
        Whether to sort closed edge loops.

    Returns
    -------
    edge_groups : np.ndarray
        Grouped edges by connectivity.
    group_types : list[str], optional
        Type of the edge group. Options are ``open`` or ``closed``.

    Notes
    -----
    This method uses an implementation of a depth-first search. For more information,
    seeDepth-first search <https://en.wikipedia.org/wiki/Depth-first_search>`_ on
    the Wikipedia site.

    Performance of this method is not tested. It might not be suitable for large arrays of edges.
    """

    def _dfs(visited, graph, node):
        # Nested in Dept-first search algorithm
        if node not in visited:
            # print(node)
            visited.add(node)
            for neighbor in graph[node]:
                _dfs(visited, graph, neighbor)

    # create adjacency list (typically referred to as "graph")
    graph = {}
    node_ids = np.unique(edges)
    for node in node_ids:
        mask = node != edges
        mask[np.all(mask, axis=1), :] = False
        connected_nodes = edges.flatten()[mask.flatten()]
        graph[node] = connected_nodes

    # check connectivity of each node using DFS
    # group connected edges
    node_ids_visited = np.zeros(node_ids.shape[0], dtype=bool)
    edge_groups = []
    while not np.all(node_ids_visited):
        # keep track of visited nodes for this group of edges
        visited = set()

        # node ID to start from (finds first un-visited node)
        start_node_id = node_ids[np.where(np.invert(node_ids_visited))[0][0]]

        # call dept first algorithm to find connectivity
        _dfs(visited, graph, start_node_id)

        # retrieve edge definitions
        edge_group = edges[np.all(np.isin(edges, list(visited)), axis=1)]
        edge_groups.append(edge_group)

        node_ids_visited[np.isin(node_ids, list(visited))] = True

    # check whether edges form a closed loop
    group_types = []
    if return_type:
        for edge_group in edge_groups:
            counts = np.unique(edge_group, return_counts=True)[1]
            if np.all(counts == 2):
                # print("Closed edge loop")
                group_types.append("closed")
            elif np.any(counts != 2):
                group_types.append("open")
                # print("Open edge tree")

    # sort any closed edge loops
    if sort_closed:
        for ii, edge_group in enumerate(edge_groups):
            if group_types[ii] == "closed":
                edges = edge_group.tolist()
                remaining_edges = edges
                next_edge = edges[0]
                sorted_edges = [next_edge]
                remaining_edges.pop(0)
                while len(remaining_edges) > 0:
                    # find connected edge of last edge
                    node = sorted_edges[-1][1]
                    mask = np.array(edges) == node
                    if np.sum(mask[:, 1]) == 1:
                        flip = True
                    elif np.sum(mask[:, 0]) == 1:
                        flip = False
                    else:
                        raise ValueError("Expecting just one match")
                    idx = np.where(np.any(mask, axis=1))[0][0]
                    if flip:
                        sorted_edges.append(np.flip(remaining_edges[idx]).tolist())
                    else:
                        sorted_edges.append(remaining_edges[idx])
                    remaining_edges.pop(idx)
                edge_groups[ii] = np.array(sorted_edges)

    if return_type:
        return edge_groups, group_types
    else:
        return edge_groups


def remove_triangle_layers_from_trimesh(triangles: np.ndarray, iters: int = 1) -> np.ndarray:
    """Remove boundary triangles.

    Parameters
    ----------
    triangles : np.ndarray
        Array of triangles.
    iters : int, default: 1
        Number of iterations.

    Returns
    -------
    np.ndarray
        Reduced set of triangles.
    """
    reduced_triangles = copy.deepcopy(triangles)
    for ii in range(0, iters, 1):
        num_triangles = reduced_triangles.shape[0]
        edges = get_edges_from_triangles(reduced_triangles)
        free_edges = get_free_edges(reduced_triangles)

        # find elements connected to the free edges
        edges = np.reshape(edges, (3, 2, num_triangles))
        free_nodes = np.unique(free_edges)

        idx_triangles_boundary = np.any(np.isin(reduced_triangles, free_nodes), axis=1)

        LOGGER.debug("Removing {0} connected triangles...".format(np.sum(idx_triangles_boundary)))

        # remove boundary triangles
        reduced_triangles = reduced_triangles[~idx_triangles_boundary, :]

    return reduced_triangles
