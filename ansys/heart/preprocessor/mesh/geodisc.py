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

"""Module for computing paths."""

import math
from typing import List, Union

import numpy as np
import vtk
from vtk.numpy_interface import dataset_adapter as dsa  # type: ignore


def get_closed_path(start_indices: Union[np.ndarray, list], surface: vtk.vtkPolyData) -> np.ndarray:
    """Get closed geodesic path on a surface from a list of start indices."""
    surf = dsa.WrapDataObject(surface)

    end_indices = np.append(start_indices[1:], start_indices[0])

    full_path = []  # node ids that form path from start indices to end indices
    for i, start_idx in enumerate(start_indices):
        end_idx = end_indices[i]

        # np.savetxt("start_end_node" + str(i) + ".csv",
        # surf.Points[ [start_idx,end_idx], : ], delimiter = ',' )

        path = vtk_geodesic(surface, start_idx, end_idx)
        full_path.extend(path[:-1])

    # np.savetxt("full_path.csv", surf.Points[ full_path, : ], delimiter = ',' )

    return np.array(full_path)


def vtk_geodesic(input: vtk.vtkPolyData, start_idx: int, end_idx: int) -> List[int]:
    """Compute the geodesic path between two vertices on a vtkPolyData surface.

    Parameters
    ----------
    input : vtk.vtkPolyData
        Input surface mesh in vtk.vtkPolyData form
    start : int
        Index of start vertex
    end : int
        Index of end vertex

    Returns
    -------
    ids : List[int]
        Array of indices which define the shortest path from start index to end index
    """
    dijkstra = vtk.vtkDijkstraGraphGeodesicPath()  # noqa
    dijkstra.SetInputData(input)
    dijkstra.SetStartVertex(start_idx)
    dijkstra.SetEndVertex(end_idx)
    dijkstra.Update()
    ids = []
    for i in range(dijkstra.GetIdList().GetNumberOfIds()):
        ids.append(dijkstra.GetIdList().GetId(i))
    return ids


def order_nodes_edgeloop(node_indices: np.ndarray, node_coords: np.ndarray) -> np.ndarray:
    """Order node indices to form closed/continuous loop.

    Parameters
    ----------
    node_indices : np.ndarray
        Array of node indices.
    node_coords : np.ndarray
        Array of node coordinates.

    Returns
    -------
    np.ndarray
        Reordered list of node indices that form a edge loop.

    Notes
    -----
    May not work if mesh density is very anisotropic and does
    not change gradually.
    """
    num_nodes = len(node_indices)

    nodes_coords_to_use = node_coords[node_indices, :]

    idx_visited = [0]
    # iteratively search closest node which has not been visited yet
    iters = 0
    while not len(idx_visited) == num_nodes:
        # for ii in range(0, num_nodes ):
        ref_node = nodes_coords_to_use[idx_visited[-1], :]
        distance = np.linalg.norm(nodes_coords_to_use - ref_node, axis=1)
        indices_sort = np.argsort(distance)
        mask = np.isin(indices_sort, idx_visited, invert=True)
        next_idx = indices_sort[mask][0]
        idx_visited.append(next_idx)
        iters = iters + 1
        if iters > num_nodes:
            raise RecursionError("More iterations needed than expected - check implementation")

    # remap to old numbering and return
    return node_indices[idx_visited]


def sort_edgeloop_anti_clockwise(points_to_sort: np.ndarray, reference_point: np.ndarray) -> bool:
    """Sort the points of an edge-loop in anti-clockwise direction.

    Parameters
    ----------
    points_to_sort : np.bdarray
        Point coordinates used for sorting.
    reference_point : np.ndarray
        Reference point.

    Returns
    -------
    bool
        Flag indicating whether the point order should be reversed or not

    Note
    ----
    This only uses the first two points in the points array, but uses all
    points to compute the center of the sorted points that make up the edge
    loop

    """
    reverse_points = False

    center_of_points = np.mean(points_to_sort, axis=0)
    p1 = points_to_sort[0, :]
    p2 = points_to_sort[1, :]

    # determine direction based on right hand rule
    normal = np.cross(p1 - center_of_points, p2 - center_of_points)
    d1 = np.linalg.norm(center_of_points + normal - reference_point)
    d2 = np.linalg.norm(center_of_points - normal - reference_point)

    if d1 >= d2:
        reverse_points = True
        # normal points away from reference point. reverse order of points

    return reverse_points


def rodrigues_rot(P, n0, n1):
    """Rodrigues rotation.

    Notes
    -----
    - Rotate given points based on a starting and ending vector
    - Axis k and angle of rotation theta given by vectors n0,n1
    P_rot = P*cos(theta) + (k x P)*sin(theta) + k*<k,P>*(1-cos(theta))
    """
    # If P is only 1d array (coords of single point), fix it to be matrix
    if P.ndim == 1:
        P = P[np.newaxis, :]

    # Get vector of rotation k and angle theta
    n0 = n0 / np.linalg.norm(n0)
    n1 = n1 / np.linalg.norm(n1)
    k = np.cross(n0, n1)
    k = k / np.linalg.norm(k)
    theta = np.arccos(np.dot(n0, n1))

    # Compute rotated points
    P_rot = np.zeros((len(P), 3))
    for i in range(len(P)):
        P_rot[i] = (
            P[i] * np.cos(theta)
            + np.cross(k, P[i]) * np.sin(theta)
            + k * np.dot(k, P[i]) * (1 - np.cos(theta))
        )

    return P_rot


def carttopolar(x, y, x0=0, y0=0):
    """Cartisian to polar coordinate system with origin shift to x0,y0."""
    x1 = x - x0
    y1 = y - y0
    r = np.sqrt(x1**2 + y1**2)
    t = np.arctan2(y1, x1) * 180 / math.pi
    if y1 < 0:
        t = 360 + t
    return r, t


def sort_aniclkwise(xy_list, x0=None, y0=None):
    """Sort points anti clockwise with x0 y0 as origin."""
    if x0 is None and y0 is None:
        (x0, y0) = np.mean(xy_list, axis=0).tolist()
    elif x0 is None:
        (x0, _) = np.mean(xy_list, axis=0).tolist()
    elif y0 is None:
        (_, y0) = np.mean(xy_list, axis=0).tolist()
    # print('origin used:', [x0, y0])

    for i in range(len(xy_list)):
        xy_list[i].append(i)

    xy_list1 = sorted(
        xy_list,
        key=lambda a_entry: carttopolar(a_entry[0], a_entry[1], x0, y0)[1],
    )

    sort_index = []
    for x in xy_list1:
        sort_index.append(x[2])
        del x[2]

    return xy_list1, sort_index


def project_3d_points(p_set):
    """Project points on representative plane.

    Notes
    -----
    Uses SVD to find representative plane:
    https://meshlogic.github.io/posts/jupyter/curve-fitting/fitting-a-circle-to-cluster-of-3d-points/
    """
    # -------------------------------------------------------------------------------
    # (1) Fitting plane by SVD for the mean-centered data
    # Eq. of plane is <p,n> + d = 0, where p is a point on plane and n is normal vector
    # -------------------------------------------------------------------------------
    P_mean = np.mean(p_set, axis=0)
    P_centered = p_set - P_mean
    _, _, V = np.linalg.svd(P_centered)
    # Normal vector of fitting plane is given by 3rd column in V
    # Note linalg.svd returns V^T, so we need to select 3rd row from V^T
    normal = V[2, :]

    # -------------------------------------------------------------------------------
    # (2) Project points to coords X-Y in 2D plane
    # -------------------------------------------------------------------------------
    P_xy = rodrigues_rot(P_centered, normal, [0, 0, 1])

    # -------------------------------------------------------------------------------
    # (2) Project points back to the original CS
    # -------------------------------------------------------------------------------
    P_project = np.zeros(p_set.shape)
    for i in range(len(P_xy)):
        point = rodrigues_rot(np.array([P_xy[i, 0], P_xy[i, 1], 0]), [0, 0, 1], normal) + P_mean
        P_project[i] = point.ravel()

    return P_project, normal, P_mean
