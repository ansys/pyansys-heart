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
import numpy as np


def rodrigues_rot(P: np.ndarray, n0: np.ndarray, n1: np.ndarray) -> np.ndarray:
    """Perform rodrigues rotation.

    Parameters
    ----------
    P : np.ndarray
        Points to rotate.
    n0 : np.ndarray
        Vector 1.
    n1 : np.ndarray
        Vector 2.

    Notes
    -----
    Rotate given points based on a starting and ending vector.
    Axis k and angle of rotation theta given by vectors n0,n1.
    P_rot = P*cos(theta) + (k x P)*sin(theta) + k*<k,P>*(1-cos(theta))

    Returns
    -------
    np.ndarray
        Rotated points.
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


def project_3d_points(p_set: np.ndarray) -> np.ndarray:
    """Project points on a representative plane.

    Parameters
    ----------
    p_set : np.ndarray
        Point set, Nx3

    Notes
    -----
    Uses SVD to find representative plane:
    https://meshlogic.github.io/posts/jupyter/curve-fitting/fitting-a-circle-to-cluster-of-3d-points/

    Returns
    -------
    np.ndarray
        Projected points onto the SVD representative plane.
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
