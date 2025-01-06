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

"""
Functions for fibers interpolation.

Note: if interpolate both fiber and sheet, make sure they should be orthogonal.
"""

import numpy as np
from scipy.spatial import cKDTree


def _slerp(v0: np.ndarray, v1: np.ndarray, t: float) -> np.ndarray:
    """Spherical Linear Interpolation between two unit vectors v0 and v1."""
    # Compute dot product and clamp to handle numerical issues
    dot = np.dot(v0, v1)
    dot = np.clip(dot, -1.0, 1.0)

    # Compute the angle between the vectors
    theta = np.arccos(dot)

    # If the angle is very small, linear interpolation suffices
    if np.isclose(theta, 0):
        return (1 - t) * v0 + t * v1

    # Compute SLERP
    sin_theta = np.sin(theta)
    v_out = (np.sin((1 - t) * theta) / sin_theta) * v0 + (np.sin(t * theta) / sin_theta) * v1
    return v_out


def interpolate_slerp(
    source_pos: np.ndarray, source_vec: np.ndarray, target_pos: np.ndarray
) -> np.ndarray:
    """Spherical linear interpolation.

    Parameters
    ----------
    source_pos : np.ndarray
        N x 3 array of source points coordinates
    source_vec : np.ndarray
        N x 3 array of source vectors
    target_pos : np.ndarray
        M x 3 array of target points coordinates

    Returns
    -------
    np.ndarray
        M x 3 array of target vectors
    """
    # legal test
    norm = np.linalg.norm(source_vec, axis=1)
    if not np.allclose(norm, 1.0):
        raise TypeError("slerp interpolation must be used for unit vectors.")

    # Build a KD-tree once
    tree = cKDTree(source_pos)

    def interpolate_with_k_nearest(query_point: np.ndarray, k: int = 4) -> np.ndarray:
        """Slerp interpolate with k nearest points.

        Parameters
        ----------
        query_point : np.ndarray
            query point coordinate
        k : int, optional
            no. of nearest points to be used, by default 4

        Returns
        -------
        np.ndarray
            vector on query point
        """
        # Find the k-nearest neighbors
        distances, indices = tree.query(query_point, k=k)

        # nearest vectors
        nearest_vectors = source_vec[indices]

        # inverse-distance weights
        weights = 1 / (distances + 1e-8)
        weights /= np.sum(weights)

        # Perform SLERP interpolation using weights
        interpolated_vector = nearest_vectors[0]
        for i in range(1, k):
            interpolated_vector = _slerp(interpolated_vector, nearest_vectors[i], weights[i])

        # Normalize
        return interpolated_vector / np.linalg.norm(interpolated_vector)

    result = np.zeros_like(target_pos)
    for i in range(target_pos.shape[0]):
        print(f"{i} / {target_pos.shape[0]}")
        result[i] = interpolate_with_k_nearest(target_pos[i])

    return result
