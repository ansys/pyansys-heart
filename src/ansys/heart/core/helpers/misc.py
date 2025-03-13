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

"""Module containing some general methods."""

import os

import numpy as np
from scipy.spatial import cKDTree

from ansys.heart.core import LOG as LOGGER
from ansys.heart.core.models import HeartModel


def clean_directory(
    directory: str,
    extensions_to_remove: list[str] = [".stl", ".vtk", ".msh.h5"],
    remove_all: bool = False,
) -> None:
    """Remove files with extension present in the working directory.

    Parameters
    ----------
    extensions_to_remove : List[str], optional
        List of extensions to remove, by default [".stl", ".vtk", ".msh.h5"]
    remove_all: bool, optional
        Flag indicating whether to remove files with any extension.
        Keeps files/folder without extension
    """
    import glob as glob

    files = []
    if not remove_all:
        for ext in extensions_to_remove:
            files += glob.glob(os.path.join(directory, "*" + ext))
    elif remove_all:
        files = glob.glob(os.path.join(directory, "*.*"))

    for file in files:
        try:
            os.remove(file)
        except Exception as e:
            LOGGER.error(f"Unable to delete: {file}. {e}")
    return


def model_summary(model: HeartModel, attributes: list = None) -> dict:
    """Generate a dictionary with model information.

    Parameters
    ----------
    model : HeartModel
        HeartModel for which to generate the summary dictionary
    attributes : list
        List of attributes to try to add to the dict.

    Returns
    -------
    dict
        Dictionary with model information.
    """
    sum_dict = {}
    sum_dict["GENERAL"] = {}

    try:
        sum_dict["GENERAL"]["total_num_tets"] = model.mesh.tetrahedrons.shape[0]
        sum_dict["GENERAL"]["total_num_nodes"] = model.mesh.points.shape[0]
    except TypeError as error:
        LOGGER.error(f"Failed to format General model information. {error}")

    sum_dict["PARTS"] = {}
    sum_dict["CAVITIES"] = {}
    for ii, part in enumerate(model.parts):
        sum_dict["PARTS"][part.name] = {}
        sum_dict["PARTS"][part.name]["num_tets"] = len(part.element_ids)

        sum_dict["PARTS"][part.name]["SURFACES"] = {}
        sum_dict["PARTS"][part.name]["CAPS"] = {}

        for surface in part.surfaces:
            sum_dict["PARTS"][part.name]["SURFACES"][surface.name] = {}
            sum_dict["PARTS"][part.name]["SURFACES"][surface.name]["num_faces"] = (
                surface.triangles.shape[0]
            )

            if attributes:
                for attribute in attributes:
                    try:
                        sum_dict["PARTS"][part.name]["SURFACES"][surface.name][attribute] = getattr(
                            surface.clean(), attribute
                        )
                    except AttributeError:
                        pass

        for cap in part.caps:
            sum_dict["PARTS"][part.name]["CAPS"][cap.name] = {}
            sum_dict["PARTS"][part.name]["CAPS"][cap.name]["num_nodes"] = len(
                cap.global_node_ids_edge
            )

            if attributes:
                for attribute in attributes:
                    try:
                        sum_dict["PARTS"][part.name]["CAPS"][cap.name][attribute] = getattr(
                            cap, attribute
                        )
                    except AttributeError:
                        pass

    for cavity in model.cavities:
        sum_dict["CAVITIES"][cavity.name] = {}
        sum_dict["CAVITIES"][cavity.name]["volume"] = cavity.surface.volume

        if attributes:
            for attribute in attributes:
                try:
                    sum_dict["CAVITIES"][cavity.name][attribute] = getattr(cavity, attribute)
                except AttributeError:
                    pass

    return sum_dict


def rodrigues_rot(points: np.ndarray, n0: np.ndarray, n1: np.ndarray) -> np.ndarray:
    """Perform rodrigues rotation.

    Parameters
    ----------
    points : np.ndarray
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
    if points.ndim == 1:
        points = points[np.newaxis, :]

    # Get vector of rotation k and angle theta
    n0 = n0 / np.linalg.norm(n0)
    n1 = n1 / np.linalg.norm(n1)
    k = np.cross(n0, n1)
    k = k / np.linalg.norm(k)
    theta = np.arccos(np.dot(n0, n1))

    # Compute rotated points
    points_rotation = np.zeros((len(points), 3))
    for i in range(len(points)):
        points_rotation[i] = (
            points[i] * np.cos(theta)
            + np.cross(k, points[i]) * np.sin(theta)
            + k * np.dot(k, points[i]) * (1 - np.cos(theta))
        )

    return points_rotation


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
    point_mean = np.mean(p_set, axis=0)
    point_centered = p_set - point_mean
    _, _, vector = np.linalg.svd(point_centered)
    # Normal vector of fitting plane is given by 3rd column in V
    # Note linalg.svd returns V^T, so we need to select 3rd row from V^T
    normal = vector[2, :]

    # -------------------------------------------------------------------------------
    # (2) Project points to coords X-Y in 2D plane
    # -------------------------------------------------------------------------------
    points_xy = rodrigues_rot(point_centered, normal, [0, 0, 1])

    # -------------------------------------------------------------------------------
    # (2) Project points back to the original CS
    # -------------------------------------------------------------------------------
    point_projected = np.zeros(p_set.shape)
    for i in range(len(points_xy)):
        point = (
            rodrigues_rot(np.array([points_xy[i, 0], points_xy[i, 1], 0]), [0, 0, 1], normal)
            + point_mean
        )
        point_projected[i] = point.ravel()

    return point_projected, normal, point_mean


def _read_orth_element_kfile(fn):
    """Read *ELEMENT_SOLID_ORTHO keywords from file."""

    def get_number_of_elements(file):
        lines = open(file).readlines()
        n = 0
        for line in lines:
            if line[0] == "*":
                n += 1
        return int((len(lines) - n) / 4)

    def generate_specific_rows(file, row_indices):
        with open(file) as f:
            lines = f.readlines()
        return [lines[i] for i in row_indices]

    nele = get_number_of_elements(fn)

    # skip first 1 row and read every 4 row
    skip_row = 1  # because the first line is *ELEMENT_SOLID_ORTHO
    every_row = 4

    # element ID and part ID
    index = np.linspace(0, nele - 1, num=nele, dtype=int) * every_row + skip_row
    data = generate_specific_rows(fn, index)
    ids = np.loadtxt(data, dtype=int)[:, :]

    # element connectivity
    index = np.linspace(0, nele - 1, num=nele, dtype=int) * every_row + skip_row + 1
    data = generate_specific_rows(fn, index)
    connect = np.loadtxt(data, dtype=int)[:, :]

    # fiber
    index = np.linspace(0, nele - 1, num=nele, dtype=int) * every_row + skip_row + 2
    data = generate_specific_rows(fn, index)
    fib = np.loadtxt(data)
    # sheet
    index = np.linspace(0, nele - 1, num=nele, dtype=int) * every_row + skip_row + 3
    data = generate_specific_rows(fn, index)
    sheet = np.loadtxt(data)

    # sort by ids
    index = np.argsort(ids[:, 0])
    elem_ids = ids[index, 0]
    part_ids = ids[index, 1]
    connect = connect[index]
    fib = fib[index]
    sheet = sheet[index]

    return elem_ids, part_ids, connect, fib, sheet


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
