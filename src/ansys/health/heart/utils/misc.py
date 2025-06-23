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

"""Module containing miscellaneous methods."""

import os

import numpy as np
from scipy.spatial import cKDTree

from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.models import HeartModel
from ansys.health.heart.parts import Chamber


def clean_directory(
    directory: str,
    extensions_to_remove: list[str] = [".stl", ".vtk", ".msh.h5"],
    remove_all: bool = False,
) -> None:
    """Remove files from the working directory with given extensions.

    Parameters
    ----------
    extensions_to_remove : List[str], default: ``['.stl', '.vtk', '.msh.h5']``
        List of extensions to remove.
    remove_all: bool, default: False
        Whether to remove files with any extension. Files without extensions are kept.
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
        Heart model to generate the summary dictionary for.
    attributes : list
        List of attributes to add to the dictionary.

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
        LOGGER.error(f"Failed to format general model information. {error}")

    sum_dict["PARTS"] = {}
    sum_dict["CAVITIES"] = {}
    for ii, part in enumerate(model.parts):
        sum_dict["PARTS"][part.name] = {}
        sum_dict["PARTS"][part.name]["num_tets"] = len(part.get_element_ids(model.mesh))

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
        if isinstance(part, Chamber):
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


def _read_orth_element_kfile(
    fn,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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
    """Spherical Linear Interpolation between two unit vectors, v0 and v1."""
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
        N x 3 array of source points coordinates.
    source_vec : np.ndarray
        N x 3 array of source vectors.
    target_pos : np.ndarray
        M x 3 array of target points coordinates.

    Returns
    -------
    np.ndarray
        M x 3 array of target vectors
    """
    # legal test
    norm = np.linalg.norm(source_vec, axis=1)
    if not np.allclose(norm, 1.0):
        raise TypeError("Slerp interpolation requires unit vectors.")

    # Build a KD-tree once
    tree = cKDTree(source_pos)

    def interpolate_with_k_nearest(query_point: np.ndarray, k: int = 4) -> np.ndarray:
        """Slerp interpolate with k nearest points.

        Parameters
        ----------
        query_point : np.ndarray
            Query point coordinate.
        k : int, default: 4
            Number of nearest points to use.

        Returns
        -------
        np.ndarray
            Vector on query point.
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
