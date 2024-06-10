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

# import copy  # noqa
# import pathlib
# from typing import List, Optional, Tuple, Union

import copy

from ansys.heart.core import LOG as LOGGER

# import ansys.heart.preprocessor.mesh.connectivity as connect
# import ansys.heart.preprocessor.mesh.geodisc as geodisc
# import ansys.heart.preprocessor.mesh.vtkmethods as vtkmethods
# from ansys.heart.simulator.settings.material.material import MechanicalMaterialModel
import numpy as np

try:
    import pyvista as pv
except ImportError:
    LOGGER.warning("Importing pyvista failed. Install with: pip install pyvista")


def _get_fill_data(
    mesh1: pv.UnstructuredGrid | pv.PolyData,
    mesh2: pv.UnstructuredGrid | pv.PolyData,
    array_name: str,
    array_association: str = "cell",
    pad_value_int: int = None,
    pad_value_float: float = None,
):
    if array_name not in mesh1.array_names:
        return

    if array_association == "cell":
        if array_name in mesh2.cell_data.keys():
            return mesh2.cell_data[array_name]

        array = mesh1.cell_data[array_name]
        n_pads = mesh2.n_cells

    elif array_association == "point":
        if array_name in mesh2.point_data.keys():
            return mesh2.point_data[array_name]

        array = mesh1.point_data[array_name]
        n_pads = mesh2.n_points

    shape = list(array.shape)
    shape[0] = n_pads
    shape = tuple(shape)

    pad_array = np.zeros(shape, dtype=array.dtype)

    if isinstance(array[0], (np.float64, np.float32)):
        if not pad_value_float:
            pad_array = pad_array * np.nan
        else:
            pad_array = pad_array * pad_value_float

    elif isinstance(array[0], (np.int32, np.int64)):
        if pad_value_int:
            pad_array = pad_array + pad_value_int

    return pad_array


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

    @property
    def triangles(self):
        """Get all triangles of the mesh."""
        return self.cells_dict[pv.CellType.TRIANGLE]

    @property
    def lines(self):
        """Get all triangles of the mesh."""
        return self.cells_dict[pv.CellType.LINE]

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

        # self.triangles: np.ndarray = None
        """Faces that make up the tetrahedrons."""
        self.face_types: np.ndarray = None
        """Type of face: 1: interior face, 2: boundary face, 3: interface face."""
        self.conn = {"c0": [], "c1": []}
        """Face-tetra connectivity array."""
        # TODO: just store used nodes in interfaces and boundaries
        # and add mapper to map from local to global (volume mesh) node ids
        # self.interfaces: List[SurfaceMesh] = []
        """List of surface meshes that make up the interface between different parts."""
        # self.boundaries: List[SurfaceMesh] = []
        """List of boundary surface meshes within the part."""
        pass

    @property
    def surface_ids(self) -> np.ndarray:
        """Unique surface ids.

        Returns
        -------
        np.ndarray
            Array with unique surface ids
        """
        try:
            mask = self.celltypes == pv.CellType.TRIANGLE
            mask1 = self.cell_data["surface-id"] > -1
            np.vstack((mask, mask1))
            mask = np.all(np.vstack((mask, mask1)), axis=0)
            return np.unique(self.cell_data["surface-id"][mask])
        except KeyError:
            return None

    @property
    def volume_ids(self) -> np.ndarray:
        """Unique volume ids.

        Returns
        -------
        np.ndarray
            Array with unique volume ids
        """
        try:
            mask = self.celltypes == pv.CellType.TETRA
            mask1 = self.cell_data["volume-id"] > -1
            np.vstack((mask, mask1))
            mask = np.all(np.vstack((mask, mask1)), axis=0)
            return np.unique(self.cell_data["volume-id"][mask])
        except KeyError:
            return None

    @property
    def beam_ids(self) -> np.ndarray:
        """Unique beam ids.

        Returns
        -------
        np.ndarray
            Array with beam surface ids
        """
        try:
            mask = self.celltypes == pv.CellType.LINE
            mask1 = self.cell_data["beam-id"] > -1
            np.vstack((mask, mask1))
            mask = np.all(np.vstack((mask, mask1)), axis=0)
            return np.unique(self.cell_data["beam-id"][mask])
        except KeyError:
            return None

    def clean(self, *args):
        """Merge duplicate points and return cleaned copy."""
        self_c = copy.deepcopy(self)
        super(Mesh, self_c).__init__(pv.UnstructuredGrid(self).clean(*args))
        return self_c

    def add_mesh(
        self,
        mesh: pv.PolyData | pv.UnstructuredGrid,
        keep_data: bool = True,
        fill_float: np.float64 = np.nan,
        fill_int: int = -1,
    ):
        """Add another mesh to this object.

        Notes
        -----
        Adding the mesh is always in-place

        Parameters
        ----------
        mesh : pv.PolyData | pv.UnstructuredGrid
            Mesh to add, either PolyData or UnstructuredGrid
        keep_data : bool, optional
            Flag specifying whether to try to keep mesh point/cell data, by default True
        """
        if keep_data:
            # add cell/point arrays in self
            cell_data_names = [k for k in mesh.cell_data.keys()]
            point_data_names = [k for k in mesh.point_data.keys()]

            for name in cell_data_names:
                self.cell_data[name] = _get_fill_data(
                    mesh, self, name, "cell", fill_int, fill_float
                )

            for name in point_data_names:
                self.point_data[name] = _get_fill_data(
                    mesh, self, name, "point", fill_int, fill_float
                )

            # add cell/point arrays mesh to be added
            cell_data_names = [k for k in self.cell_data.keys()]
            point_data_names = [k for k in self.point_data.keys()]

            for name in cell_data_names:
                mesh.cell_data[name] = _get_fill_data(self, mesh, name, "cell")

            for name in point_data_names:
                mesh.point_data[name] = _get_fill_data(self, mesh, name, "point")

        # super().__init__(self + mesh)
        merged = pv.merge((self, mesh), merge_points=False)
        super().__init__(merged)
        return self
