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

"""Module that defines classes used in the heart model."""

from __future__ import annotations

import copy
from enum import Enum
import json
import os
import pathlib
from typing import List, Literal, Union

from deprecated import deprecated
import numpy as np
import pyvista as pv

from ansys.health.heart import LOG as LOGGER
import ansys.health.heart.utils.vtk_utils as vtk_utils

_SURFACE_CELL_TYPES = [pv.CellType.QUAD, pv.CellType.TRIANGLE]
_VOLUME_CELL_TYPES = [pv.CellType.HEXAHEDRON, pv.CellType.TETRA]


def _get_fill_data(
    mesh1: Union[pv.UnstructuredGrid, pv.PolyData],
    mesh2: Union[pv.UnstructuredGrid, pv.PolyData],
    array_name: str,
    array_association: str = "cell",
    pad_value_int: int = None,
    pad_value_float: float = None,
) -> np.ndarray:
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


def _get_global_cell_ids(mesh: pv.UnstructuredGrid, celltype: pv.CellType) -> np.ndarray:
    """Get the global cell iID of a given cell type.

    Parameters
    ----------
    mesh : pv.UnstructuredGrid
        Unstructured grid to obtain the global cell IDs from.
    celltype : pv.CellType
        Cell type to get global cell IDs of.

    Returns
    -------
    np.ndarray
        Array with global cell IDs.
    """
    return np.argwhere(np.isin(mesh.celltypes, celltype)).flatten()


def _invert_dict(dictionary: dict) -> dict:
    """Invert a dictionary.

    Parameters
    ----------
    dict : dict
        Dictionary to invert.

    Returns
    -------
    dict
        Inverted dictionary.

    """
    if dictionary == {}:
        return {}
    else:
        return {v: k for k, v in dictionary.items()}


def _convert_int64_to_int32(
    mesh: Mesh | pv.UnstructuredGrid | pv.PolyData | SurfaceMesh, array_names: list[str] = None
):
    """Change the datatype of cell and point arrays to int32.

    Parameters
    ----------
    mesh : Mesh | pv.UnstructuredGrid | pv.PolyData | SurfaceMesh
        The ``PyVista`` mesh object containing the cell and point data arrays to convert.
    array_names : list[str], default: None
        List of specific array names to convert. If not provided, all arrays in
        the mesh will be checked and converted if necessary.

    Returns
    -------
    Mesh | pv.UnstructuredGrid | pv.PolyData | SurfaceMesh
        The input mesh with cell and point data arrays converted to ``int32``.

    Notes
    -----
    ``PyVista`` uses ``int64`` by default, which is not compatible with
    the ``PyVista`` sphinx plot directive for interactive plots in the documentation.
    This method changes the datatype of the cell and point arrays to ``int32``.
    """
    if array_names is None:
        array_names = mesh.array_names

    cell_data_keys = mesh.cell_data.keys()
    point_data_keys = mesh.point_data.keys()

    for array_name in mesh.array_names:
        try:
            # convert cell data to int32
            if array_name in cell_data_keys and mesh.cell_data[array_name].dtype == np.int64:
                mesh.cell_data[array_name] = mesh.cell_data[array_name].astype(np.int32)

            # convert point data to int32
            if array_name in point_data_keys and mesh.point_data[array_name].dtype == np.int64:
                mesh.point_data[array_name] = mesh.point_data[array_name].astype(np.int32)

        except Exception as err:
            LOGGER.debug(f"Failed to convert {array_name} to int32. {err}")
    return mesh


# TODO: Deprecate
class Feature:
    """Feature class."""

    def __init__(self, name: str = None) -> None:
        #! This class can be deprecated.
        self.name = name
        """Name of feature."""
        self.type = None
        """Type of feature."""
        self._node_set_id: int = None
        """Nodeset ID associated with feature."""
        self._seg_set_id: int = None
        """Segment set ID associated with feature."""
        self.pid: int = None
        """Part ID associated with feature."""

        pass


class SurfaceMesh(pv.PolyData):
    """Surface class."""

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
                    "Assigning less nodes than the original. Not implemented yet."
                )

        except Exception as e:
            LOGGER.error(f"Failed to set nodes. {e}")
            return

    @property
    def triangles(self):
        """Triangular faces of the surface ``num_faces`` x 3."""
        faces = np.reshape(self.faces, (self.n_cells, 3 + 1))[:, 1:]
        return faces

    @triangles.setter
    def triangles(self, value: np.ndarray):
        # sets faces of PolyData
        try:
            num_faces = value.shape[0]
            faces = np.hstack([np.full((num_faces, 1), 3, dtype=np.int8), value])
            self.faces = faces
        except Exception:
            return

    @property
    def triangles_global(self):
        """Global triangle IDs.

        Returns
        -------
        Tries to use ``point_data["_global-point-ids"]`` to retrieve
        triangle definitions in global IDs.
        """
        return self.point_data["_global-point-ids"][self.triangles]

    @property
    def boundary_edges(self):
        """Boundary edges of self."""
        boundary_edges = vtk_utils.get_boundary_edge_loops(self, remove_open_edge_loops=False)
        boundary_edges = np.vstack(list(boundary_edges.values()))
        return boundary_edges

    @property
    def boundary_edges_global(self):
        """Global point IDs of boundary edges."""
        return self.point_data["_global-point-ids"][self.boundary_edges]

    def __init__(
        self,
        var_inp: Union[pv.PolyData, np.ndarray, list, str, pathlib.Path] = None,
        name: str = None,
        triangles: np.ndarray = None,
        nodes: np.ndarray = None,
        id: int = None,
        **kwargs,
    ) -> None:
        # *NOTE: pv.PolyData supports variable input through the first argument (var_inp)
        # * the following is to make sure this object behaves similar to pv.PolyData
        # * https://github.com/pyvista/pyvista/blob/release/0.44/pyvista/core/pointset.py#L500-L1693
        # * https://docs.pyvista.org/version/stable/api/core/_autosummary/pyvista.polydata#pyvista.PolyData # noqa E501

        if isinstance(var_inp, (pv.PolyData, np.ndarray, list, str, pathlib.Path)):
            kwargs["var_inp"] = var_inp

        super(SurfaceMesh, self).__init__(**kwargs)

        self.name = name
        """Name of the surface."""

        self.id: int = id
        """ID of the surface."""

        self.triangles = triangles
        """Triangular faces of the surface ``num_faces`` x 3."""
        self.nodes = nodes
        """Node coordinates."""
        self._seg_set_id: int = None
        """Segment set ID."""
        self._node_set_id: int = None
        """Nodeset ID."""

    @property
    def node_ids_triangles(self) -> np.ndarray:
        """Local node IDs sorted by earliest occurrence."""
        _, idx = np.unique(self.triangles.flatten(), return_index=True)
        node_ids = self.triangles.flatten()[np.sort(idx)]
        return node_ids

    @property
    def global_node_ids_triangles(self):
        """Global node IDs from point data."""
        return self.point_data["_global-point-ids"][self.node_ids_triangles]

    @property
    def _boundary_nodes(self) -> np.ndarray:
        """Global node IDs of nodes on the boundary of the mesh (if any)."""
        _, idx = np.unique(self.boundary_edges.flatten(), return_index=True)
        node_ids = self.boundary_edges.flatten()[np.sort(idx)]
        return node_ids

    def force_normals_inwards(self):
        """Force the cell ordering of the closed surface such that normals point inward."""
        if not self.is_manifold:
            LOGGER.warning("Surface is non-manifold.")

        #! Flip normals and consistent normals should enforce that normals are pointing
        #! inwards for a manifold surface. See:
        #! https://docs.pyvista.org/api/core/_autosummary/pyvista.polydatafilters.compute_normals
        #! With 0.44.1 we may need to remove the normals prior to computing them. With earlier
        #! versions this seems to have unintentionally worked.
        try:
            self.cell_data.remove("Normals")
        except KeyError:
            pass
        try:
            self.point_data.remove("Normals")
        except KeyError:
            pass

        self.compute_normals(inplace=True, auto_orient_normals=True, flip_normals=True)
        return self


class Cavity(Feature):
    """Cavity class."""

    def __init__(self, surface: SurfaceMesh = None, centroid: np.ndarray = None, name=None) -> None:
        super().__init__(name)

        #! that that if we don't do a deepcopy the associated algorithms may
        #! modify the cells/points in the original object!!
        self.surface: SurfaceMesh = copy.deepcopy(surface)
        """Surface mesh making up the cavity."""
        self.centroid: np.ndarray = centroid
        """Centroid of the cavity."""

    @property
    def volume(self):
        """Volume of the cavity."""
        self.surface.force_normals_inwards()
        return self.surface.volume

    def compute_centroid(self):
        """Compute the centroid of the cavity."""
        # self.centroid = np.mean(self.surface.nodes[np.unique(self.surface.triangles), :], axis=0)
        self.centroid = self.surface.center
        return self.centroid


# Naming convention of caps.
class CapType(Enum):
    """Enumeration tracking cap names."""

    MITRAL_VALVE = "mitral-valve"
    """Cap representing mitral valve region."""
    AORTIC_VALVE = "aortic-valve"
    """Cap representing aortic valve region."""
    MITRAL_VALVE_ATRIUM = "mitral-valve-atrium"
    """Cap representing mitral valve region on the atrial side."""
    COMBINED_MITRAL_AORTIC_VALVE = "combined-mitral-aortic-valve"
    """Combined mitral aortic valve. Valid for truncated models."""
    PULMONARY_VALVE = "pulmonary-valve"
    """Cap representing pulmonary valve region."""
    TRICUSPID_VALVE = "tricuspid-valve"
    """Cap representing tricuspid valve region."""
    TRICUSPID_VALVE_ATRIUM = "tricuspid-valve-atrium"
    """Cap representing tricuspid valve region on the atrial side."""

    LEFT_ATRIUM_APPENDAGE = "left-atrium-appendage"
    """Cap representing left atrium appendage region."""
    LEFT_SUPERIOR_PULMONARY_VEIN = "left-superior-pulmonary-vein"
    """Cap representing left superior pulmonary vein region."""
    LEFT_INFERIOR_PULMONARY_VEIN = "left-inferior-pulmonary-vein"
    """Cap representing left inferior pulmonary vein region."""
    RIGHT_INFERIOR_PULMONARY_VEIN = "right-inferior-pulmonary-vein"
    """Cap representing right inferior pulmonary vein region."""
    RIGHT_SUPERIOR_PULMONARY_VEIN = "right-superior-pulmonary-vein"
    """Cap representing right superior pulmonary vein region."""
    SUPERIOR_VENA_CAVA = "superior-vena-cava"
    """Cap representing superior vena cava region."""
    INFERIOR_VENA_CAVA = "inferior-vena-cava"
    """Cap representing inferior vena cava region."""
    UNKNOWN = "unknown-cap"
    """Cap with unknown association."""


class Cap(Feature):
    """Cap class."""

    @property
    def _local_node_ids_edge(self):
        """Local node IDs of the cap edge."""
        edges = vtk_utils.get_boundary_edge_loops(self._mesh)
        edge_local_ids = np.unique(np.array([np.array(edge) for edge in edges.values()]))
        return edge_local_ids

    @property
    def global_node_ids_edge(self):
        """Global node IDs of the edge of the cap."""
        return self._mesh.point_data["_global-point-ids"][self._local_node_ids_edge]

    @property
    def _local_centroid_id(self):
        """Local ID of the centroid."""
        centroid_id = np.setdiff1d(np.arange(0, self._mesh.n_points), self._local_node_ids_edge)
        if len(centroid_id) != 1:
            LOGGER.error("Failed to identify single centroid node.")
            return None

        return centroid_id[0]

    @property
    def global_centroid_id(self):
        """Global centroid ID."""
        return self._mesh.point_data["_global-point-ids"][self._local_centroid_id]

    @property
    def centroid(self):
        """Centroid of the cap."""
        return self._mesh.points[self._local_centroid_id, :]

    @property
    def cap_normal(self):
        """Compute mean normal of the cap."""
        return np.mean(self._mesh.compute_normals().cell_data["Normals"], axis=0)

    def __init__(
        self,
        name: str = None,
        cap_type: CapType = None,
    ) -> None:
        super().__init__(name)
        """Centroid of the cap ID (in case centroid node is created)."""
        self._mesh: SurfaceMesh = None

        if cap_type is None or isinstance(cap_type, CapType):
            self.type = cap_type
        else:
            LOGGER.warning(f"Failed to set cap type for {name}, {cap_type}.")

        return


class Point(Feature):
    """Point class, which can be used to collect relevant points in the mesh."""

    def __init__(self, name: str = None, xyz: np.ndarray = None, node_id: int = None) -> None:
        super().__init__(name)

        self.xyz: np.ndarray = xyz
        """XYZ coordinates of the point."""
        self.node_id: int = node_id
        """Global node ID of the point."""


class Mesh(pv.UnstructuredGrid):
    """Mesh class, which inherits from the PyVista unstructured grid object.

    Notes
    -----
    This class inherits from the ``pyvista.UnstructuredGrid`` object and adds additional
    attributes and convenience methods for enhanced functionality. The ``_volume_id``,
    ``_surface_id``, and ``_line_id`` cell arrays keep track of *labeled* selections of
    cells. The ``_volume_id`` cell array is used to group 3D volume cells together.
    Any non-3D volume cell is labeled as ``numpy.nan``. Similarly 2D and 1D cells are tracked
    through the ``_surface_id`` and ``_line_id`` cell arrays respectively.
    """

    @property
    def tetrahedrons(self):
        """Tetrahedrons ``num_tetra`` x 4."""
        return self.cells_dict[pv.CellType.TETRA]

    @property
    def triangles(self):
        """All triangles of the mesh."""
        return self.cells_dict[pv.CellType.TRIANGLE]

    @property
    def lines(self):
        """Get all triangles of the mesh."""
        return self.cells_dict[pv.CellType.LINE]

    @property
    def _surfaces(self) -> List[SurfaceMesh]:
        """List of surfaces in the mesh."""
        if self.surface_ids is None:
            return []
        surfaces = []
        for sid in self.surface_ids:
            surface = SurfaceMesh(self.get_surface(sid))
            surface.id = sid
            try:
                surface.name = self._surface_id_to_name[sid]
            except KeyError as error:
                LOGGER.debug(f"Failed to give surface with ID {sid} a name. {error}")
            surfaces.append(surface)
        return surfaces

    @property
    def _volumes(self):
        """List of volumes in the mesh."""
        if self.volume_ids is None:
            return []
        return [self.get_volume(volume_id) for volume_id in self.volume_ids]

    @property
    def _global_triangle_ids(self):
        """Global IDs of triangular cells."""
        return _get_global_cell_ids(self, pv.CellType.TRIANGLE)

    @property
    def _global_tetrahedron_ids(self):
        """Global IDs of tetrahedral cells."""
        return _get_global_cell_ids(self, pv.CellType.TETRA)

    @property
    def surface_ids(self) -> np.ndarray:
        """Unique surface IDs.

        Returns
        -------
        np.ndarray
            NumPy array with unique surface IDs.
        """
        try:
            mask = np.isin(self.celltypes, _SURFACE_CELL_TYPES)
            mask1 = np.invert(np.isnan(self.cell_data["_surface-id"]))
            mask = np.all(np.vstack((mask, mask1)), axis=0)
            return np.unique(self.cell_data["_surface-id"][mask])
        except KeyError:
            LOGGER.debug(f"Failed to extract one of {_SURFACE_CELL_TYPES}.")
            return []

    @property
    def surface_names(self) -> List[str]:
        """List of surface names."""
        return [v for k, v in self._surface_id_to_name.items()]

    @property
    def volume_ids(self) -> np.ndarray:
        """NumPy array with unique volume IDs.

        Returns
        -------
        np.ndarray
            NumPy array with unique volume IDs.
        """
        try:
            mask = np.isin(self.celltypes, _VOLUME_CELL_TYPES)
            mask1 = np.invert(np.isnan(self.cell_data["_volume-id"]))
            mask = np.all(np.vstack((mask, mask1)), axis=0)
            return np.unique(self.cell_data["_volume-id"][mask])
        except KeyError:
            LOGGER.debug(f"Failed to extract one of {_VOLUME_CELL_TYPES}.")
            return None

    @property
    def _unused_volume_id(self) -> int:
        """Get unused volume ID."""
        if self.volume_ids is None:
            return 1
        return int(np.max(self.volume_ids) + 1)

    @property
    def _unused_surface_id(self) -> int:
        """Get unused surface ID."""
        if self.surface_ids is None:
            return 1
        return int(np.max(self.surface_ids) + 1)

    @property
    def _unused_line_id(self) -> int:
        """Get unused line ID."""
        if self.line_ids is None:
            return 1
        return int(np.max(self.line_ids) + 1)

    @property
    def volume_names(self) -> List[str]:
        """List of volume names."""
        return [v for k, v in self._volume_id_to_name.items()]

    @property
    def line_ids(self) -> np.ndarray:
        """NumPy array with unique line IDs.

        Returns
        -------
        np.ndarray
            NumPy array with unique line IDs.
        """
        try:
            mask = self.celltypes == pv.CellType.LINE
            mask1 = np.invert(np.isnan(self.cell_data["_line-id"]))
            mask = np.all(np.vstack((mask, mask1)), axis=0)
            return np.unique(self.cell_data["_line-id"][mask])
        except KeyError:
            return None

    @property
    def line_names(self) -> List[str]:
        """List of volume names."""
        return [v for k, v in self._line_id_to_name.items()]

    @property
    def _surface_name_to_id(self):
        return _invert_dict(self._surface_id_to_name)

    @property
    def _volume_name_to_id(self):
        return _invert_dict(self._volume_id_to_name)

    @property
    def _line_name_to_id(self):
        return _invert_dict(self._line_id_to_name)

    @property
    def _global_cell_ids(self):
        """Global cell IDs."""
        self._set_global_ids()
        return self.cell_data["_global-cell-ids"]

    @property
    def _global_point_ids(self):
        """Global point IDs."""
        self._set_global_ids()
        return self.point_data["_global-point-ids"]

    def __init__(self, *args):
        super().__init__(*args)

        self._surface_id_to_name: dict = {}
        """Surface ID to name map."""
        self._volume_id_to_name: dict = {}
        """Volume ID to name map."""
        self._line_id_to_name: dict = {}
        """Line ID to name map."""
        pass

    def _add_mesh(
        self,
        mesh_input: Union[pv.PolyData, pv.UnstructuredGrid],
        keep_data: bool = True,
        fill_float: np.float64 = np.nan,
        fill_int: int = -1,
        merge_points: bool = False,
    ):
        """Add another mesh to this object.

        Notes
        -----
        Adding the mesh is always in place.

        Parameters
        ----------
        mesh_input : pv.PolyData | pv.UnstructuredGrid
            Mesh to add, either ``PolyData`` or ``UnstructuredGrid``.
        keep_data : bool, default: True
            Whether to try to keep mesh point/cell data.
        merge_points : bool, default: False
            Flag specifying whether to merge the points.
        """
        mesh = copy.copy(mesh_input)
        # NOTE: PyVista 0.45.0 sometimes has more data cell/point data arrays than number of
        # cells/points. This seems to happen mostly in PolyData objects, casting to an unstructured
        # grid seems to fix this. Aternatively we can call clean and deactivate all flags.
        # However, that may have other side effects.
        if isinstance(mesh, pv.PolyData):
            mesh = mesh.cast_to_unstructured_grid()

        if keep_data:
            # add cell/point arrays in self
            cell_data_names = [k for k in mesh.cell_data.keys()]
            point_data_names = [k for k in mesh.point_data.keys()]

            for name in cell_data_names:
                fill_data = _get_fill_data(mesh, self, name, "cell", fill_int, fill_float)
                if isinstance(fill_data, np.ndarray) and fill_data.shape[0] == 0:
                    continue
                self.cell_data[name] = fill_data

            for name in point_data_names:
                fill_data = _get_fill_data(mesh, self, name, "point", fill_int, fill_float)
                if isinstance(fill_data, np.ndarray) and fill_data.shape[0] == 0:
                    continue
                self.point_data[name] = fill_data

            # add cell/point arrays mesh to be added
            cell_data_names = [k for k in self.cell_data.keys()]
            point_data_names = [k for k in self.point_data.keys()]

            for name in cell_data_names:
                fill_data = _get_fill_data(self, mesh, name, "cell")
                if isinstance(fill_data, np.ndarray) and fill_data.shape[0] == 0:
                    continue
                mesh.cell_data[name] = fill_data

            for name in point_data_names:
                fill_data = _get_fill_data(self, mesh, name, "point")
                if isinstance(fill_data, np.ndarray) and fill_data.shape[0] == 0:
                    continue
                mesh.point_data[name] = fill_data

        merged = pv.merge((self, mesh), merge_points=merge_points, main_has_priority=False)
        super().__init__(merged)
        return self

    def _set_global_ids(self):
        """Add global cell and point IDs as cell and point data array."""
        self.cell_data["_global-cell-ids"] = np.array(np.arange(0, self.n_cells), dtype=int)
        self.point_data["_global-point-ids"] = np.array(np.arange(0, self.n_points), dtype=int)
        return

    def _get_submesh(
        self, sid: int, scalar: Literal["_surface-id", "_line-id", "_volume-id"]
    ) -> pv.UnstructuredGrid:
        # NOTE: extract_cells cleans the object, removing any unused points.
        if scalar not in self.cell_data.keys():
            LOGGER.debug(f"{scalar} does not exist in 'cell_data'.")
            return None
        mask = np.isin(self.cell_data[scalar], sid)
        self._set_global_ids()
        return self.extract_cells(mask)

    def _get_duplicate_surface_names(self):
        names, counts = np.unique(self.surface_names, return_counts=True)
        return names[counts > 1]

    def _get_duplicate_volume_names(self):
        names, counts = np.unique(self.volume_names, return_counts=True)
        return names[counts > 1]

    def _get_duplicate_line_names(self):
        names, counts = np.unique(self.line_names, return_counts=True)
        return names[counts > 1]

    def _get_unmapped_volumes(self):
        if self.volume_ids is None:
            return []
        unmapped_ids = self.volume_ids[
            np.invert(np.isin(self.volume_ids, list(self._volume_id_to_name.keys())))
        ]
        return unmapped_ids

    def _get_unmapped_surfaces(self):
        if self.surface_ids is None:
            return []
        unmapped_ids = self.surface_ids[
            np.invert(np.isin(self.surface_ids, list(self._surface_id_to_name.keys())))
        ]
        return unmapped_ids

    def _get_unmapped_lines(self):
        if self.line_ids is None:
            return []
        unmapped_ids = self.line_ids[
            np.invert(np.isin(self.line_ids, list(self._line_id_to_name.keys())))
        ]
        return unmapped_ids

    def get_unused_id_in_range(
        self, id_type: Literal["volume", "surface", "line"], start: int = 1, end: int = 1000
    ) -> int:
        """
        Get an unused ID within a specified range.

        Parameters
        ----------
        id_type : Literal["volume", "surface", "line"]
            The type of ID for which to retrieve an unused value.
        start : int, default: 1
            Start of the range (inclusive).
        end : int, default: 1000
            End of the range (inclusive).

        Returns
        -------
        int
            An unused and unique ID within the specified range.

        Raises
        ------
        ValueError
            If no unused ID is found in the specified range.

        """
        if id_type == "volume":
            ids = self.volume_ids
        elif id_type == "surface":
            ids = self.surface_ids
        elif id_type == "line":
            ids = self.line_ids
        else:
            raise ValueError(f"Invalid id_type: {id_type}. Must be 'volume', 'surface', or 'line'.")

        if ids is None:
            return start

        # Subtract existing IDs from IDs in range to get candidate IDs
        candidate_ids = set(range(start, end + 1)) - set(ids)
        if candidate_ids != {}:
            return min(candidate_ids)
        else:
            raise ValueError(f"No unused ID found in the range {start} to {end}.")

    def save(self, filename: Union[str, pathlib.Path], **kwargs):
        """Save mesh."""
        super(Mesh, self).save(filename, **kwargs)
        extension = pathlib.Path(filename).suffix
        self._save_id_to_name_map(filename.replace(extension, ".namemap.json"))
        return

    def load_mesh(self, filename: Union[str, pathlib.Path]):
        """Load an existing mesh.

        Notes
        -----
        This method tries to read a JSON file with the volume/surface ID to a name map
        with extension ``.namemap.json`` in the same directory as the file. Alternatively,
        you can read the name map manually by calling ``._load_id_to_name_map(filename)``.

        Parameters
        ----------
        filename : Union[str, pathlib.Path]
            Full path to the mesh file.
        """
        super(Mesh, self).__init__(filename)
        extension = pathlib.Path(filename).suffix
        filename_map = filename.replace(extension, ".namemap.json")
        try:
            self._load_id_to_name_map(filename_map)
        except FileNotFoundError:
            if not os.path.isfile(filename_map):
                LOGGER.warning(
                    f"""{filename_map} not found. Set 'id_to_name' map manually with
                               'mesh._load_id_to_name_map(filename)'."""
                )
            else:
                LOGGER.error(
                    f"""Failed to read surface/volume ID to name map from {filename_map}.
                    Set 'id_to_name' map manually with
                    'mesh._load_id_to_name_map(filename)'."""
                )
        return

    def _save_id_to_name_map(self, filename: Union[str, pathlib.Path]):
        """Save the ID to name map.

        Parameters
        ----------
        filename : Union[str, pathlib.Path]
            Path to file.
        """
        id_to_name = {
            "_surface_id_to_name": self._surface_id_to_name,
            "_volume_id_to_name": self._volume_id_to_name,
        }
        with open(filename, "w") as f:
            json.dump(id_to_name, f, indent=4)

    def _load_id_to_name_map(self, filename: Union[str, pathlib.Path]):
        """Load the ID to name map for volumes and surfaces.

        Parameters
        ----------
        filename : Union[str, pathlib.Path]
            Filename of the ID to the name map (JSON).
        """
        with open(filename, "r") as f:
            data = json.load(
                f,
                object_hook=lambda d: {
                    int(k) if k.lstrip("-").isdigit() else k: v for k, v in d.items()
                },
            )
            self._surface_id_to_name = data["_surface_id_to_name"]
            self._volume_id_to_name = data["_volume_id_to_name"]

        # check whether map is valid, and print info to logger.
        self.validate_ids_to_name_map()
        return

    def validate_ids_to_name_map(self):
        """Check whether there are any duplicate or unmapped surfaces/volumes."""
        # TODO: Ensure there are no duplicate names.
        unmapped_volumes = self._get_unmapped_volumes()
        unmapped_surfaces = self._get_unmapped_surfaces()
        unmapped_lines = self._get_unmapped_lines()

        duplicate_volume_names = self._get_duplicate_volume_names()
        duplicate_surface_names = self._get_duplicate_surface_names()
        duplicate_line_names = self._get_duplicate_line_names()

        if len(unmapped_volumes) > 0 or len(unmapped_surfaces) > 0 or len(unmapped_lines) > 0:
            LOGGER.debug(f"Volume IDs {unmapped_volumes} are not associated with a volume name.")
            LOGGER.debug(f"Surface IDs {unmapped_surfaces} are not associated with a surface name.")
            LOGGER.debug(f"Line IDs {unmapped_lines} are not associated with a surface name.")
            return False
        if (
            len(duplicate_surface_names) > 0
            or len(duplicate_volume_names) > 0
            or len(duplicate_line_names) > 0
        ):
            LOGGER.debug(f"Volume names {duplicate_volume_names} occur more than once.")
            LOGGER.debug(f"Surface names {duplicate_surface_names} occur more than once.")
            LOGGER.debug(f"Line names {duplicate_line_names} occur more than once.")
            return False
        else:
            return True

    def clean(self, ignore_nans_in_point_average: bool = False, **kwargs):
        """Merge duplicate points and return a cleaned copy.

        Parameters
        ----------
        ignore_nans_in_point_average : bool, default: False
            Whether to ignore nan values when averaging point data.

        Returns
        -------
        Mesh
            Cleaned copy of self.
        """
        self_c = copy.deepcopy(self)

        # Compute point data average ignoring nan values.
        if ignore_nans_in_point_average:
            if "produce_merge_map" not in list(kwargs.keys()):
                kwargs["produce_merge_map"] = True

            super(Mesh, self_c).__init__(pv.UnstructuredGrid(self).clean(**kwargs))

            merge_map = self_c["PointMergeMap"]
            for key, data in self.point_data.items():
                non_nan_avg = [
                    np.nanmean(data[merge_map == merge_id]) for merge_id in np.unique(merge_map)
                ]
                self_c.point_data[key] = non_nan_avg
        else:
            super(Mesh, self_c).__init__(pv.UnstructuredGrid(self).clean(**kwargs))

        return self_c

    def add_volume(self, volume: pv.UnstructuredGrid, id: int = None, name: str = None):
        """Add a volume.

        Parameters
        ----------
        volume : pv.PolyData
            PolyData representation of the volume to add.
        id : int
            ID of the volume to add. This ID is tracked as ``_volume-id``.
        name : str, default: None
            Name of the added volume. The added volume is not tracked by default.
        """
        if not id:
            if "_volume-id" not in volume.cell_data.keys():
                LOGGER.debug("Failed to set '_volume-id'.")
                return None
        else:
            if not isinstance(id, int):
                LOGGER.debug("'sid' should be an integer.")
                return None
            volume.cell_data["_volume-id"] = np.ones(volume.n_cells, dtype=float) * id

        if name:
            self._volume_id_to_name[id] = name

        self_copy = self._add_mesh(volume, keep_data=True, fill_float=np.nan)
        return self_copy

    def add_surface(
        self,
        surface: pv.PolyData,
        id: int = None,
        name: str = None,
        overwrite_existing: bool = False,
    ):
        """Add a surface.

        Parameters
        ----------
        surface : pv.PolyData
            PolyData representation of the surface to add.
        sid : int
            ID of the surface to add. This ID is tracked as ``_surface-id``.
        name : str, default: None
            Name of the added surface. The added surface is not tracked by default.
        overwrite_existing : bool, default: False
            Whether to overwrite a surface with the same ID. If ``False``, the added
            surface is appended.
        """
        if not id:
            if "_surface-id" not in surface.cell_data.keys():
                LOGGER.error("Failed to set '_surface-id'.")
                return None
        else:
            if not isinstance(id, int):
                LOGGER.error("'sid' should be an integer.")
                return None
            surface.cell_data["_surface-id"] = np.ones(surface.n_cells, dtype=float) * id

        if not overwrite_existing:
            if id in self.surface_ids:
                LOGGER.error(f"{id} is already used. Pick any ID other than {self.surface_ids}.")
                return None

        self_copy = self._add_mesh(surface, keep_data=True, fill_float=np.nan)

        if name:
            self._surface_id_to_name[id] = name

        return self_copy

    def add_lines(self, lines: pv.PolyData, id: int = None, name: str = None):
        """Add lines.

        Parameters
        ----------
        lines : pv.PolyData
            PolyData representation of the lines to add.
        id : int
            ID of the surface to add. This ID is tracked as ``_line-id``.
        name : str, optional
            Name of the added lines. The added lines are not tracked by default.
        """
        if not id:
            if "_line-id" not in lines.cell_data.keys():
                LOGGER.error("Failed to set '_line-id'")
                return None
        else:
            if not isinstance(id, int):
                LOGGER.error("'sid' should be an integer.")
                return None
            lines.cell_data["_line-id"] = np.ones(lines.n_cells, dtype=float) * id

        self_copy = self._add_mesh(lines, keep_data=True, fill_float=np.nan, merge_points=False)

        if name:
            self._line_id_to_name[id] = name

        return self_copy

    def get_volume(self, sid: int) -> pv.UnstructuredGrid:
        """Get a volume as a PyVista unstructured grid object."""
        return self._get_submesh(sid, scalar="_volume-id")

    def get_volume_by_name(self, name: str) -> pv.UnstructuredGrid:
        """Get the surface associated with a given name."""
        if name not in list(self._volume_name_to_id.keys()):
            LOGGER.error(f"No volume is associated with {name}.")
            return None
        volume_id = self._volume_name_to_id[name]
        return self.get_volume(volume_id)

    def get_surface(self, sid: int) -> Union[pv.PolyData, SurfaceMesh]:
        # ?: Return SurfaceMesh instead of PolyData?
        """Get a surface as a PyVista polydata object.

        Notes
        -----
        This method tries to return a ``SurfaceMesh`` object that also contains a name, ID,
        and additional convenience properties.
        """
        if sid in list(self._surface_id_to_name.keys()):
            return SurfaceMesh(
                self._get_submesh(sid, scalar="_surface-id").extract_surface(),
                name=self._surface_id_to_name[sid],
                id=sid,
            )
        else:
            return self._get_submesh(sid, scalar="_surface-id").extract_surface()

    def get_surface_by_name(self, name: str) -> Union[pv.PolyData, SurfaceMesh]:
        # ?: Return SurfaceMesh instead of PolyData?
        """Get the surface associated with a given name."""
        if name not in list(self._surface_name_to_id.keys()):
            LOGGER.error(f"No surface is associated with {name}.")
            return None
        surface_id = self._surface_name_to_id[name]
        return self.get_surface(surface_id)

    def get_lines(self, sid: int) -> pv.PolyData:
        """Get lines as a PyVista polydata object."""
        return self._get_submesh(sid, scalar="_line-id").extract_surface()

    def get_lines_by_name(self, name: str) -> pv.PolyData:
        """Get the lines associated with a given name."""
        if name not in list(self._line_name_to_id.keys()):
            LOGGER.error(f"No lines associated with {name}")
            return None
        line_id = self._line_name_to_id[name]
        return self.get_lines(line_id)

    def remove_surface(self, sid: int):
        """Remove a surface with a given ID.

        Parameters
        ----------
        sid : int
            ID of the surface to remove.
        """
        mask = self.cell_data["_surface-id"] == sid
        return self.remove_cells(mask, inplace=True)

    def remove_volume(self, vid: int):
        """Remove a volume with a given ID.

        Parameters
        ----------
        vid : int
            ID of the volume to remove.
        """
        mask = self.cell_data["_volume-id"] == vid
        return self.remove_cells(mask, inplace=True)

    def remove_lines(self, lid: int):
        """Remove a set of lines with a given ID.

        Parameters
        ----------
        lid : int
            ID of the lines to remove.
        """
        mask = self.cell_data["_line-id"] == lid
        return self.remove_cells(mask, inplace=True)

    @staticmethod
    def _get_shifted_id(solid_mesh: Mesh, conduction_mesh: Mesh) -> np.ndarray:
        """Get the shifted IDs of the conduction mesh.

        Parameters
        ----------
        solid_mesh : Mesh
            Solid mesh.
        conduction_mesh : Mesh
            Path mesh with "_is-connected" cell data.

        Returns
        -------
        np.ndarray
            Shifted node IDs of the conduction mesh.
        """
        from scipy import spatial

        kdtree = spatial.cKDTree(solid_mesh.points)

        is_connected = conduction_mesh["_is-connected"].astype(bool)
        querry_points = conduction_mesh.points[is_connected]
        dst, solid_id = kdtree.query(querry_points)
        LOGGER.info(f"Maximal distance from solid-beam connected node:{np.max(dst)}")

        shifted_ids = np.linspace(
            0, conduction_mesh.n_points - 1, num=conduction_mesh.n_points, dtype=int
        )
        # for connected nodes, replace by solid mesh ID
        shifted_ids[is_connected] = solid_id
        # for beam-only nodes, shift their IDs
        for i in range(conduction_mesh.n_points):
            if not is_connected[i]:
                shifted_ids[i] += solid_mesh.n_points - np.sum(is_connected[:i])

        return shifted_ids

    @staticmethod
    def _safe_line_merge(base: Mesh, add_mesh: Mesh, merge_id: list, target_id: list) -> Mesh:
        """Safely merge two line meshes by specify the node ID to be merged.

        Parameters
        ----------
        base : Mesh
            Base mesh to merge into.
        add_mesh : Mesh
            Mesh to add to merge into the base mesh.
        merged_id : list
            List of node IDs to be merged.
        target_id : list
            List of node IDs in the base to be merged.

        Returns
        -------
        Mesh
            Merged line mesh.

        Notes
        -----
            point_data["_is-connected"] and cell_data["_line-id"] will be merged.
        """

        def get_lines(m: pv.UnstructuredGrid | pv.PolyData):
            if m.GetNumberOfCells() == 0:
                return np.empty(shape=(0, 2), dtype=np.int_)
            if isinstance(m, pv.UnstructuredGrid):
                return m.cells.reshape(-1, 3)[:, 1:]
            elif isinstance(m, pv.PolyData):
                return m.lines.reshape(-1, 3)[:, 1:]

        base_points = base.points
        base_lines = get_lines(base)

        if base.GetNumberOfCells() == 0:
            point_data = np.empty(shape=(0,))
            cell_data = np.empty(shape=(0,))
        else:
            point_data = base.point_data["_is-connected"]
            cell_data = base.cell_data["_line-id"]

        if merge_id == []:
            # no merge
            new_points = add_mesh.points
            new_point_data = add_mesh.point_data["_is-connected"]
            new_lines = get_lines(add_mesh) + len(base_points)

        elif merge_id == [0]:
            new_points = add_mesh.points[1:]
            # first node is merged, lead to an offset of all lines
            new_lines = get_lines(add_mesh) + len(base_points) - 1
            # replace first node
            new_lines[0, 0] = target_id[0]

            # point data
            new_point_data = add_mesh.point_data["_is-connected"][1:]

        elif merge_id == [0, -1]:
            # first node is merged, lead to an offset of all lines
            # last node is just dropped
            new_points = add_mesh.points[1:-1]
            new_lines = get_lines(add_mesh) + len(base_points) - 1
            # replace first node
            new_lines[0, 0] = target_id[0]
            # replace last node
            new_lines[-1, 1] = target_id[1]

            # point data
            new_point_data = add_mesh.point_data["_is-connected"][1:-1]
        else:
            NotImplementedError("Do not handle this merge lines.")

        merged = pv.PolyData()
        merged.points = np.vstack((base_points, new_points))
        merged_lines = np.vstack((base_lines, new_lines))
        merged.lines = np.hstack(
            (2 * np.ones(len(merged_lines), dtype=int)[:, np.newaxis], merged_lines)
        )
        merged.cell_data["_line-id"] = np.hstack((cell_data, add_mesh.cell_data["_line-id"]))
        merged.point_data["_is-connected"] = np.hstack((point_data, new_point_data))

        return Mesh(merged)


@deprecated(
    reason="""Importing Part class with ``from ansys.health.heart.objects import Part``
    is deprecated. Import with ``from ansys.health.heart.parts import Part`` instead.""",
)
class Part:
    """Part class for backward compatibility."""

    def __init__(self, *args, **kwargs):
        from ansys.health.heart.parts import Part as RealPart

        self.__class__ = RealPart
        self.__init__(*args, **kwargs)
