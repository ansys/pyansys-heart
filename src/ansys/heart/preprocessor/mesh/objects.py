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

import copy
from enum import Enum
import json
import os
import pathlib
from typing import List, Literal, Union

import numpy as np

from ansys.heart.core import LOG as LOGGER
import ansys.heart.preprocessor.mesh.vtkmethods as vtkmethods
from ansys.heart.simulator.settings.material.material import MechanicalMaterialModel

try:
    import pyvista as pv
except ImportError:
    LOGGER.warning("Importing pyvista failed. Install with: pip install pyvista")

SURFACE_CELL_TYPES = [pv.CellType.QUAD, pv.CellType.TRIANGLE]
VOLUME_CELL_TYPES = [pv.CellType.HEXAHEDRON, pv.CellType.TETRA]


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
    """Get the global cell ids of a particular cell type.

    Parameters
    ----------
    mesh : pv.UnstructuredGrid
        Unstructured grid from which to obtain the global cell ids
    celltype : pv.CellType
        Cell type to get global cell ids of.

    Returns
    -------
    np.ndarray
        Array with global cell ids.
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


# TODO: Deprecate
class Feature:
    """Feature class."""

    def __init__(self, name: str = None) -> None:
        LOGGER.error(DeprecationWarning("Deprecated"))
        self.name = name
        """Name of feature."""
        self.type = None
        """Type of feature."""
        self.nsid: int = None
        """Node set id associated with feature."""
        self.pid: int = None
        """Part id associated with the feature."""

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
                    "Assigning less nodes than the original not implemented yet."
                )

        except:
            LOGGER.warning("Failed to set nodes.")
            return

    @property
    def triangles(self):
        """Triangular faces of the surface num_faces x 3."""
        faces = np.reshape(self.faces, (self.n_cells, 3 + 1))[:, 1:]
        return faces

    @triangles.setter
    def triangles(self, value: np.ndarray):
        # sets faces of PolyData
        try:
            num_faces = value.shape[0]
            faces = np.hstack([np.full((num_faces, 1), 3, dtype=np.int8), value])
            self.faces = faces
        except:
            return

    @property
    def boundary_edges(self):
        """Get boundary edges of self."""
        boundary_edges = vtkmethods.get_boundary_edge_loops(self, remove_open_edge_loops=False)
        boundary_edges = np.vstack(list(boundary_edges.values()))
        return boundary_edges

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
        # **********************

        self.name = name
        """Name of the surface."""

        self.id: int = id
        """ID of surface."""
        self.nsid: int = None
        """ID of corresponding set of nodes."""

        self.triangles = triangles
        """Triangular faces of the surface num_faces x 3."""
        self.nodes = nodes
        """Node coordinates."""

    @property
    def node_ids(self) -> np.ndarray:
        """Global node ids - sorted by earliest occurrence."""
        _, idx = np.unique(self.triangles.flatten(), return_index=True)
        node_ids = self.triangles.flatten()[np.sort(idx)]
        return node_ids

    @property
    def _boundary_nodes(self) -> np.ndarray:
        """Global node ids of nodes on the boundary of the mesh (if any)."""
        _, idx = np.unique(self.boundary_edges.flatten(), return_index=True)
        node_ids = self.boundary_edges.flatten()[np.sort(idx)]
        return node_ids

    def write_to_stl(self, filename: pathlib.Path = None) -> None:
        """Write the surface to a vtk file."""
        if not filename:
            filename = "_".join(self.name.lower().split()) + ".stl"
        if filename[-4:] != ".stl":
            filename = filename + ".stl"

        # NOTE: The below should yield the same stls, but somehow fluent meshing
        # produces a slightly different mesh. Should still be valid though
        # cleaned = self.clean()
        # cleaned.save(filename)
        # vtkmethods.add_solid_name_to_stl(filename, self.name, file_type="binary")

        vtk_surface = vtkmethods.create_vtk_surface_triangles(self.nodes, self.triangles)
        vtkmethods.vtk_surface_to_stl(vtk_surface, filename, self.name)
        return


class BeamMesh(pv.UnstructuredGrid, Feature):
    """Beam class."""

    all_beam_nodes = []
    # beam nodes array

    @property
    def nodes(self):
        """Node coordinates."""
        return np.array(self.points)

    @nodes.setter
    def nodes(self, array: np.ndarray):
        if isinstance(array, type(None)):
            return
        try:
            self.points = array
        except:
            LOGGER.warning("Failed to set nodes.")
            return

    @property
    def edges(self):
        """Tetrahedrons num_tetra x 4."""
        return self.cells_dict[pv.CellType.LINE]

    @edges.setter
    def edges(self, value: np.ndarray):
        # sets lines of UnstructuredGrid
        try:
            points = self.points
            celltypes = np.full(value.shape[0], pv.CellType.LINE, dtype=np.int8)
            lines = np.hstack([np.full(len(celltypes), 2)[:, None], value])
            super().__init__(lines, celltypes, points)
        except:
            LOGGER.warning("Failed to set lines.")
            return

    def __init__(
        self,
        name: str = None,
        edges: np.ndarray = None,
        nodes: np.ndarray = None,
        beam_nodes_mask: np.ndarray = None,
        pid: int = None,
        nsid: int = -1,
    ) -> None:
        super().__init__(self)
        Feature.__init__(self, name)

        self.edges = edges
        """Beams edges."""

        self.nodes = nodes
        """Node coordinates."""

        self.beam_nodes_mask = beam_nodes_mask
        """True for beam nodes, False for solid nodes."""

        self.pid = pid
        """Part id associated with the network."""

        self.nsid: int = nsid
        """Surface id associated with the network."""

        self._all_beam_nodes: np.ndarray = np.empty((0, 3))
        """Temporary attribute to save all previously created beam nodes."""


class Cavity(Feature):
    """Cavity class."""

    def __init__(self, surface: SurfaceMesh = None, centroid: np.ndarray = None, name=None) -> None:
        super().__init__(name)

        self.surface: SurfaceMesh = surface
        """Surface mesh making up the cavity."""
        self.centroid: np.ndarray = centroid
        """Centroid of the cavity."""

    @property
    def volume(self):
        """Volume of the cavity."""
        return self.surface.volume

    def compute_centroid(self):
        """Compute the centroid of the cavity."""
        # self.centroid = np.mean(self.surface.nodes[np.unique(self.surface.triangles), :], axis=0)
        self.centroid = self.surface.center
        return self.centroid


class Cap(Feature):
    """Cap class."""

    def __init__(self, name: str = None, node_ids: Union[List[int], np.ndarray] = []) -> None:
        super().__init__(name)
        self.node_ids = node_ids
        """(Global) node ids of the cap."""
        self.triangles = None
        """Triangulation of cap."""
        self.normal = None
        """Normal of cap."""
        self.centroid = None
        """Centroid of cap."""
        self.centroid_id = None
        """Centroid of cap ID (in case centroid node is created)."""
        return


class Point(Feature):
    """Point class. Can be used to collect relevant points in the mesh."""

    def __init__(self, name: str = None, xyz: np.ndarray = None, node_id: int = None) -> None:
        super().__init__(name)

        self.xyz: np.ndarray = xyz
        """XYZ Coordinates of point."""
        self.node_id: int = node_id
        """Global node id of point."""


class Mesh(pv.UnstructuredGrid):
    """Mesh class: inherits from pyvista UnstructuredGrid.

    Notes
    -----
    This class inherits from pyvista.UnstructuredGrid and adds additional
    attributes and convenience methods for enhanced functionality. E.g. we use _volume_id,
    _surface_id and _line_id cell arrays to keep track of "labeled" selections of
    cells. _volume_id is used to group 3D volume cells together.
    Any non 3D volume cell is labeled as numpy.nan. Similarly 2D and 1D cells are tracked
    through _surface_id and _line_id respectively.
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
        # TODO: manage cell data
        # TODO: could deprecate now that there is an add_volume method?
        try:
            points = self.points
            celltypes = np.full(value.shape[0], pv.CellType.TETRA, dtype=np.int8)
            tetra = np.hstack([np.full(len(celltypes), 4)[:, None], value])
            super().__init__(tetra, celltypes, points)
        except:
            LOGGER.warning("Failed to set tetrahedrons.")
            return

    @property
    def _surfaces(self) -> List[SurfaceMesh]:
        """List of boundaries in the mesh."""
        if self.surface_ids is None:
            return []
        return [self.get_surface(surface_id) for surface_id in self.surface_ids]

    @property
    def _volumes(self):
        """List of volumes in the mesh."""
        if self.volume_ids is None:
            return []
        return [self.get_volume(volume_id) for volume_id in self.volume_ids]

    @property
    def _global_triangle_ids(self):
        """Global ids of triangular cells."""
        return _get_global_cell_ids(self, pv.CellType.TRIANGLE)

    @property
    def _global_tetrahedron_ids(self):
        """Global ids of tetrahedral cells."""
        return _get_global_cell_ids(self, pv.CellType.TETRA)

    # TODO: deprecate
    @property
    def part_ids(self) -> np.ndarray:
        """Array of part ids indicating to which part the tetrahedron belongs.

        Notes
        -----
        This is derived from the "part-id" field in cell data
        """
        try:
            value = self.cell_data["tags"].astype(int)
            return value
        except (KeyError, NameError):
            LOGGER.warning("'tags' field not found in self.cell_data")
            value = None
        try:
            value = self.cell_data["part-id"].astype(int)
            return value
        except (KeyError, NameError):
            LOGGER.warning("'part-id' field not found in self.cell_data")
            value = None
        return value

    # TODO: This needs to be refactored
    @property
    def boundary_names(self) -> List[str]:
        """Iterate over boundaries and returns their names."""
        return [b.name for b in self.boundaries]

    # TODO: deprecate. This is redundant with the _add_mesh and
    # TODO corresponding add_surface, add_volume, add_line methods
    def _sync_nodes_of_surfaces(self):
        """Synchronize the node array of each associated surface.

        Notes
        -----
        Temporary until this module is refactored.
        """
        for b in self.boundaries:
            b.nodes = self.nodes

        return

    # TODO: This method needs to be refactored.
    def _get_surface_from_name(self, name: str = None):
        """Return a list of surfaces that match the given list of names.

        Notes
        -----
        Returns single surface. When multiple matches are found returns list of surfaces
        """
        surfaces_search = self.boundaries
        surfaces = [s for s in surfaces_search if s.name == name]
        if len(surfaces) == 0:
            return None
        if len(surfaces) == 1:
            return surfaces[0]
        else:
            return surfaces

    @property
    def surface_ids(self) -> np.ndarray:
        """Unique surface ids.

        Returns
        -------
        np.ndarray
            Array with unique surface ids
        """
        try:
            mask = np.isin(self.celltypes, SURFACE_CELL_TYPES)
            mask1 = np.invert(np.isnan(self.cell_data["_surface-id"]))
            mask = np.all(np.vstack((mask, mask1)), axis=0)
            return np.unique(self.cell_data["_surface-id"][mask])
        except KeyError:
            LOGGER.debug(f"Failed to extrect one of {SURFACE_CELL_TYPES}")
            return None

    @property
    def surface_names(self) -> List[str]:
        """List of surface names."""
        return [v for k, v in self._surface_id_to_name.items()]

    @property
    def volume_ids(self) -> np.ndarray:
        """Unique volume ids.

        Returns
        -------
        np.ndarray
            Array with unique volume ids
        """
        try:
            mask = np.isin(self.celltypes, VOLUME_CELL_TYPES)
            mask1 = np.invert(np.isnan(self.cell_data["_volume-id"]))
            mask = np.all(np.vstack((mask, mask1)), axis=0)
            return np.unique(self.cell_data["_volume-id"][mask])
        except KeyError:
            LOGGER.debug(f"Failed to extrect one of {VOLUME_CELL_TYPES}")
            return None

    @property
    def volume_names(self) -> List[str]:
        """List of volume names."""
        return [v for k, v in self._volume_id_to_name.items()]

    @property
    def line_ids(self) -> np.ndarray:
        """Unique line ids.

        Returns
        -------
        np.ndarray
            Array with unique line ids
        """
        try:
            mask = self.celltypes == pv.CellType.LINE
            mask1 = np.invert(np.isnan(self.cell_data["_line-id"]))
            mask = np.all(np.vstack((mask, mask1)), axis=0)
            return np.unique(self.cell_data["_line-id"][mask])
        except KeyError:
            return None

    @property
    def _surface_name_to_id(self):
        return _invert_dict(self._surface_id_to_name)

    @property
    def _volume_name_to_id(self):
        return _invert_dict(self._volume_id_to_name)

    @property
    def _global_cell_ids(self):
        """Global cell ids."""
        self._set_global_ids()
        return self.cell_data["_global-cell-ids"]

    @property
    def _global_point_ids(self):
        """Global point ids."""
        self._set_global_ids()
        return self.point_data["_global-point-ids"]

    def __init__(self, *args):
        super().__init__(*args)

        # ! replace this list by read-only property
        self.boundaries: List[SurfaceMesh] = []
        """List of boundary surface meshes within the part."""

        self._surface_id_to_name: dict = {}
        """Surface id to name map."""
        self._volume_id_to_name: dict = {}
        """Volume id to name map."""
        pass

    def _add_mesh(
        self,
        mesh_input: Union[pv.PolyData, pv.UnstructuredGrid],
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
        mesh_input : pv.PolyData | pv.UnstructuredGrid
            Mesh to add, either PolyData or UnstructuredGrid
        keep_data : bool, optional
            Flag specifying whether to try to keep mesh point/cell data, by default True
        """
        mesh = copy.copy(mesh_input)
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

        merged = pv.merge((self, mesh), merge_points=False)
        super().__init__(merged)
        return self

    def _set_global_ids(self):
        """Add global cell and point ids as cell and point data array."""
        self.cell_data["_global-cell-ids"] = np.array(np.arange(0, self.n_cells), dtype=int)
        self.point_data["_global-point-ids"] = np.array(np.arange(0, self.n_points), dtype=int)
        return

    def _get_submesh(
        self, sid: int, scalar: Literal["_surface-id", "_line-id", "_volume-id"]
    ) -> pv.UnstructuredGrid:
        # NOTE: extract_cells cleans the object, removing any unused points.
        if not scalar in self.cell_data.keys():
            LOGGER.debug(f"{scalar} does not exist in cell_data")
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

    def _get_unmapped_volumes(self):
        unmapped_ids = self.volume_ids[
            np.invert(np.isin(self.volume_ids, list(self._volume_id_to_name.keys())))
        ]
        return unmapped_ids

    def _get_unmapped_surfaces(self):
        unmapped_ids = self.surface_ids[
            np.invert(np.isin(self.surface_ids, list(self._surface_id_to_name.keys())))
        ]
        return unmapped_ids

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
        This tries to read a JSON file with the volume/surface id to name map
        with extension .namemap.json in the same directory as the file. Alternatively,
        you can read the name map manually by calling `._load_id_to_name_map()`

        Parameters
        ----------
        filename : Union[str, pathlib.Path]
            Path to filename.
        """
        super(Mesh, self).__init__(filename)
        extension = pathlib.Path(filename).suffix
        filename_map = filename.replace(extension, ".namemap.json")
        try:
            self._load_id_to_name_map(filename_map)
        except:
            if not os.path.isfile(filename_map):
                LOGGER.warning(f"{filename_map} not found.")
            else:
                LOGGER.error(f"Failed to read surface/volume id to name map from {filename_map}")
        return

    def _save_id_to_name_map(self, filename: Union[str, pathlib.Path]):
        """Save the id to name map.

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
        """Load the id to name map for volumes and surfaces.

        Parameters
        ----------
        filename : Union[str, pathlib.Path]
            Filename of the id to name map (JSON).
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
        # TODO Ensure there are no duplicate names.
        unmapped_volumes = self._get_unmapped_volumes()
        unmapped_surfaces = self._get_unmapped_surfaces()

        duplicate_volume_names = self._get_duplicate_volume_names()
        duplicate_surface_names = self._get_duplicate_surface_names()

        if len(unmapped_volumes) > 0 or len(unmapped_surfaces) > 0:
            LOGGER.debug(f"Volume ids {unmapped_volumes} not associated with a volume name.")
            LOGGER.debug(f"Surface ids {unmapped_surfaces} not associated with a surface name.")
            return False
        if len(duplicate_surface_names) > 0 or len(duplicate_volume_names) > 0:
            LOGGER.debug(f"Volume names {duplicate_volume_names} occur more than once")
            LOGGER.debug(f"Surface names {duplicate_surface_names} occur more than once")
            return False
        else:
            return True

    def clean(self, ignore_nans_in_point_average: bool = False, **kwargs):
        """Merge duplicate points and return cleaned copy.

        Parameters
        ----------
        ignore_nans_in_point_average : bool, optional
            Flag indicating whether to ignore nan values when averaging point data, by default False

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

    def add_volume(self, volume: pv.UnstructuredGrid, id: int = None):
        """Add a volume.

        Parameters
        ----------
        volume : pv.PolyData
            PolyData representation of the volume to add
        id : int
            ID of the volume to be added. This id will be tracked as "_volume-id"
        """
        if not id:
            if "_volume-id" not in volume.cell_data.keys():
                LOGGER.debug("Failed to set _volume-id")
                return None
        else:
            if not isinstance(id, int):
                LOGGER.debug("sid should by type int.")
                return None
            volume.cell_data["_volume-id"] = np.ones(volume.n_cells, dtype=float) * id

        self_copy = self._add_mesh(volume, keep_data=True, fill_float=np.nan)
        return self_copy

    def add_surface(self, surface: pv.PolyData, id: int = None):
        """Add a surface.

        Parameters
        ----------
        surface : pv.PolyData
            PolyData representation of the surface to add
        sid : int
            ID of the surface to be added. This id will be tracked as "_surface-id"
        """
        if not id:
            if "_surface-id" not in surface.cell_data.keys():
                LOGGER.debug("Failed to set _surface-id")
                return None
        else:
            if not isinstance(id, int):
                LOGGER.debug("sid should by type int.")
                return None
            surface.cell_data["_surface-id"] = np.ones(surface.n_cells, dtype=float) * id

        self_copy = self._add_mesh(surface, keep_data=True, fill_float=np.nan)
        return self_copy

    def add_lines(self, lines: pv.PolyData, id: int = None):
        """Add lines.

        Parameters
        ----------
        lines : pv.PolyData
            PolyData representation of the lines to add
        id : int
            ID of the surface to be added. This id will be tracked as "_line-id"
        """
        if not id:
            if "_line-id" not in lines.cell_data.keys():
                LOGGER.debug("Failed to set _surface-id")
                return None
        else:
            if not isinstance(id, int):
                LOGGER.debug("sid should by type int.")
                return None
            lines.cell_data["_line-id"] = np.ones(lines.n_cells, dtype=float) * id

        self_copy = self._add_mesh(lines, keep_data=True, fill_float=np.nan)
        return self_copy

    def get_volume(self, sid: int) -> pv.UnstructuredGrid:
        """Get a volume as a UnstructuredGrids object."""
        return self._get_submesh(sid, scalar="_volume-id")

    def get_volume_by_name(self, name: str) -> pv.UnstructuredGrid:
        """Get the surface associated with `name`."""
        if name not in list(self._volume_name_to_id.keys()):
            LOGGER.debug(f"No volume associated with {name}")
            return None
        volume_id = self._volume_name_to_id[name]
        return self.get_volume(volume_id)

    def get_surface(self, sid: int) -> pv.PolyData:
        # ?: Return SurfaceMesh instead of PolyData?
        """Get a surface as PolyData object."""
        return self._get_submesh(sid, scalar="_surface-id").extract_surface()

    def get_surface_by_name(self, name: str) -> pv.PolyData:
        # ?: Return SurfaceMesh instead of PolyData?
        """Get the surface associated with `name`."""
        if name not in list(self._surface_name_to_id.keys()):
            LOGGER.debug(f"No surface associated with {name}")
            return None
        surface_id = self._surface_name_to_id[name]
        return self.get_surface(surface_id)

    def get_lines(self, sid: int) -> pv.PolyData:
        """Get lines as a PolyData object."""
        return self._get_submesh(sid, scalar="_line-id").extract_surface()

    def remove_surface(self, sid: int):
        """Remove a surface with id.

        Parameters
        ----------
        sid : int
            Id of surface to remove.
        """
        mask = self.cell_data["_surface-id"] == sid
        return self.remove_cells(mask, inplace=True)

    def remove_volume(self, vid: int):
        """Remove a volume with id.

        Parameters
        ----------
        vid : int
            Id of volume to remove.
        """
        mask = self.cell_data["_volume-id"] == vid
        return self.remove_cells(mask, inplace=True)

    def remove_lines(self, lid: int):
        """Remove a set of lines with id.

        Parameters
        ----------
        lid : int
            Id of lines to remove.
        """
        mask = self.cell_data["_volume-id"] == lid
        return self.remove_cells(mask, inplace=True)


class PartType(Enum):
    """Stores valid part types."""

    VENTRICLE = "ventricle"
    ATRIUM = "atrium"
    SEPTUM = "septum"
    ARTERY = "artery"
    MYOCARDIUM = "myocardium"
    UNDEFINED = "undefined"


class Part:
    """Part class."""

    @property
    def surfaces(self) -> List[SurfaceMesh]:
        """List of surfaces belonging to part."""
        surfaces = []
        for key, value in self.__dict__.items():
            if isinstance(value, SurfaceMesh):
                surfaces.append(value)
        return surfaces

    @property
    def surface_names(self) -> List[str]:
        """List of surface names belonging to part."""
        surface_names = []
        for key, value in self.__dict__.items():
            if isinstance(value, SurfaceMesh):
                surface_names.append(value.name)
        return surface_names

    def get_point(self, pointname: str) -> Point:
        """Get point from part."""
        for point in self.points:
            if point.name == pointname:
                return point
        LOGGER.error("Cannot find point {0:s}.".format(pointname))
        return None

    def __init__(self, name: str = None, part_type: PartType = PartType.UNDEFINED) -> None:
        self.name = name
        """Name of the part."""
        self.pid = None
        """Part ID."""
        self.mid = None
        """Material id associated with part."""
        self.part_type: PartType = part_type
        """Type of the part."""
        self.element_ids: np.ndarray = np.empty((0, 4), dtype=int)
        """Array holding element ids that make up this part."""
        self.points: List[Point] = []
        """Points of interest belonging to the part."""
        self.caps: List[Cap] = []
        """List of caps belonging to the part."""
        self.cavity: Cavity = None

        self.fiber: bool = False
        """If this part has fiber/sheet data."""
        self.active: bool = False
        """If active stress will be established."""

        self.meca_material: MechanicalMaterialModel = MechanicalMaterialModel.DummyMaterial()
        """Material model will be assiggned in Simulator."""

        """Cavity belonging to the part."""
        if self.part_type in [PartType.VENTRICLE]:
            self.apex_points: List[Point] = []
            """Points on apex."""

        self._add_surfaces()

    def _add_surfaces(self):
        """Add surfaces to the part."""
        if self.part_type in [PartType.VENTRICLE, PartType.ATRIUM]:
            self.endocardium = SurfaceMesh(name="{0} endocardium".format(self.name))
            """Endocardium."""
            self.epicardium = SurfaceMesh(name="{0} epicardium".format(self.name))
            """Epicardium."""
            if self.part_type == PartType.VENTRICLE:
                self.septum = SurfaceMesh(name="{0} septum".format(self.name))
                """Septum surface."""
        elif self.part_type in [PartType.ARTERY]:
            self.wall = SurfaceMesh(name="{0} wall".format(self.name))
            """Wall."""
        return

    def _add_myocardium_part(self):
        self.myocardium = Part(name="myocardium", part_type=PartType.MYOCARDIUM)
        return

    def _add_septum_part(self):
        self.septum = Part(name="septum", part_type=PartType.SEPTUM)
        return
