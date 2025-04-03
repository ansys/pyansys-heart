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

import copy
from enum import Enum
import json
import os
import pathlib
from typing import List, Literal, Union

from deprecated import deprecated
import numpy as np
import pyvista as pv

from ansys.heart.core import LOG as LOGGER
import ansys.heart.core.utils.vtk_utils as vtk_utils
from ansys.heart.simulator.settings.material.ep_material import EPMaterial
from ansys.heart.simulator.settings.material.material import MechanicalMaterialModel

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
        #! This class can be deprecated.
        self.name = name
        """Name of feature."""
        self.type = None
        """Type of feature."""
        self._node_set_id: int = None
        """Node set id associated with feature."""
        self._seg_set_id: int = None
        """Segment set id associated with feature."""
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
                    "Assigning less nodes than the original, not implemented yet."
                )

        except Exception as e:
            LOGGER.error(f"Failed to set nodes. {e}")
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
        except Exception:
            return

    @property
    def triangles_global(self):
        """Global triangle ids.

        Returns
        -------
        Tries to use point_data["_global-point-ids"] to retrieve
        triangle definitions in global ids.
        """
        return self.point_data["_global-point-ids"][self.triangles]

    @property
    def boundary_edges(self):
        """Get boundary edges of self."""
        boundary_edges = vtk_utils.get_boundary_edge_loops(self, remove_open_edge_loops=False)
        boundary_edges = np.vstack(list(boundary_edges.values()))
        return boundary_edges

    @property
    def boundary_edges_global(self):
        """Global point ids of boundary edges."""
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
        """ID of surface."""

        self.triangles = triangles
        """Triangular faces of the surface num_faces x 3."""
        self.nodes = nodes
        """Node coordinates."""
        self._seg_set_id: int = None
        """Segment set id."""
        self._node_set_id: int = None
        """Node set id."""

    @property
    def node_ids_triangles(self) -> np.ndarray:
        """Local node ids - sorted by earliest occurrence."""
        _, idx = np.unique(self.triangles.flatten(), return_index=True)
        node_ids = self.triangles.flatten()[np.sort(idx)]
        return node_ids

    @property
    def global_node_ids_triangles(self):
        """Retrieve the global node ids from point data."""
        return self.point_data["_global-point-ids"][self.node_ids_triangles]

    @property
    def _boundary_nodes(self) -> np.ndarray:
        """Global node ids of nodes on the boundary of the mesh (if any)."""
        _, idx = np.unique(self.boundary_edges.flatten(), return_index=True)
        node_ids = self.boundary_edges.flatten()[np.sort(idx)]
        return node_ids

    def force_normals_inwards(self):
        """Force the cell ordering of a the closed surface such that normals point inward."""
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


# TODO: Refactor BeamMesh: why is this different from "Mesh"?
@deprecated(reason="BeamMesh is replaced by new class.")
class _BeamMesh(pv.UnstructuredGrid, Feature):
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
        except Exception as e:
            LOGGER.error(f"Failed to set nodes. {e}")
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
        except Exception as e:
            LOGGER.error(f"Failed to set lines. {e}")
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

        self.pid = pid
        """Part id associated with the network."""

        self.nsid: int = nsid
        """Surface id associated with the network."""

        self._all_beam_nodes: np.ndarray = np.empty((0, 3))
        """Temporary attribute to save all previously created beam nodes."""

        self.ep_material: EPMaterial = EPMaterial.DummyMaterial()
        """Initialize dummy ep material model."""


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
        """Local node ids of cap edge."""
        edges = vtk_utils.get_boundary_edge_loops(self._mesh)
        edge_local_ids = np.unique(np.array([np.array(edge) for edge in edges.values()]))
        return edge_local_ids

    @property
    def global_node_ids_edge(self):
        """Global node ids of the edge of the cap."""
        return self._mesh.point_data["_global-point-ids"][self._local_node_ids_edge]

    @property
    def _local_centroid_id(self):
        """Local id of centroid."""
        centroid_id = np.setdiff1d(np.arange(0, self._mesh.n_points), self._local_node_ids_edge)
        if len(centroid_id) != 1:
            LOGGER.error("Failed to identify single centroid node.")
            return None

        return centroid_id[0]

    @property
    def global_centroid_id(self):
        """Global centroid id."""
        return self._mesh.point_data["_global-point-ids"][self._local_centroid_id]

    @property
    def centroid(self):
        """Centroid of cap."""
        return self._mesh.points[self._local_centroid_id, :]

    @property
    def cap_normal(self):
        """Compute mean normal of cap."""
        return np.mean(self._mesh.compute_normals().cell_data["Normals"], axis=0)

    def __init__(
        self,
        name: str = None,
        cap_type: CapType = None,
    ) -> None:
        super().__init__(name)
        """Centroid of cap ID (in case centroid node is created)."""
        self._mesh: SurfaceMesh = None

        if cap_type is None or isinstance(cap_type, CapType):
            self.type = cap_type
        else:
            LOGGER.warning(f"Failed to set cap type for {name}, {cap_type}")

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
                LOGGER.debug(f"Failed to give surface with id {sid} a name. {error}")
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
        """Global ids of triangular cells."""
        return _get_global_cell_ids(self, pv.CellType.TRIANGLE)

    @property
    def _global_tetrahedron_ids(self):
        """Global ids of tetrahedral cells."""
        return _get_global_cell_ids(self, pv.CellType.TETRA)

    @property
    def surface_ids(self) -> np.ndarray:
        """Unique surface ids.

        Returns
        -------
        np.ndarray
            Array with unique surface ids
        """
        try:
            mask = np.isin(self.celltypes, _SURFACE_CELL_TYPES)
            mask1 = np.invert(np.isnan(self.cell_data["_surface-id"]))
            mask = np.all(np.vstack((mask, mask1)), axis=0)
            return np.unique(self.cell_data["_surface-id"][mask])
        except KeyError:
            LOGGER.debug(f"Failed to extract one of {_SURFACE_CELL_TYPES}")
            return []

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
            mask = np.isin(self.celltypes, _VOLUME_CELL_TYPES)
            mask1 = np.invert(np.isnan(self.cell_data["_volume-id"]))
            mask = np.all(np.vstack((mask, mask1)), axis=0)
            return np.unique(self.cell_data["_volume-id"][mask])
        except KeyError:
            LOGGER.debug(f"Failed to extrect one of {_VOLUME_CELL_TYPES}")
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

        merged = pv.merge((self, mesh), merge_points=False, main_has_priority=False)
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
        if scalar not in self.cell_data.keys():
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
        you can read the name map manually by calling `._load_id_to_name_map(filename)`

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
        except FileNotFoundError:
            if not os.path.isfile(filename_map):
                LOGGER.warning(
                    f"""{filename_map} not found. Please set id_to_name map manually by
                               mesh._load_id_to_name_map(filename)"""
                )
            else:
                LOGGER.error(
                    f"""Failed to read surface/volume id to name map from {filename_map}.
                    Please set id_to_name map manually by
                    mesh._load_id_to_name_map(filename)"""
                )
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
        # TODO: Ensure there are no duplicate names.
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

    def add_volume(self, volume: pv.UnstructuredGrid, id: int = None, name: str = None):
        """Add a volume.

        Parameters
        ----------
        volume : pv.PolyData
            PolyData representation of the volume to add
        id : int
            ID of the volume to be added. This id will be tracked as "_volume-id"
        name : str, optional
            Name of the added volume, by default None (not tracked)
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
            PolyData representation of the surface to add
        sid : int
            ID of the surface to be added. This id will be tracked as "_surface-id"
        name : str, optional
            Name of the added surface, by default None (not tracked)
        overwrite_existing : bool, optional
            Flag indicating whether to overwrite/append a surface with the same id, by default False
        """
        if not id:
            if "_surface-id" not in surface.cell_data.keys():
                LOGGER.error("Failed to set _surface-id")
                return None
        else:
            if not isinstance(id, int):
                LOGGER.error("sid should by type int.")
                return None
            surface.cell_data["_surface-id"] = np.ones(surface.n_cells, dtype=float) * id

        if not overwrite_existing:
            if id in self.surface_ids:
                LOGGER.error(
                    f"{id} already used. Please pick any id other than {self.surface_ids}."
                )
                return None

        self_copy = self._add_mesh(surface, keep_data=True, fill_float=np.nan)

        if name:
            self._surface_id_to_name[id] = name

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
                LOGGER.error("Failed to set _surface-id")
                return None
        else:
            if not isinstance(id, int):
                LOGGER.error("sid should by type int.")
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
            LOGGER.error(f"No volume associated with {name}")
            return None
        volume_id = self._volume_name_to_id[name]
        return self.get_volume(volume_id)

    def get_surface(self, sid: int) -> Union[pv.PolyData, SurfaceMesh]:
        # ?: Return SurfaceMesh instead of PolyData?
        """Get a surface as PolyData object.

        Notes
        -----
        Tries to return a SurfaceMesh object that also contains a name and id.
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
        """Get the surface associated with `name`."""
        if name not in list(self._surface_name_to_id.keys()):
            LOGGER.error(f"No surface associated with {name}")
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


class _ConductionType(Enum):
    """Enum containing type of conduction system."""

    LEFT_PURKINJE = "Left-purkinje"
    """Left Purkinje network."""
    RIGHT_PURKINJE = "Right-purkinje"
    """Right Purkinje network."""
    SAN_AVN = "SAN_to_AVN"
    """Sino-atrial node to atrio-ventricular node."""
    LEFT_BUNDLE_BRANCH = "Left bundle branch"
    """Left bundle branch."""
    RIGHT_BUNDLE_BRANCH = "Right bundle branch"
    """Right bundle branch."""
    HIS = "His"
    """His segment."""
    BACHMANN_BUNDLE = "Bachman bundle"
    """Bachmann bundle."""


class _BeamsMesh(Mesh):
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

    def __init__(self, *args):
        super().__init__(*args)

        self._line_id_to_name: dict = {}
        """line id to name map."""
        self.ep_material: dict = {}
        """Ep material map."""
        self._line_id_to_pid: dict = {}
        """line id to part id map."""
        pass

    def _get_submesh(
        self, sid: int, scalar: Literal["_surface-id", "_line-id", "_volume-id"]
    ) -> pv.PolyData:
        # NOTE: extract_cells cleans the object, removing any unused points.
        if scalar not in self.cell_data.keys():
            LOGGER.debug(f"{scalar} does not exist in cell_data")
            return None
        mask = np.isin(self.cell_data[scalar], sid)
        self._set_global_ids()
        return self.extract_cells(mask)

    def _add_mesh(
        self,
        mesh_input: pv.PolyData,
        keep_data: bool = True,
        fill_float: np.float64 = np.nan,
        fill_int: int = 0,
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

        merged = pv.merge((self, mesh), merge_points=True, main_has_priority=False)
        super().__init__(merged)
        return self

    def get_unique_lines_id(self) -> int:
        """Get unique lines id."""
        new_id: int
        if "_line-id" not in self.cell_data.keys():
            new_id = 1
        else:
            new_id = np.max(np.unique(self.cell_data["_line-id"])) + 1
        return int(new_id)

    def add_lines(self, lines: pv.PolyData, id: int = None, name: str = None):
        """Add lines.

        Parameters
        ----------
        lines : pv.PolyData
            PolyData representation of the lines to add
        id : int
            ID of the surface to be added. This id will be tracked as "_line-id"
        """
        if not id:
            return None
        else:
            if not isinstance(id, int):
                LOGGER.debug("sid should by type int.")
                return None
            lines.cell_data["_line-id"] = np.ones(lines.n_cells, dtype=float) * id
            if "_is-connected" not in lines.point_data.keys():
                lines.point_data["_is-connected"] = np.zeros(lines.n_points, dtype=int)
        self_copy = self._add_mesh(lines, keep_data=True, fill_float=np.nan)
        if name:
            self._line_id_to_name[id] = name
            self.ep_material[id] = EPMaterial.DummyMaterial()
        return self_copy

    def get_line_id_from_name(self, name: str) -> int:
        """Get line id from name using the `_line_id_to_name` attribute."""
        position_in_list = list(self._line_id_to_name.values()).index(name)
        line_id = list(self._line_id_to_name.keys())[position_in_list]
        return line_id

    def get_lines_by_name(self, name: str) -> pv.PolyData:
        # ?: Return SurfaceMesh instead of PolyData?
        """Get the lines associated with `name`."""
        if name not in list(self._line_id_to_name.values()):
            LOGGER.error(f"No lines associated with {name}")
            return None
        line_id = self.get_line_id_from_name(name)
        return self.get_lines(line_id)

    def get_lines(self, sid: int) -> pv.PolyData:
        """Get lines as a PolyData object."""
        return self._get_submesh(sid, scalar="_line-id").extract_surface()


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

        self.ep_material: EPMaterial = EPMaterial.DummyMaterial()
        """EP Material model will be assiggned in Simulator."""

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
                self.septum = SurfaceMesh(name="{0} endocardium septum".format(self.name))
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

    def _get_info(self):
        """Get part info in order to reconstruct from a mesh file."""
        info = {
            self.name: {
                "part-id": self.pid,
                "part-type": self.part_type.value,
                "surfaces": {},
                "caps": {},
                "cavity": {},
            }
        }

        info2 = {}
        info2["surfaces"] = {}
        info2["caps"] = {}
        info2["cavity"] = {}
        for surface in self.surfaces:
            if isinstance(surface, SurfaceMesh):
                if surface.id:
                    info2["surfaces"][surface.name] = surface.id

        for cap in self.caps:
            info2["caps"][cap.name] = cap._mesh.id

        if self.cavity:
            info2["cavity"][self.cavity.surface.name] = self.cavity.surface.id

        info[self.name].update(info2)

        return info
