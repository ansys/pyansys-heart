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

"""Module containing functions to read/write fluent meshes in HDF5 format."""

import numpy as np
import pyvista as pv

import ansys.dpf.core as dpf

num_points_to_volume_cell_type = {
    4: pv.CellType.TETRA,
    8: pv.CellType.HEXAHEDRON,
    5: pv.CellType.PYRAMID,
}


def _infer_cellzone_celltypes_from_cells(cells):
    """Simplify polyhedral cells into their simple 3D shape."""
    len_cells = len(cells)
    index = 0
    celltypes = []
    while index < len_cells:
        n_points_cell = cells[index]
        index = n_points_cell + index + 1
        celltypes += [num_points_to_volume_cell_type[n_points_cell]]
    return celltypes


def _get_cell_zones_from_meshes_container(
    meshes_container: dpf.MeshesContainer, mesh_info: dpf.MeshInfo
) -> list[pv.UnstructuredGrid]:
    """Get all cell zones as pyvista unstructured grid from the meshes container.

    Parameters
    ----------
    meshes_container : dpf.MeshesContainer
        Meshes container containing all meshes.

    Returns
    -------
    list[pv.UnstructuredGrid]
        List of unstructured grids.
    """
    cell_zones = []
    cell_zone_ids = [int(cell_zone_id) for cell_zone_id in mesh_info.cell_zones.keys()]
    available_ids = meshes_container.get_available_ids_for_label("zone")
    for mesh, zone_id in zip(meshes_container, available_ids):
        if zone_id in cell_zone_ids:
            # NOTE: volume elements are stored as pv.CellType.POLYHEDRON's, so convert to their
            # "simple" 3d shape. Currently just TETRA and HEXAHEDRON
            celltypes = _infer_cellzone_celltypes_from_cells(mesh.grid.cells)
            grid = pv.UnstructuredGrid(mesh.grid.cells, celltypes, mesh.grid.points)
            grid.cell_data["cell-zone-id"] = zone_id
            cell_zones += [grid]
    return cell_zones


def _get_grids_from_zones(
    meshes_container: dpf.MeshesContainer, mesh_info: dpf.MeshInfo
) -> list[pv.UnstructuredGrid]:
    """Get pyvista grids from each zone.

    Parameters
    ----------
    meshes_container : dpf.MeshesContainer
        Meshes container containing all meshes.
    mesh_info : dpf.MeshInfo
        Mesh info describing zone ids/zone names

    Returns
    -------
    list[pv.UnstructuredGrid]
        List of unstructured grids.
    """
    face_zone_ids = [int(face_zone_id) for face_zone_id in mesh_info.face_zones.keys()]
    cell_zone_ids = [int(cell_zone_id) for cell_zone_id in mesh_info.cell_zones.keys()]
    available_ids = meshes_container.get_available_ids_for_label("zone")
    cell_zones = []
    face_zones = []
    for mesh, zone_id in zip(meshes_container, available_ids):
        grid = mesh.grid
        if zone_id in cell_zone_ids:
            celltypes = _infer_cellzone_celltypes_from_cells(mesh.grid.cells)
            grid = pv.UnstructuredGrid(mesh.grid.cells, celltypes, mesh.grid.points)
            grid.cell_data["cell-zone-id"] = float(zone_id)
            grid.cell_data["face-zone-id"] = np.nan
            cell_zones += [grid]

        elif zone_id in face_zone_ids:
            grid.cell_data["face-zone-id"] = float(zone_id)
            grid.cell_data["cell-zone-id"] = np.nan
            face_zones += [grid]

    return cell_zones, face_zones


class _FluentMesh:
    """Class that stores the Fluent mesh."""

    @property
    def cell_zone_names(self):
        """List of cell zone names of non-empty cell zones."""
        return [key for key in self.meshinfo.cell_zones.values()]

    @property
    def face_zone_names(self):
        """List of cell zone names of non-empty cell zones."""
        return [val for val in self.meshinfo.face_zones.values()]

    @property
    def face_zone_ids(self):
        return [int(key) for key in self.meshinfo.face_zones.key()]

    @property
    def cell_zone_ids(self):
        return [int(key) for key in self.meshinfo.cell_zones.key()]

    def __init__(self, filename: str = None) -> None:
        self.filename: str = filename
        """Path to file."""

        self.data_source = dpf.DataSources()
        """DPF Datasource."""
        self.data_source.set_result_file_path(filename, key="cas")
        """DPF results path."""

        model = dpf.Model(self.data_source)
        self.meshinfo: dpf.MeshInfo = model.metadata.mesh_info
        """mesh info."""
        return

    def _load_mesh(self) -> None:
        """Load the mesh from the hdf5 file."""
        streams = dpf.operators.metadata.streams_provider(data_sources=self.data_source)
        # reads all meshes in one go.
        self._dpf_mesh = dpf.operators.mesh.meshes_provider(
            streams_container=streams, region_scoping=None
        ).eval()
        return

    def _to_vtk(self, add_cells: bool = True, add_faces: bool = False):
        """Convert mesh to vtk unstructured grid or polydata.

        Parameters
        ----------
        add_cells : bool, optional
            Whether to add cells to the vtk object, by default True
        add_faces : bool, optional
            Whether to add faces to the vtk object, by default False

        Returns
        -------
        pv.UnstructuredGrid
            Unstructured grid representation of the fluent mesh.
        """
        self._load_mesh()
        # cell_zones = _get_cell_zones_from_meshes_container(self._dpf_mesh, self.meshinfo)
        cell_zones, face_zones = _get_grids_from_zones(self._dpf_mesh, self.meshinfo)

        if add_cells and add_faces:
            return pv.merge(cell_zones + face_zones)

        if add_faces:
            return pv.merge(face_zones)

        if add_cells:
            return pv.merge(cell_zones)
