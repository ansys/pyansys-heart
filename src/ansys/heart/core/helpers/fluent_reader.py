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

from typing import List, Tuple

import h5py
import numpy as np


class FluentCellZone:
    """Class that stores information of the cell zone."""

    def __init__(
        self, min_id: int = None, max_id: int = None, name: str = None, cid: int = None
    ) -> None:
        self.min_id: int = min_id
        """Min cell id of the cell zone: indexing starts at 0."""
        self.max_id: int = max_id
        """Max cell id of the cell zone: indexing starts at 0."""
        self.name: str = name
        """Name of the cell zone."""
        self.id: int = cid
        """Id of the cell zone."""
        self.cells: np.ndarray = None
        """Array of cells for this cell zone."""

        return

    def get_cells(self, all_cells: np.ndarray) -> None:
        """Select the cells between min and max id.

        Notes
        -----
        Requires list of all cells.

        """
        self.cells = all_cells[self.min_id : self.max_id + 1, :]
        return


class FluentFaceZone:
    """Class that stores information of the face zone."""

    def __init__(
        self,
        min_id: int = None,
        max_id: int = None,
        name: str = None,
        zone_id: int = None,
        zone_type: str = None,
        hdf5_id: int = None,
        faces: np.ndarray = None,
        c0c1: np.ndarray = None,
    ) -> None:
        self.min_id: int = min_id
        """Min face id of the face zone: indexing starts at 0."""
        self.max_id: int = max_id
        """Max face id of the face zone: indexing starts at 0."""
        self.name: str = name
        """Name of the face zone."""
        self.id: int = zone_id
        """Id of the face zone."""
        self.zone_type: str = zone_type
        """Type of face zone."""
        self.faces: np.ndarray = faces
        """Array of faces for this face zone."""
        self.c0c1: np.ndarray = c0c1
        """Array that stores connected cell-ids."""
        self.hdf5_id = hdf5_id
        """Id of face zone in hdf5 file."""

        return


class FluentMesh:
    """Class that stores the Fluent mesh."""

    @property
    def cell_zone_names(self):
        """List of cell zone names of non-empty cell zones."""
        return [cz.name for cz in self.cell_zones if cz is not None]

    def __init__(self, filename: str = None) -> None:
        self.filename: str = filename
        """Path to file."""
        self.fid: h5py.File = None
        """File id to h5py file."""
        self.nodes: np.ndarray = None
        """All nodes of the mesh."""
        self.faces: np.ndarray = None
        """All faces."""
        self.cells: np.ndarray = None
        """All cells."""
        self.cell_ids: np.ndarray = None
        """Array of cell ids use to define the cell zones."""
        self.cell_zones: List[FluentCellZone] = []
        """List of cell zones."""
        self.face_zones: List[FluentFaceZone] = []
        """List of face zones."""
        self._unique_map: np.ndarray = None
        """Map to go from full node list to node-list without duplicates."""

        pass

    def load_mesh(self, filename: str = None, reconstruct_tetrahedrons: bool = True) -> None:
        """Load the mesh from the hdf5 file."""
        if not filename and not self.filename:
            raise FileNotFoundError("Please specify a file to read")

        if self.filename:
            filename = self.filename

        self._open_file(filename)

        self._read_nodes()
        self._read_face_zone_info()
        self._read_all_faces_of_face_zones()
        self._remove_duplicate_nodes()

        if reconstruct_tetrahedrons:
            self._read_cell_zone_info()
            self._read_c0c1_of_face_zones()
            self._convert_interior_faces_to_tetrahedrons()
            self._set_cells_in_cell_zones()
            self._remove_empty_cell_zones()

        self._close_file()
        return

    def clean(self):
        """Remove all unused nodes."""
        used_node_ids1 = np.unique(
            np.array(np.vstack([fz.faces for fz in self.face_zones]), dtype=int)
        )
        used_node_ids2 = np.unique(np.array(self.cells, dtype=int))
        used_node_ids = np.unique(np.append(used_node_ids1, used_node_ids2))
        mask = np.zeros((self.nodes.shape[0]), dtype=bool)
        mask[used_node_ids] = True

        old_to_new_indices = np.zeros(self.nodes.shape[0], dtype=int) - 1
        old_to_new_indices[used_node_ids] = np.arange(0, used_node_ids.shape[0])

        # remove unused nodes
        self.nodes = self.nodes[used_node_ids, :]
        # reorder cell and face zones accordingly
        self.cells = old_to_new_indices[self.cells]
        for cz in self.cell_zones:
            cz.cells = old_to_new_indices[cz.cells]
        for fz in self.face_zones:
            fz.faces = old_to_new_indices[fz.faces]

        return

    def _remove_duplicate_nodes(self) -> None:
        """Remove duplicate nodes and remaps the face zone definitions."""
        self._unique_nodes, _, self._unique_map = np.unique(
            self.nodes,
            axis=0,
            return_index=True,
            return_inverse=True,
        )
        for fz in self.face_zones:
            fz.faces = self._unique_map[fz.faces - 1]
        self.nodes = self._unique_nodes
        return

    def _set_cells_in_cell_zones(self) -> List[FluentCellZone]:
        """Iterate over the cell zones and assigns cells to them."""
        for cell_zone in self.cell_zones:
            zone_cell_ids = np.arange(cell_zone.min_id, cell_zone.max_id + 1, 1)
            mask = np.isin(self.cell_ids, zone_cell_ids)
            cell_zone.cells = self.cells[mask, :]

        return self.cell_zones

    def _open_file(self, filename: str = None) -> h5py.File:
        """Open the file for reading."""
        if not filename:
            raise ValueError("Please specify input file")

        if filename[-7:] != ".msh.h5":
            raise FileNotFoundError("File does not have extension '.msh.h5'")

        self.fid = h5py.File(filename, "r")
        return self.fid

    def _close_file(self) -> None:
        """Close file."""
        self.fid.close()
        return

    def _read_nodes(self) -> None:
        """Read the node field(s)."""
        self.nodes = np.zeros((0, 3), dtype=float)
        for ii in np.array(self.fid["meshes/1/nodes/coords"]):
            self.nodes = np.vstack([self.nodes, np.array(self.fid["meshes/1/nodes/coords/" + ii])])
        return

    def _read_cell_zone_info(self) -> List[FluentCellZone]:
        """Initialize the list of cell zones."""
        cell_zone_names = (
            np.chararray.tobytes(np.array(self.fid["meshes/1/cells/zoneTopology/name"]))
            .decode()
            .split(";")
        )
        cell_zone_ids = np.array(self.fid["meshes/1/cells/zoneTopology/id"], dtype=int)
        min_ids = np.array(self.fid["meshes/1/cells/zoneTopology/minId"], dtype=int)
        max_ids = np.array(self.fid["meshes/1/cells/zoneTopology/maxId"], dtype=int)
        cell_zones: List[FluentCellZone] = []

        for ii in range(0, len(cell_zone_names), 1):
            cell_zones.append(
                FluentCellZone(
                    name=cell_zone_names[ii],
                    cid=cell_zone_ids[ii],
                    min_id=min_ids[ii],
                    max_id=max_ids[ii],
                )
            )
        self.cell_zones = cell_zones
        return cell_zones

    def _read_face_zone_info(self) -> List[FluentFaceZone]:
        """Initialize the list of face zones."""
        ids = np.array(self.fid["meshes/1/faces/zoneTopology/id"], dtype=int)
        max_ids = np.array(self.fid["meshes/1/faces/zoneTopology/maxId"], dtype=int)
        min_ids = np.array(self.fid["meshes/1/faces/zoneTopology/minId"], dtype=int)
        names = (
            np.chararray.tobytes(np.array(self.fid["meshes/1/faces/zoneTopology/name"]))
            .decode()
            .split(";")
        )
        zone_types = np.array(self.fid["meshes/1/faces/zoneTopology/zoneType"], dtype=int)
        num_face_zones = len(ids)
        face_zones: List[FluentFaceZone] = []

        for ii in range(0, num_face_zones, 1):
            face_zones.append(
                FluentFaceZone(
                    min_id=min_ids[ii],
                    max_id=max_ids[ii],
                    name=names[ii],
                    zone_id=ids[ii],
                    zone_type=zone_types[ii],
                    hdf5_id=ii + 1,
                )
            )
        self.face_zones = face_zones
        return face_zones

    def _read_all_faces_of_face_zones(self) -> List[FluentFaceZone]:
        """Read the faces of the face zone."""
        for face_zone in self.face_zones:
            subdir = "meshes/1/faces/nodes/" + str(face_zone.hdf5_id) + "/nodes"
            subdir2 = "meshes/1/faces/nodes/" + str(face_zone.hdf5_id) + "/nnodes"
            nnodes = np.array(self.fid[subdir2], dtype=int)
            if not np.all(nnodes == 3):
                raise ValueError("Only triangular meshes supported")

            node_ids = np.array(self.fid[subdir], dtype=int)
            num_triangles = int(len(node_ids) / 3)
            face_zone.faces = np.reshape(node_ids, (num_triangles, 3))

        return self.face_zones

    def _read_c0c1_of_face_zones(self) -> List[FluentFaceZone]:
        """Read the cell connectivity of the face zone. Only do for interior cells."""
        for face_zone in self.face_zones:
            subdir0 = "meshes/1/faces/c0/" + str(face_zone.hdf5_id)
            subdir1 = "meshes/1/faces/c1/" + str(face_zone.hdf5_id)
            c0c1 = np.array([self.fid[subdir0], self.fid[subdir1]], dtype=int).T
            face_zone.c0c1 = c0c1

        return self.face_zones

    def _convert_interior_faces_to_tetrahedrons(self) -> Tuple[np.ndarray, np.ndarray]:
        """Use c0c1 matrix to get tetrahedrons.

        Notes
        -----
        f1: n1 n2 n3 c0 c1
        f2: n3 n1 n4 c0 c1

        If f1 and f2 are connected to same face - extract node not occurring in
        f1. The resulting four nodes will make up the tetrahedron

        Do this for all faces.

        """
        self.cells = np.zeros((0, 4), dtype=int)
        self.cell_ids = np.zeros(0, dtype=int)

        # collect all faces
        faces = np.empty((0, 3), dtype=int)
        c0c1 = np.empty((0, 2), dtype=int)
        for face_zone in self.face_zones:
            faces = np.vstack([faces, face_zone.faces])
            c0c1 = np.vstack([c0c1, face_zone.c0c1])

        c0c1 = c0c1.T.ravel()

        cell_ids1, idx1, counts1 = np.unique(c0c1, return_index=True, return_counts=True)
        cell_ids2, idx2, counts2 = np.unique(np.flipud(c0c1), return_index=True, return_counts=True)

        faces_temp = np.vstack([faces, faces])

        # remove the
        if cell_ids1[0] == 0:
            cell_ids1 = cell_ids1[1:]
            idx1 = idx1[1:]
            idx2 = idx2[1:]

        f1 = faces_temp[idx1, :]
        f2 = np.flipud(faces_temp)[idx2, :]

        # Find node in connected face which completes tetrahedron
        mask = (f2[:, :, None] == f1[:, None, :]).any(-1)
        mask = np.invert(mask)

        if not np.all(np.sum(mask, axis=1) == 1):
            raise ValueError("The two faces do not seem to be connected with two nodes")

        tetrahedrons = np.hstack([f1, f2[mask][:, None]])

        cell_ids = cell_ids1

        self.cells = np.vstack([self.cells, tetrahedrons])
        self.cell_ids = np.append(self.cell_ids, cell_ids)

        return tetrahedrons, cell_ids

    # NOTE: no typehint due to lazy import of pyvista
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
        try:
            import pyvista as pv
        except ImportError:
            print("Failed to import pyvista. Try installing pyvista with `pip install pyvista`.")

        if add_cells and add_faces:
            add_both = True
        else:
            add_both = False

        if add_cells:
            # get cell zone ids.
            cell_zone_ids = np.concatenate(
                [[cz.id] * cz.cells.shape[0] for cz in self.cell_zones], dtype=int
            )

            cells = np.empty((0, self.cell_zones[0].cells.shape[1]), dtype=int)
            for cz in self.cell_zones:
                cells = np.vstack([cells, cz.cells])

            num_cells = cells.shape[0]

            cells = np.hstack([np.ones((num_cells, 1), dtype=int) * 4, cells])
            celltypes = [pv.CellType.TETRA] * num_cells
            grid = pv.UnstructuredGrid(cells.flatten(), celltypes, self.nodes)

            grid.cell_data["cell-zone-ids"] = cell_zone_ids

        if add_faces:
            # add faces.
            grid_faces = pv.UnstructuredGrid()
            grid_faces.nodes = self.nodes

            face_zone_ids = np.concatenate([[fz.id] * fz.faces.shape[0] for fz in self.face_zones])
            faces = np.array(np.concatenate([fz.faces for fz in self.face_zones]), dtype=int)
            faces = np.hstack([np.ones((faces.shape[0], 1), dtype=int) * 3, faces])

            grid_faces = pv.UnstructuredGrid(
                faces.flatten(), [pv.CellType.TRIANGLE] * faces.shape[0], self.nodes
            )
            grid_faces.cell_data["face-zone-ids"] = face_zone_ids

        if add_both:
            # ensure same cell arrays are present but with dummy values.
            grid_faces.cell_data["cell-zone-ids"] = np.ones(grid_faces.n_cells, dtype=int) * -1
            grid.cell_data["face-zone-ids"] = np.ones(grid.n_cells, dtype=int) * -1

            return grid + grid_faces

        if add_faces:
            return grid_faces

        if add_cells:
            return grid

    def _fix_negative_cells(self):
        """Rorder base face in cells that have a negative cell volume.

        Notes
        -----
        For a positive volume the base face (n1, n2, n3) needs to point in direction
        of n4. Hence, swapping the order to (n3, n2, n1) fixes negative cell volumes.
        """
        grid = self._to_vtk(add_cells=True, add_faces=False)
        grid = grid.compute_cell_sizes(length=False, area=False, volume=True)
        negative_cells = np.argwhere(grid.cell_data["Volume"] <= 0).flatten()

        if negative_cells.shape[0] == 0:
            return

        reordered_cells = self.cells[:, [2, 1, 0, 3]]
        self.cells[negative_cells, :] = reordered_cells[negative_cells, :]

        # update cell zones with reordered cells
        self._set_cells_in_cell_zones()
        return

    def _remove_empty_cell_zones(self):
        """Remove empty cell zones from cell zone list."""
        self.cell_zones = [cz for cz in self.cell_zones if cz.cells.shape[0] > 0]
        return


if __name__ == "__main__":
    print("protected")
