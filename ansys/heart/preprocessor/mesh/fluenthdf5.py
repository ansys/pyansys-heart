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

        Note
        ----
        Requires list of all cells.

        """
        self.cells = all_cells[self.min_id : self.max_id]
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
        return [cz.name for cz in self.cell_zones if cz != None]

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
        self.cell_zones: List[FluentCellZone] = None
        """List of cell zones."""
        self.face_zones: List[FluentFaceZone] = None
        """List of face zones."""

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
        self._read_cell_zone_info()

        self._read_all_faces_of_face_zones()

        if reconstruct_tetrahedrons:
            self._read_c0c1_of_face_zones()
            self._convert_interior_faces_to_tetrahedrons2()
            self._set_cells_in_cell_zones()

        self._update_indexing()

        self._close_file()
        return

    def _update_indexing(self) -> None:
        """Update the indexing of all cells and faces: index should start from 0."""
        for cell_zone in self.cell_zones:
            cell_zone.cells = cell_zone.cells - 1
        for face_zone in self.face_zones:
            face_zone.faces = face_zone.faces - 1
        return

    def _set_cells_in_cell_zones(self) -> List[FluentCellZone]:
        """Iterate over the cell zones and assigns cells to them."""
        for cell_zone in self.cell_zones:
            zone_cell_ids = np.arange(cell_zone.min_id, cell_zone.max_id, 1)
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
        num_cell_zones = len(cell_zone_ids)

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

        Do this for all interior faces.

        """
        self.cells = np.zeros((0, 4), dtype=int)
        self.cell_ids = np.zeros(0, dtype=int)

        for face_zone in self.face_zones:
            if "interior" not in face_zone.name:
                continue

            faces = face_zone.faces
            c0c1 = face_zone.c0c1.T.ravel()

            cell_ids1, idx1, counts1 = np.unique(c0c1, return_index=True, return_counts=True)
            cell_ids2, idx2, counts2 = np.unique(
                np.flipud(c0c1), return_index=True, return_counts=True
            )

            faces_temp = np.vstack([faces, faces])

            f1 = faces_temp[idx1, :]
            f2 = np.flipud(faces_temp)[idx2, :]

            # Find node in connected face which completes tetrahedron
            mask = (f2[:, :, None] == f1[:, None, :]).any(-1)
            mask = np.invert(mask)

            if not np.all(np.sum(mask, axis=1) == 1):
                # note: in an edge case we may have an interior face that
                # is connected to a cell where all the other faces are actually
                # not of the interior.
                # as a solution we can collect all faces first and consequently construct
                # the tetrahedrons
                face_ids = np.where(np.sum(mask, axis=1) != 1)[0]
                for face_id in face_ids:
                    f1[face_id]

                raise ValueError("The two faces do not seem to be connected with two nodes")

            tetrahedrons = np.hstack([f1, f2[mask][:, None]])

            # element ids?
            cell_ids = cell_ids1

            self.cells = np.vstack([self.cells, tetrahedrons])
            self.cell_ids = np.append(self.cell_ids, cell_ids)

        return tetrahedrons, cell_ids

    def _convert_interior_faces_to_tetrahedrons2(self) -> Tuple[np.ndarray, np.ndarray]:
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


if __name__ == "__main__":
    print("protected")
