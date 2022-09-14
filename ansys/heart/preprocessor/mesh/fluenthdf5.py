import time
from typing import List, Optional

from ansys.heart.custom_logging import LOGGER
import h5py
import meshio
import numpy as np

"""Module containing functions to read/write fluent meshes in HDF5 format
"""


class FluentCellZone:
    """Class that stores information of the cell zone"""

    def __init__(
        self, min_id: int = None, max_id: int = None, name: str = None, cid: int = None
    ) -> None:

        self.min_id: int = min_id
        """Min cell id of the cell zone: indexing starts at 0"""
        self.max_id: int = max_id
        """Max cell id of the cell zone: indexing starts at 0"""
        self.name: str = name
        """Name of the cell zone"""
        self.id: int = cid
        """Id of the cell zone"""
        self.cells: np.ndarray = None
        """Array of cells for this cell zone"""

        return

    def get_cells(self, all_cells: np.ndarray) -> None:
        """Selects the cells between min and max id from the list
        of all cells"""
        self.cells = all_cells[self.min_id : self.max_id]
        return


class FluentFaceZone:
    """Class that stores information of the face zone"""

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
        """Min face id of the face zone: indexing starts at 0"""
        self.max_id: int = max_id
        """Max face id of the face zone: indexing starts at 0"""
        self.name: str = name
        """Name of the face zone"""
        self.id: int = zone_id
        """Id of the face zone"""
        self.zone_type: str = zone_type
        """Type of face zone"""
        self.faces: np.ndarray = faces
        """Array of faces for this face zone"""
        self.c0c1: np.ndarray = c0c1
        """Array that stores connected cell-ids"""
        self.hdf5_id = hdf5_id
        """Id of face zone in hdf5 file"""

        return


class FluentMesh:
    """Class that stores the Fluent mesh"""

    @property
    def cell_zone_names(self):
        """List of cell zone names of non-empty cell zones"""
        return [cz.name for cz in self.cell_zones if cz != None]

    def __init__(self, filename: str = None) -> None:

        self.filename: str = filename
        """Path to file"""
        self.fid: h5py.File = None
        """File id to h5py file"""
        self.nodes: np.ndarray = None
        """All nodes of the mesh"""
        self.faces: np.ndarray = None
        """All faces"""
        self.cells: np.ndarray = None
        """All cells"""
        self.cell_ids: np.ndarray = None
        """Array of cell ids use to define the cell zones"""
        self.cell_zones: List[FluentCellZone] = None
        """List of cell zones"""
        self.face_zones: List[FluentFaceZone] = None
        """List of face zones"""

        pass

    def load_mesh(self, filename: str = None, reconstruct_tetrahedrons: bool = True):
        """Loads the mesh from the hdf5 file"""
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

    def _update_indexing(self):
        """Updates the indexing of all cells and faces: indexshould start from 0"""
        for cell_zone in self.cell_zones:
            cell_zone.cells = cell_zone.cells - 1
        for face_zone in self.face_zones:
            face_zone.faces = face_zone.faces - 1
        return

    def _set_cells_in_cell_zones(self):
        """Iterates over the cell zones and assigns cells to them"""
        for cell_zone in self.cell_zones:
            zone_cell_ids = np.arange(cell_zone.min_id, cell_zone.max_id, 1)
            mask = np.isin(self.cell_ids, zone_cell_ids)
            cell_zone.cells = self.cells[mask, :]

        return self.cell_zones

    def _open_file(self, filename: str = None):
        """Opens the file for reading"""
        if not filename:
            raise ValueError("Please specify input file")

        if filename[-7:] != ".msh.h5":
            raise FileNotFoundError("File does not have extension '.msh.h5'")

        self.fid = h5py.File(filename, "r")
        return self.fid

    def _close_file(self):
        """Closes file"""
        self.fid.close()

    def _read_nodes(self):
        """Reads the node field(s)"""
        self.nodes = np.zeros((0, 3), dtype=float)
        for ii in np.array(self.fid["meshes/1/nodes/coords"]):
            self.nodes = np.vstack([self.nodes, np.array(self.fid["meshes/1/nodes/coords/" + ii])])

    def _read_cell_zone_info(self) -> List[FluentCellZone]:
        """Initializes the list of cell zones"""

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
        """Initializes the list of face zones"""

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

    def _read_all_faces_of_face_zones(self):
        """Reads the faces of the face zone"""

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

    def _read_c0c1_of_face_zones(self):
        """Reads the cell connectivity of the face zone. Only do for interior cells"""

        for face_zone in self.face_zones:
            # if face_zone.zone_type != 2 or "interior" not in face_zone.name:
            #     continue

            subdir0 = "meshes/1/faces/c0/" + str(face_zone.hdf5_id)
            subdir1 = "meshes/1/faces/c1/" + str(face_zone.hdf5_id)

            c0c1 = np.array([self.fid[subdir0], self.fid[subdir1]], dtype=int).T

            face_zone.c0c1 = c0c1

        return self.face_zones

    def _convert_interior_faces_to_tetrahedrons(self):
        """Uses c0c1 matrix to get tetrahedrons

        Notes
        -----

        f1: n1 n2 n3 c0 c1
        f2: n3 n1 n4 c0 c1

        If f1 and f2 are connected to same face - extract node not occuring in
        f1. The resulting four nodes will make up the tetrahedron

        Do this for all interior faces

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

    def _convert_interior_faces_to_tetrahedrons2(self):
        """Uses c0c1 matrix to get tetrahedrons

        Notes
        -----

        f1: n1 n2 n3 c0 c1
        f2: n3 n1 n4 c0 c1

        If f1 and f2 are connected to same face - extract node not occuring in
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


def _deprecated_fluenthdf5_to_vtk(hdf5_filename: str, vtk_filename: str):
    """Converts Fluent hdf5 to vtk format

    Parameters
    ----------
    hdf5filename : str
        path to hdf5 file
    vtkfilename : str
        path to vtk file
    """
    add_surface = False

    fid = h5py.File(hdf5_filename, "r")
    mesh_group = _deprecated_get_mesh_group(fid)
    num_cell_zones = len(mesh_group)

    if num_cell_zones > 1:
        raise NotImplementedError("Writing multiple cell zones not implemented yet")

    for key in mesh_group.keys():
        cell_zone_id = int(key)

        face_group = _deprecated_get_face_group(mesh_group[key])
        face_zone_info = _deprecated_get_face_zone_info(face_group)
        tetrahedrons = _deprecated_face_group_to_tetrahedrons(face_group, face_zone_info[0])
        face_zones = _deprecated_face_group_to_face_zones(
            face_group, face_zone_info[0], get_interior=False
        )

        points = _deprecated_get_nodes_from_mesh_group(mesh_group[key])[0]

        cells = []
        cell_data = []
        if add_surface:
            for _, value in face_zones.items():
                # logger.info(value)
                cells.append(("triangle", value["faces"] - 1))
                zoneids = np.ones(cells[-1][-1].shape[0], dtype=int) * value["zone-id"]
                cell_data.append(zoneids)

        cells.append(("tetra", tetrahedrons))

        cell_data.append(np.ones(tetrahedrons.shape[0] * cell_zone_id))
        cell_data = {"zone-id": cell_data}

        # add triangles: disable
        add_triangles = False  # this may cause some issues
        if add_triangles:
            # triangle_zones = ("triangle", [[0, 1, 2], [1, 3, 2]])
            triangles = np.empty((0, 3))
            zone_ids_all_triangles = np.empty(0)
            for zone, value in face_zones.items():
                num_faces = value["faces"].shape[0]
                zone_ids = np.ones(num_faces, dtype=int) * value["zone-id"]
                zone_ids_all_triangles = np.append(zone_ids_all_triangles, zone_ids)
                triangles = np.vstack([triangles, value["faces"]])

            # put in right format
            triangles = ("triangle", triangles)
            # append to cells
            cells.append(triangles)

            cell_data["zone-id"].append(zone_ids_all_triangles)

        # write file
        mesh = meshio.Mesh(
            points=points,
            cells=cells,
            # ("triangle", tris-1 ),
            # ("tetra", tetrahedrons - 1)
            # ],
            cell_data=cell_data,
        )
        mesh.write(vtk_filename)

    fid.close()

    return tetrahedrons, face_zones, points


def _deprecated_get_mesh_group(fid: h5py.File) -> h5py.Group:
    """Gets meshgroup from hdf5 file handle"""
    keys = list(fid.keys())
    if keys[0] != "meshes":
        raise KeyError("No 'meshes' field found")
    return fid["meshes"]


def _deprecated_get_face_group(group: h5py.Group) -> h5py.Group:
    """Gets face zone group"""

    def find_faces(name):
        if "faces" in name:
            return name

    face_group = group[group.visit(find_faces)]
    return face_group


def _deprecated_get_face_zone_info(face_group: h5py.Group) -> List:
    """Gets list of face zone names"""
    # get face zone ids
    face_zone_ids = np.array(face_group["zoneTopology/id"])
    # get face zones
    face_zone_names = np.array(face_group["zoneTopology/name"])
    face_zone_names = np.chararray.tobytes(face_zone_names).decode()
    face_zone_names = face_zone_names.split(";")

    face_zone_types = np.array(face_group["zoneTopology/zoneType"])
    return face_zone_names, face_zone_ids, face_zone_types


def _deprecated_get_nodes_from_mesh_group(mesh_group: h5py.Group):
    """Gets nodes from the mesh group"""

    nodes = np.empty((0, 3), dtype=float)
    for key in mesh_group["nodes/coords"].keys():
        node_group = np.array(mesh_group["nodes/coords/" + key])
        nodes = np.append(nodes, node_group, axis=0)

    node_ids = np.arange(0, np.shape(nodes)[0], 1, dtype=int)
    return nodes, node_ids


def _deprecated_face_group_to_face_zones(
    face_group: h5py.Group,
    face_zone_names: List[str],
    get_all_zones: bool = False,
    get_interior: Optional[bool] = False,
) -> dict:
    """extracts face zones from face group object

    Parameters
    ----------
    face_group : h5py.Group
        Face group object
    face_zone_names : List[str]
        List of face zone names to extract
    get_all_zones : bool, optional
        Flag to indicate whether to extract all face zones, by default False
    get_interior : Optional[bool], optional
        Flag to indicate whether to extract all interior face zones, by default False

    Returns
    -------
    dict
        dictionary with all face zone definitions
    """

    # get face definitions of all face zones
    faces_all = np.empty((0, 3))
    faces_zoneid = np.empty((0, 0), dtype=int)
    faces_label = []

    # checks whether provided list of face zones is present in the face group
    face_zone_names_1 = _deprecated_get_face_zone_info(face_group)[0]
    face_zone_ids = _deprecated_get_face_zone_info(face_group)[1]

    idx_face_zones = []
    for ii, face_zone in enumerate(face_zone_names_1):
        if face_zone in face_zone_names and "interior" not in face_zone:
            idx_face_zones.append(ii)

    faces_out = {}
    for ii in idx_face_zones:
        subdir1 = "nodes/" + str(ii + 1) + "/nodes"
        faces = np.array(face_group[subdir1])
        # concatenate all faces into single array
        faces = np.reshape(faces, (int(len(faces) / 3), 3))

        faces_zoneid = np.append(faces_zoneid, np.ones(faces.shape[0], dtype=int) * (ii + 1))

        faces_all = np.append(faces_all, faces, axis=0)

        num_faces = faces.shape[0]

        faces_label.extend([face_zone_names_1[ii]] * num_faces)
        faces_out[face_zone_names_1[ii]] = {
            "faces": faces - 1,
            "zone-id": face_zone_ids[ii],
        }

    return faces_out


def _deprecated_face_group_to_tetrahedrons(
    face_group: h5py.Group, face_zone_names: List[str]
) -> np.array:
    """ "Converts Fluent's face connectivity matrix to tetrahedron elements.

    Notes
    -----
    Format Fluent face:
    [n0 n1 n2] [c0 c1]
    n0 n1 n2 : node indices that make up face
    c0 c1    : cell indices to which face is connected. Max 2.

    Note: this finds n3 to construct a tetrahedron
    based on the connectivity of the faces. I.e.: the tetrahedron
    is defined by four faces as:
    f1: n1 n2 n3
    f2: n1 n2 n4
    f3: n2 n3 n4
    f3: n1 n3 n4

    Hence the tetrahedron will be composed of n1 n2 n3 n4.
    This function exploits this characteristic

    Parameters
    ----------
    face_group : h5py.Group
        handle to face group
    face_zone_names : List
        list of face zone names

    Returns
    -------
    np.array
        array with tetrahedron definitions

    Raises
    ------
    Exception
        Non-triangular faces detected
    """

    # get face definitions of all face zones
    faces_all = np.empty((0, 3))
    faces_zoneid = np.empty((0, 0), dtype=int)
    for ii, facezone in enumerate(face_zone_names):
        subdir1 = "nodes/" + str(ii + 1) + "/nodes"
        subdir2 = "nodes/" + str(ii + 1) + "/nnodes"
        faces = np.array(face_group[subdir1])
        nnodes = np.array(face_group[subdir2])

        if not np.all(nnodes == 3):
            raise Exception("Functionality for non-triangular faces not implemented...")

        # concatenate all faces into single array
        faces = np.reshape(faces, (int(len(faces) / 3), 3))
        faces_zoneid = np.append(faces_zoneid, np.ones(faces.shape[0]) * (ii + 1))
        faces_all = np.append(faces_all, faces, axis=0)

    # cell connectivity arrays
    c0c1_all = np.empty((0, 2), dtype=int)

    # get cell connectivity of all faces
    for ii, facezone in enumerate(face_zone_names):
        if "interior" not in facezone:
            continue
        subdir1 = "c0/" + str(ii + 1)
        subdir2 = "c1/" + str(ii + 1)
        c0c1 = np.column_stack((np.array(face_group[subdir1]), np.array(face_group[subdir2])))
        c0c1_all = np.append(c0c1_all, c0c1, axis=0)

    uniq_elem_id, uniq_count = np.unique(c0c1_all[c0c1_all > 0], return_counts=True)
    num_elem = uniq_elem_id.size

    tetrahedrons = np.empty((num_elem, 4), dtype=int)

    LOGGER.info("Generating {0} tetrahedrons...".format(num_elem))

    t0 = time.time()

    # exploit unique function for fast construction of tets
    face_idx = np.arange(0, c0c1_all.shape[0], 1, dtype=int)
    face_idx_r = np.hstack((face_idx, face_idx))
    face_idx_r_flip = np.flipud(face_idx_r)

    a1, b1 = np.unique(np.ravel(c0c1_all.T), return_index=True)

    a2, b2 = np.unique(np.flipud(np.ravel(c0c1_all.T)), return_index=True)
    a1 = np.delete(a1, 0)
    a2 = np.delete(a2, 0)
    b1 = np.delete(b1, 0)
    b2 = np.delete(b2, 0)

    # get node indices that construct all tetrahedrons
    A = faces_all[face_idx_r[b1], :]
    B = faces_all[face_idx_r_flip[b2], :]

    # find node index which to append to the reference triangle to complete tetrahedron
    C = (B[:, :, None] == A[:, None, :]).any(-1)
    BB = B
    BB[C] = -1
    node_idx_append = BB[BB > -1]

    tetrahedrons = np.column_stack((A, node_idx_append))

    t1 = time.time()
    # logger.info( '** Time elapsed: {:.2f} s **'.format ( t1-t0 ) )

    tetrahedrons = np.ndarray.astype(tetrahedrons, dtype=int)
    return tetrahedrons - 1


def _deprecated_face_group_to_tetrahedrons2(
    face_group: h5py.Group, face_zone_names: List[str]
) -> np.array:
    """ "Converts Fluent's face connectivity matrix to tetrahedron elements.

    Notes
    -----
    Format Fluent face:
    [n0 n1 n2] [c0 c1]
    n0 n1 n2 : node indices that make up face
    c0 c1    : cell indices to which face is connected. Max 2.

    Note: this finds n3 to construct a tetrahedron
    based on the connectivity of the faces. I.e.: the tetrahedron
    is defined by four faces as:
    f1: n1 n2 n3
    f2: n1 n2 n4
    f3: n2 n3 n4
    f3: n1 n3 n4

    Hence the tetrahedron will be composed of n1 n2 n3 n4.
    This function exploits this characteristic

    Parameters
    ----------
    face_group : h5py.Group
        handle to face group
    face_zone_names : List
        list of face zone names

    Returns
    -------
    np.array
        array with tetrahedron definitions

    Raises
    ------
    Exception
        Non-triangular faces detected
    """

    # get face definitions of all face zones
    faces_all = np.empty((0, 3))
    faces_zoneids = np.empty((0, 0), dtype=int)
    cell_zones = {}
    for ii, facezone_name in enumerate(face_zone_names):
        if "interior" not in facezone_name:
            continue
        print("Processing: {}".format(facezone_name))
        subdir1 = "nodes/" + str(ii + 1) + "/nodes"
        subdir2 = "nodes/" + str(ii + 1) + "/nnodes"
        faces = np.array(face_group[subdir1], dtype=int)
        nnodes = np.array(face_group[subdir2], dtype=int)

        if not np.all(nnodes == 3):
            raise Exception("Functionality for non-triangular faces not implemented...")

        # concatenate all faces into single array
        faces = np.reshape(faces, (int(len(faces) / 3), 3))
        num_faces = faces.shape[0]

        faces_zoneids = np.append(faces_zoneids, np.ones(faces.shape[0], dtype=int) * (ii + 1))

        # get c0c1
        subdir_c0 = "c0/" + str(ii + 1)
        subdir_c1 = "c1/" + str(ii + 1)

        c0c1 = np.array([face_group[subdir_c0], face_group[subdir_c1]]).T

        # np.savetxt(
        #     "c0c1_func1_{}.txt".format(
        # facezone_name.replace("interior-", "")),c0c1, delimiter=","
        # )

        # exploit unique to find indices of face pairs
        c0c1 = c0c1.T.ravel()
        cell_ids1, idx1, counts1 = np.unique(c0c1, return_index=True, return_counts=True)

        cell_ids2, idx2, counts2 = np.unique(np.flipud(c0c1), return_index=True, return_counts=True)

        faces_temp = np.vstack([faces, faces])

        num_elements = cell_ids1.shape[0]

        f1 = faces_temp[idx1, :]
        f2 = np.flipud(faces_temp)[idx2, :]

        node_ids_in_tet = np.hstack([f1, f2])

        # only select non-reocurring node-id
        unique = np.sort(node_ids_in_tet)
        duplicates = unique[:, 1:] == unique[:, :-1]
        unique[:, 1:][duplicates] = 0
        mask = unique != 0

        # np.savetxt(
        #     "faces_func1_{}.txt".format(facezone_name.replace("interior-", "")),
        #     faces,
        #     delimiter=",",
        # )
        # np.savetxt(
        #     "node_ids_in_tet_func1_{}.txt".format(facezone_name.replace("interior-", "")),
        #     node_ids_in_tet,
        #     delimiter=",",
        # )

        if not np.all(np.sum(mask, axis=1) == 4):
            raise ValueError("Tetrahedron does not consist of 4 node ids")

        tetrahedrons = np.reshape(node_ids_in_tet[mask], (num_elements, 4))
        cell_zones[facezone_name] = {"tetra": tetrahedrons, "tetra-ids": cell_ids1}

    t1 = time.time()
    # logger.info( '** Time elapsed: {:.2f} s **'.format ( t1-t0 ) )

    return cell_zones


def add_solid_name_to_stl(filename, solid_name, file_type: str = "ascii"):
    """Adds name of solid to stl file. Supports only single block"""
    if file_type == "ascii":
        start_str = "solid"
        end_str = "endsolid"
        f = open(filename, "r")
        list_of_lines = f.readlines()
        f.close()
        list_of_lines[0] = "{0} {1}\n".format(start_str, solid_name)
        list_of_lines[-1] = "{0} {1}\n".format(end_str, solid_name)

        f = open(filename, "w")
        f.writelines(list_of_lines)
        f.close()
    # replace part name in binary file
    elif file_type == "binary":
        fid = open(filename, "r+b")
        fid.seek(0)
        data = fid.read(40)
        fid.seek(0)
        string_replace = "{:<40}".format(solid_name).encode()
        fid.write(string_replace)
        fid.close()

    return


if __name__ == "__main__":
    print("protected")
