import h5py
import numpy as np
import time
import meshio
from typing import List, Optional

from ansys.heart.custom_logging import logger


"""Module containing functions to read/write fluent meshes in HDF5 format
"""

def fluenthdf5_to_vtk(hdf5_filename: str, vtk_filename: str):
    """Converts Fluent hdf5 to vtk format

    Parameters
    ----------
    hdf5filename : str
        path to hdf5 file
    vtkfilename : str
        path to vtk file
    """

    fid = h5py.File(hdf5_filename, "r")

    mesh_group = get_mesh_group(fid)

    num_cell_zones = len(mesh_group)

    # get faces / tetrahedrons from all mesh groups
    for key in mesh_group.keys():
        cell_zone_id = int(key)

        face_group = get_face_group(mesh_group[key])
        face_zone_info = get_face_zone_info(face_group)
        tetrahedrons = face_group_to_tetrahedrons(face_group, face_zone_info[0])
        face_zones = face_group_to_face_zones(face_group, face_zone_info[0], get_interior=False)

        points = get_nodes_from_mesh_group(mesh_group[key])[0]

        cells = []
        cell_data = []
        for _, value in face_zones.items():
            # logger.info(value)
            cells.append(("triangle", value["faces"] - 1))
            zoneids = np.ones(cells[-1][-1].shape[0], dtype=int) * value["zone-id"]
            cell_data.append(zoneids)

        cells.append(("tetra", tetrahedrons - 1))

        cell_data.append(np.ones(tetrahedrons.shape[0] * cell_zone_id))
        cell_data = {"face-zone-id": cell_data}

        mesh = meshio.Mesh(
            points=points,
            cells=cells,
            # ("triangle", tris-1 ),
            # ("tetra", tetrahedrons - 1)
            # ],
            cell_data=cell_data,
        )
        mesh.write("fluent_volume_mesh.vtk")

    fid.close()

    return


#     return
def get_mesh_group(fid: h5py.File) -> h5py.Group:
    keys = list(fid.keys())
    if keys[0] != "meshes":
        raise KeyError("No 'meshes' field found")
    return fid["meshes"]


def get_face_group(group: h5py.Group) -> h5py.Group:
    """Gets face zone group """

    def find_faces(name):
        if "faces" in name:
            return name

    face_group = group[group.visit(find_faces)]
    return face_group


def get_face_zone_info(face_group: h5py.Group) -> List:
    """Gets list of face zone names """
    # get face zone ids
    face_zone_ids = np.array(face_group["zoneTopology/id"])
    # get face zones
    face_zone_names = np.array(face_group["zoneTopology/name"])
    face_zone_names = np.chararray.tobytes(face_zone_names).decode()
    face_zone_names = face_zone_names.split(";")

    face_zone_types = np.array(face_group["zoneTopology/zoneType"])
    return face_zone_names, face_zone_ids, face_zone_types


def get_nodes_from_mesh_group(mesh_group: h5py.Group):
    """Gets nodes from the mesh group"""

    nodes = np.empty((0, 3), dtype=float)
    for key in mesh_group["nodes/coords"].keys():
        node_group = np.array(mesh_group["nodes/coords/" + key])
        nodes = np.append(nodes, node_group, axis=0)

    # boundary_nodes = np.array( mesh_group['nodes/coords/1'] )
    # interior_nodes = np.array( mesh_group['nodes/coords/2'] )
    # nodes          = np.append( boundary_nodes, interior_nodes, axis=0 )
    node_ids = np.arange(0, np.shape(nodes)[0], 1, dtype=int)
    return nodes, node_ids


def face_group_to_face_zones(
    face_group: h5py.Group,
    face_zone_names: List[str],
    get_all_zones: bool = False,
    get_interior: Optional[bool] = False,
) -> dict:

    # TODO: write these in a more sensible structure, e.g. a dictionary?

    # get face definitions of all face zones
    faces_all = np.empty((0, 3))
    faces_zoneid = np.empty((0, 0), dtype=int)
    faces_label = []

    # checks whether provided list of face zones is present in the face group
    face_zone_names_1 = get_face_zone_info(face_group)[0]
    face_zone_ids = get_face_zone_info(face_group)[1]

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
            "faces": faces,
            "zone-id": face_zone_ids[ii],
        }

    return faces_out


def face_group_to_tetrahedrons(face_group: h5py.Group, face_zone_names: List[str]) -> np.array:
    """"Converts Fluent's face connectivity matrix to tetrahedron elements.
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
        subdir1 = "c0/" + str(ii + 1)
        subdir2 = "c1/" + str(ii + 1)
        c0c1 = np.column_stack((np.array(face_group[subdir1]), np.array(face_group[subdir2])))
        c0c1_all = np.append(c0c1_all, c0c1, axis=0)

    uniq_elem_id, uniq_count = np.unique(c0c1_all[c0c1_all > 0], return_counts=True)
    num_elem = uniq_elem_id.size

    tetrahedrons = np.empty((num_elem, 4), dtype=int)

    logger.info("Generating {0} tetrahedrons...".format(num_elem))

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
    logger.info("** Time elapsed: {:.2f} s **".format(t1 - t0))

    tetrahedrons = np.ndarray.astype(tetrahedrons, dtype=int)
    return tetrahedrons


if __name__ == "__main__":
    print("Protected")
