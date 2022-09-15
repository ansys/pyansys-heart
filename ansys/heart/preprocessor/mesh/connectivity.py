"""Module containing methods for mesh connectivity"""

import copy

from ansys.heart.custom_logging import LOGGER
import numpy as np

# import tqdm as tqdm


def get_faces_tetra(tetra):
    """Gets faces that make up the tetra"""
    num_tetra = tetra.shape[0]
    faces = np.zeros((num_tetra, 3, 4), dtype=int)
    masks = np.array(
        [
            [True, True, True, False],
            [True, True, False, True],
            [True, False, True, True],
            [False, True, True, True],
        ]
    )
    for ii, mask in enumerate(masks):
        faces[:, :, ii] = tetra[:, mask]

    return faces


def tetra_to_faces(tetra):
    """Creates list of unique faces from tetrahedrons and returns tetra_face_map"""
    faces = get_faces_tetra(tetra)

    # reshape to (4*NumTetra, 3)
    num_tetra = tetra.shape[0]
    faces_1 = np.reshape(faces.transpose(0, 2, 1), (4 * num_tetra, 3))

    # construct tetra_face_map. Constructs tetrahedron such that it points to face array
    tetra_face_map = np.reshape(np.arange(0, num_tetra * 4, 1), (num_tetra, 4))

    # sort faces to facilitate finding uniques
    faces_1 = np.sort(faces_1, axis=1)

    # find unique faces
    unique_faces, index, inverse, counts = np.unique(
        faces_1, axis=0, return_index=True, return_inverse=True, return_counts=True
    )

    # update tetra_face_map accordingly
    tetra_face_map = inverse[tetra_face_map]

    # find connectivity c0 an c1 (cell indices that are connected to the face)
    tet_ids = np.arange(0, num_tetra, 1)
    face_ids = np.arange(0, len(unique_faces), 1)

    # this gives faces_1 again:
    # unique_faces[tetra_face_map]

    return tetra_face_map, unique_faces


def face_tetra_connectivity(tetra: np.array):
    """Computes the tetra-face connectivity tables"""

    import time as time

    LOGGER.debug("Establishing tetra-face connectivity...")
    t0 = time.time()

    # num_tetra = tetra.shape[0]
    # faces = np.zeros( (num_tetra, 3, 4), dtype=int )
    # masks = np.array([
    #     [ True, True, True, False],
    #     [ True, True, False, True],
    #     [ True, False, True, True],
    #     [ False, True, True, True]
    #      ] )
    # for ii, mask in enumerate( masks ):
    #     faces[:,:,ii] = tetra[:, mask]
    faces = get_faces_tetra(tetra)
    # reshape faces such that shape is (NumTetra*4, 3) with following structure:
    # i n1 n2 n3
    # 0 tet1_face1
    # 1 tet1_face2
    # 2 tet1_face3
    # 3 tet1_face4
    # 4 tet2_face1
    # 5 tet2_face2
    # 6 tet2_face3
    # 7 tet2_face4
    # ...
    num_tetra = tetra.shape[0]
    faces_1 = np.reshape(faces.transpose(0, 2, 1), (4 * num_tetra, 3))

    # sort faces in order to find duplicates
    faces_sorted = np.sort(faces_1, axis=1)
    np.sort(faces_sorted, axis=0)

    # find duplicate rows
    faces_unique, index, inverse, counts = np.unique(
        faces_sorted, return_index=True, return_inverse=True, return_counts=True, axis=0
    )

    # find duplicate rows in reversed order
    faces_sorted_flip = np.flipud(faces_sorted)
    faces_unique_r, index_r, inverse_r, counts_r = np.unique(
        faces_sorted_flip, return_index=True, return_inverse=True, return_counts=True, axis=0
    )

    # number of cells connected to each face
    num_cells_connected = counts[inverse]

    tetra_ids = np.repeat(np.arange(0, num_tetra, 1), 4)
    tetra_ids_flip = np.flipud(tetra_ids)

    # get connected tetra id for each face (two for interior, one for boundary face)
    c0 = tetra_ids[index][inverse]
    c1 = np.flip(tetra_ids_flip[index_r][inverse_r])

    # removing any duplicate faces
    mapper = np.sort(index)
    faces_1 = faces_1[mapper, :]
    c0 = c0[mapper]
    c1 = c1[mapper]

    # print("num cells connected")
    # print("  " + str(num_cells_connected))
    # print("c0" + str(c0))
    # print("c1" + str(c1))

    # NOTE: cell ordering makes sense in the following, but way slower.
    # c0 = tetra_ids
    # c1 = tetra_ids[inverse]
    # c1 = copy.deepcopy(c0)
    # c1[ num_cells_connected == 1 ] = c0[num_cells_connected == 1]
    # # populate c1 with a loop
    # for ii in np.where(num_cells_connected == 2)[0]:
    #     mask = np.all( faces_sorted[ii,:] == faces_sorted, axis = 1)
    #     mask[ii] = False
    #     c1[mask] = c0[ii]

    # # for all interior faces find matching cell
    # # TODO: How to get c1 efficiently without using a loop
    # print("num cells connected")
    # print("  " + str(num_cells_connected))
    # print("c0" + str(c0))
    # print("c1" + str(c1))
    t1 = time.time()
    LOGGER.debug("Time elapsed: {:.1f} s".format(t1 - t0))

    return faces_1, c0, c1


def get_face_type(faces: np.array, face_cell_connectivity: np.array) -> np.array:
    """Establishes face type. Either boundary faces or interior faces

    Parameters
    ----------
    faces : np.array
        Array with face definitions
    face_cell_connectivity : np.array
        Array describing to which cells each of the faces is connected to, e.g. np.array([c0, c1])

    Returns
    -------
    np.array
        Type of face. Either interior (face_type = 1) or boundary (face_type = 2)
    """

    interior_face_ids = face_cell_connectivity[:, 0] != face_cell_connectivity[:, 1]
    boundary_face_ids = face_cell_connectivity[:, 0] == face_cell_connectivity[:, 1]
    face_types = np.zeros((faces.shape[0]), dtype=int)
    face_types[interior_face_ids] = 1
    face_types[boundary_face_ids] = 2
    num_assigned = np.sum(boundary_face_ids) + np.sum(interior_face_ids)
    assert num_assigned == faces.shape[0], "Not all faces assigned as either interior or boundary"
    return face_types


def get_edges_from_triangles(triangles: np.array) -> np.array:
    """Generates an array of edges from a array of triangles"""
    num_triangles = triangles.shape[0]
    num_edges = num_triangles * 3
    edges = np.repeat(triangles, 3, axis=0)
    mask = np.tile(
        np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]], dtype=bool),
        (num_triangles, 1),
    )
    edges = np.reshape(edges[mask], (num_edges, 2))

    return edges


def get_free_edges(triangles: np.array, return_free_triangles: bool = False) -> np.array:
    """Gets the boundary edges that are only referenced once"""

    edges = get_edges_from_triangles(triangles)

    edges_sort = np.sort(edges, axis=1)

    unique_edges, idx, counts = np.unique(edges_sort, axis=0, return_counts=True, return_index=True)
    free_edges = edges[idx, :][counts == 1, :]

    if not return_free_triangles:
        return free_edges

    elif return_free_triangles:
        # get free triangles
        free_triangles = triangles[
            np.argwhere(np.sum(np.isin(triangles, free_edges), axis=1) == 2).flatten(), :
        ]

        return free_edges, free_triangles


def edge_connectivity(
    edges: np.array, return_type: bool = False, sort_closed: bool = False
) -> np.ndarray:
    """Group edges by connectivity

    Parameters
    ----------
    edges : np.array
        NumEdges x 2 Numpy array with edge definitions
    return_type : bool, optional
        Flag indicating whether to return the type of the edge group, by default False:
            "open": edge group is open-ended
            "closed": edge group forms closed edge loop
    sort_closed : bool, optional
        Flag indicating whether to sort any closed edge loops, by default False

    Note
    ----
    Uses an implementation of Dept-first search:
    https://en.wikipedia.org/wiki/Depth-first_search
    https://www.educative.io/answers/how-to-implement-depth-first-search-in-python

    Performance is not tested so may not be suitable for large arrays of edges
    """
    # Nested in Dept-first search algorithm
    def _dfs(visited, graph, node):
        if node not in visited:
            # print(node)
            visited.add(node)
            for neighbor in graph[node]:
                _dfs(visited, graph, neighbor)

    # create adjacency list (typically referred to as "graph")
    graph = {}
    node_ids = np.unique(edges)
    for node in node_ids:
        mask = node != edges
        mask[np.all(mask, axis=1), :] = False
        connected_nodes = edges.flatten()[mask.flatten()]
        graph[node] = connected_nodes

    # check connectivity of each node using DFS.
    # Group connected edges
    node_ids_visited = np.zeros(node_ids.shape[0], dtype=bool)
    edge_groups = []
    while not np.all(node_ids_visited):
        # keep track of visited nodes for this group of edges
        visited = set()

        # node id to start from (finds first un-visited node)
        start_node_id = node_ids[np.where(np.invert(node_ids_visited))[0][0]]

        # call dept first algorithm to find connectivity
        _dfs(visited, graph, start_node_id)

        # retrieve edge definitions
        edge_group = edges[np.all(np.isin(edges, list(visited)), axis=1)]
        edge_groups.append(edge_group)

        node_ids_visited[np.isin(node_ids, list(visited))] = True

    # check whether edges form a closed loop
    group_types = []
    if return_type:
        for edge_group in edge_groups:
            counts = np.unique(edge_group, return_counts=True)[1]
            if np.all(counts == 2):
                # print("Closed edge loop")
                group_types.append("closed")
            elif np.any(counts != 2):
                group_types.append("open")
                # print("Open edge tree")

    # sort any closed edge loops
    if sort_closed:
        for ii, edge_group in enumerate(edge_groups):
            if group_types[ii] == "closed":
                edges = edge_group.tolist()
                remaining_edges = edges
                next_edge = edges[0]
                sorted_edges = [next_edge]
                remaining_edges.pop(0)
                while len(remaining_edges) > 0:
                    # find connected edge of last edge
                    node = sorted_edges[-1][1]
                    mask = np.array(edges) == node
                    if np.sum(mask[:, 1]) == 1:
                        flip = True
                    elif np.sum(mask[:, 0]) == 1:
                        flip = False
                    else:
                        raise ValueError("Expecting just one match")
                    idx = np.where(np.any(mask, axis=1))[0][0]
                    if flip:
                        sorted_edges.append(np.flip(remaining_edges[idx]).tolist())
                    else:
                        sorted_edges.append(remaining_edges[idx])
                    remaining_edges.pop(idx)
                edge_groups[ii] = np.array(sorted_edges)

    if return_type:
        return edge_groups, group_types
    else:
        return edge_groups


def remove_triangle_layers_from_trimesh(triangles: np.array, iters: int = 1) -> np.array:
    """Identifies triangles connected to the boundary, and removes these from the triangle list

    Parameters
    ----------
    triangles : np.array
        Array of triangles
    iters : int, optional
        Number of iterations, by default 1

    Returns
    -------
    np.array
        Reduced set of triangles
    """
    reduced_triangles = copy.deepcopy(triangles)
    for ii in range(0, iters, 1):
        num_triangles = reduced_triangles.shape[0]
        edges = get_edges_from_triangles(reduced_triangles)
        free_edges = get_free_edges(reduced_triangles)

        # find elements connected to the free edges
        edges = np.reshape(edges, (3, 2, num_triangles))
        free_nodes = np.unique(free_edges)

        idx_triangles_boundary = np.any(np.isin(reduced_triangles, free_nodes), axis=1)
        # idx_triangles_boundary = np.any(
        #     np.all( np.isin(edges, free_edges),
        #         axis = 1),
        #     axis = 0 )

        LOGGER.debug("Removing {0} connected triangles".format(np.sum(idx_triangles_boundary)))

        # remove boundary triangles
        reduced_triangles = reduced_triangles[~idx_triangles_boundary, :]

    return reduced_triangles


if __name__ == "__main__":

    # ************** Simple 2 tetrahedron example *******
    nodes = np.array([[1, 0, 0], [-1, 0, 0], [0, 0, 1], [0, -1, 0], [0, 1, 0]], dtype=float)
    tetra = np.array([[0, 1, 2, 3], [0, 1, 2, 4]])
    part_ids = np.array([1, 2])

    # get face-cell connectivity table
    faces, c0, c1 = face_tetra_connectivity(tetra)
    face_types = get_face_type(faces, np.array([c0, c1]).T)

    interior_faces = faces[face_types == 1, :]
    boundary_faces = faces[face_types == 2, :]
    interface_faces = faces[part_ids[c0] != part_ids[c1], :]

    LOGGER.debug("Protected")
