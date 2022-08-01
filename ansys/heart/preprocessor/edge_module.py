"""Module for extracting edge connectivity"""
import numpy as np


def edge_connectivity(edges: np.array, return_type: bool = False, sort_closed: bool = False):
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
            for neighbour in graph[node]:
                _dfs(visited, graph, neighbour)

    # create adjacency list (typically refered to as "graph")
    graph = {}
    node_ids = np.unique(edges)
    for node in node_ids:
        mask = node != edges
        mask[np.all(mask, axis=1), :] = False
        connected_nodes = edges.flatten()[mask.flatten()]
        graph[node] = connected_nodes

    # check connecivity of each node using DFS.
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

    # check whether edges form a closed loopg
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
            if group_types[ii] != "closed":
                continue
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


if __name__ == "__main__":
    # create some example data
    edges = np.array(
        [[1, 4], [1, 5], [5, 6], [6, 7], [12, 10], [10, 8], [12, 11], [8, 3], [3, 11]], dtype=int
    )
    edge_loops, edge_loop_types = edge_connectivity(edges, return_type=True, sort_closed=True)
