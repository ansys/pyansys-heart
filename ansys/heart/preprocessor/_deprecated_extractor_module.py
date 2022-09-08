from ansys.heart.preprocessor.mesh.geodisc import (
    get_closed_path,
    project_3d_points,
    rodrigues_rot,
    sort_aniclkwise,
)
import numpy as np
from scipy.spatial import KDTree
import scipy.spatial.distance
import vtk
from vtk.numpy_interface import dataset_adapter as dsa  # this is an improved numpy integration
from vtk.util import numpy_support as VN  # noqa


def find_superior_cap(p_set, mesh):
    """
    find node ID of surface mesh
    :return:
    """
    P_project, normal, P_mean = project_3d_points(p_set)
    points = mesh["node"]

    a = np.dot(points - P_mean, normal)
    set = np.where([a <= 0])[1]
    tmp = np.where([a > 0])[1]
    if len(set) > len(tmp):
        set = tmp
    # np.savetxt('bc_.csv', set, delimiter=',')
    return set


def get_nodes_cap_edge(
    points_close_to_surface: np.array,
    surface: vtk.vtkPolyData,
    surface_points_coord=None,
):
    """Gets the nodes which close each of the cavities. Returns the node
    indices that close the cap. Local node numbering is used

    """
    surface_np = dsa.WrapDataObject(surface)

    nodes_surface = surface_np.Points

    # project the intersection points
    p_project, cap_normal, p_mean = project_3d_points(points_close_to_surface)

    # flip normal direction if necessary, it should point to the cavity center
    center = np.mean(nodes_surface, axis=0)
    if np.dot(cap_normal, center - p_mean) < 0:
        cap_normal = -cap_normal

    # node_set = get_close_nodes(points, p_project)

    s_tree = KDTree(p_project)
    dst, _ = s_tree.query(nodes_surface, 1)
    # todo: this is based on the size of mesh, I choose 1 mm
    # Idea is to find a dozens of nodes
    node_set = np.argwhere(dst < 1).ravel()

    # sort node set in anti clock direction
    p_xy = rodrigues_rot(nodes_surface[node_set] - p_mean, cap_normal, [0, 0, 1])
    _, b = sort_aniclkwise(p_xy[:, 0:2].tolist())
    nodeset_ordered = node_set[b]

    # Valve nodes must on surface, this can not be ensured in case of multi-cavity
    # NOTE Martijn: Why is this the case?
    if surface_points_coord is not None:
        dst = scipy.spatial.distance.cdist(nodes_surface[nodeset_ordered], surface_points_coord)
        dd = np.min(dst, axis=1)
        # we eliminate the nodes if they are not on the surface
        nodeset_ordered = nodeset_ordered[np.where(dd == 0)]

    # np.savetxt("y.csv", points[nodeset_ordered])
    # exit()

    # keep more or less 6 nodes to find geodesic path
    n = int(len(nodeset_ordered) / 6)
    nodeset_ordered = nodeset_ordered[::n]

    # np.savetxt("nodes_surface.csv", nodes_surface[nodeset_ordered,:], delimiter=",")

    # id on the surface mesh
    local_node_ids = get_closed_path(nodeset_ordered, surface)

    # reorder nodes such that it forms a continuous loop.
    # NOTE: not sure why this not necessary before

    # get global ids of nodes that close the cavity

    return cap_normal, local_node_ids


def extract_endo_epi_coords(s, split_ele, parts=2):
    """ """
    # delete cells
    data = vtk.vtkPolyData()
    data.DeepCopy(s.GetOutput())
    data.BuildLinks()

    for i in split_ele:
        data.DeleteCell(i)
    data.RemoveDeletedCells()

    # # debug
    # writer = vtk.vtkPolyDataWriter()
    # writer.SetInputData(data)
    # writer.SetFileName("tmp.vtk")
    # writer.SetFileTypeToBinary()
    # writer.Write()
    # # end debug

    dct = {}
    for i in range(parts):
        connectivity0 = vtk.vtkPolyDataConnectivityFilter()
        connectivity0.SetExtractionModeToSpecifiedRegions()
        connectivity0.AddSpecifiedRegion(i)
        connectivity0.SetInputData(data)
        connectivity0.Update()
        pdata0 = connectivity0.GetOutput()

        nb_cell = pdata0.GetNumberOfCells()
        connect = np.zeros((nb_cell, 3), dtype=int)
        for icell in range(nb_cell):
            ids = pdata0.GetCell(icell).GetPointIds()
            connect[icell, 0] = ids.GetId(0)
            connect[icell, 1] = ids.GetId(1)
            connect[icell, 2] = ids.GetId(2)

        # local ID of endo/epi
        id = np.unique(np.ravel(connect))

        # get their coordinate
        coords = VN.vtk_to_numpy(pdata0.GetPoints().GetData())
        dct[str(i)] = coords[id]

    return dct


def get_endo_epi_node_id(ventricle, valve_nodes, parts=2):
    """
    deprecation method
    """

    local_id = np.nonzero(valve_nodes[:, None] == ventricle.surface["ids_to_volume"])[1]
    #  select elements if more than one node is in the list
    split_ele = np.where(np.any(np.isin(ventricle.surface["connect"], local_id), axis=1))[0]

    # get endo/epi node coords
    dct = extract_endo_epi_coords(ventricle.surface_polydata, split_ele, parts=parts)

    # get their ID from volume mesh
    coords = VN.vtk_to_numpy(ventricle.volume.GetPoints().GetData())
    tree = KDTree(coords)
    _, l1 = tree.query(dct["0"], 1)
    _, l2 = tree.query(dct["1"], 1)

    if parts == 2:
        # endo part is smaller
        if len(l1) > len(l2):
            return l2, l1
        else:
            return l1, l2
    elif parts == 3:  # decide the type of 3 surface later
        _, l3 = tree.query(dct["2"], 1)
        return l1, l2, l3


def get_segment(nodelist, ventricle, edges=None):
    """ """

    connect = ventricle.surface["connect"]
    global_id = ventricle.surface["ids_to_volume"]

    # p,e = ventricle.get_tetra_mesh()
    # np.savetxt("lv.csv",p[nodelist])
    # exit()
    # get the local Id

    localid = np.nonzero(nodelist[:, None] == global_id)[1]
    # np.savetxt("lv.Csv",ventricle.surface['points'][localid])
    # exit()

    # select elements if more than one node is in the list
    selected_elem = np.where(np.any(np.isin(connect, localid), axis=1))[0]

    # special case at the edge, only for epi segment
    if edges is not None:
        # localId
        xx = np.nonzero(edges[:, None] == global_id)[1]
        # select element if all three nodes are in the edge
        selected_elem2 = np.where(np.all(np.isin(connect, xx), axis=1))[0]
        selected_elem = np.hstack((selected_elem, selected_elem2))

    # get the connectivity table
    connect = connect[selected_elem]

    # connectivity table with  global ID
    connect2 = global_id[connect]

    return connect2


def get_endo_epi_nodeset(ventricle, valve_nodes, parts=2):

    split_element = get_elements_attached_to_nodes(valve_nodes, ventricle)

    # polydata to be split
    data = vtk.vtkPolyData()
    data.DeepCopy(ventricle.surface_polydata.GetOutput())
    data.BuildLinks()

    for i in split_element:
        data.DeleteCell(i)
    data.RemoveDeletedCells()

    # assign a array to save ids of volume mesh
    vtk_array = vtk.vtkIntArray()
    vtk_array.SetName("meshId")
    n_points = len(ventricle.surface["ids_to_volume"])

    vtk_array.SetNumberOfValues(n_points)
    for ii in range(n_points):
        vtk_array.SetValue(ii, ventricle.surface["ids_to_volume"][ii])
    data.GetPointData().AddArray(vtk_array)

    # # debug
    # writer = vtk.vtkPolyDataWriter()
    # writer.SetInputData(data)
    # writer.SetFileName("tmp2.vtk")
    # writer.SetFileTypeToBinary()
    # writer.Write()
    # # end debug

    node_set = [None] * parts
    count = 0
    for i in range(parts):
        connectivity = vtk.vtkPolyDataConnectivityFilter()
        connectivity.SetInputData(data)
        connectivity.SetExtractionModeToSpecifiedRegions()
        connectivity.AddSpecifiedRegion(i)
        connectivity.Update()

        pdata = connectivity.GetOutput()
        nb_cell = pdata.GetNumberOfCells()
        connect = np.zeros((nb_cell, 3), dtype=int)
        for icell in range(nb_cell):
            ids = pdata.GetCell(icell).GetPointIds()
            connect[icell, 0] = ids.GetId(0)
            connect[icell, 1] = ids.GetId(1)
            connect[icell, 2] = ids.GetId(2)

        # get node ids
        local_id = np.unique(np.ravel(connect))
        count += len(local_id)

        # find their global(volume mesh) id
        global_id = VN.vtk_to_numpy(pdata.GetPointData().GetArray("meshId"))
        node_set[i] = global_id[local_id]

    # verify results
    err = len(ventricle.surface["points"]) - len(valve_nodes) - count
    print(f"{err} nodes on surface are either endocardium or epicardium")

    # a = valve_nodes
    # for i in range(5):
    #     a = np.hstack((a, node_set[i]))
    # is_ok = np.array_equal(a.sort(), ventricle.surface["ids_to_volume"].sort())
    # assert is_ok==True

    # if err !=0:
    #     gb_ids = valve_nodes
    #     for i in range(parts):
    #         gb_ids = np.hstack((gb_ids, node_set[i]))
    #     surface_ids = np.nonzero(gb_ids[:, None] == ventricle.surface["ids_to_volume"])[1]
    #     nnode = len(ventricle.surface["points"])
    #     diff_ids = np.setdiff1d(np.linspace(0, nnode - 1, nnode, dtype=int), surface_ids)
    #     print(diff_ids) #surface id
    #     np.savetxt("p.csv",ventricle.surface['points'][diff_ids])
    #     exit()
    #     print(ventricle.surface["ids_to_volume"][diff_ids]) # gloabl id
    #     print(len(diff_ids))

    return node_set


def get_elements_attached_to_nodes(node_ids, all_nodes):
    # valve node id on surface mesh
    local_id = np.nonzero(node_ids[:, None] == all_nodes)[1]
    # select elements if more than one node is in the list
    split_element = np.where(np.any(np.isin(all_nodes.surface["connect"], local_id), axis=1))[0]

    return split_element
