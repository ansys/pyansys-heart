from aifc import Error
import os
import copy
import warnings
import numpy as np
import vtk
from typing import List, Union

from vtk.numpy_interface import dataset_adapter as dsa  # this is an improved numpy integration

# NOTE: do more specific imports!
# from ansys.heart.preprocessor.mesh_module import *
from ansys.heart.preprocessor.vtk_module import (
    create_vtk_surface_triangles,
    vtk_surface_filter,
    threshold_vtk_data,
    get_tri_info_from_polydata,
    add_vtk_array,
    vtk_surface_to_stl,
    write_vtkdata_to_vtkfile,
    compute_volume_stl,
)

# from ansys.heart.preprocessor.vtk_module import compute_volume_stl

# from ansys.heart.preprocessor.fluenthdf5_module import fluenthdf5_to_vtk
from ansys.heart.preprocessor.geodisc_module import (
    order_nodes_edgeloop,
    sort_edgeloop_anti_clockwise,
)

from ansys.heart.preprocessor.extractor_module import get_nodes_cap_edge

# these import the "old" files
# from preprocessing.extractor_module import *
from ansys.heart.preprocessor.model_information import ModelInformation
from ansys.heart.custom_logging import LOGGER


class ClosingCap:
    """Definition of the closing cap"""

    def __init__(self, name, capid, cap_type=None):
        """Initializes closing cap"""

        self.name = name
        """Name of the cap"""
        self.id = capid
        """"ID of the cap/valve"""
        self.cap_type = cap_type
        """Type of cap/valve"""
        self.node_ids_cap_edge = None
        """Node ids that make up the cap edge"""
        self.nodes_cap_edge = None
        """Node coordinates of the cap edge"""
        self.edge_loop = None
        """Array of nodes that form a closed loop around the cap edge"""
        self.closing_triangles = np.empty((0, 3))
        """Triangles that close the cap/valve"""

        self.nodeset_id: int = 0  # Only one nodeset id per cap?
        """LS-DYNA node set id"""
        self.segset_id: int = 0  # only one segment set per cap?
        """LS-DYNA segment set id"""
        self.part_id: int = 0
        """LS-DYNA Part id of the cap/valve"""
        self.center: np.empty((0, 3), dtype=float)
        """Center of the valve"""
        self.center_raw: np.empty((0, 3), dtype=float)
        """Original center of the valve"""

    @property
    def global_node_ids_source_mesh(self) -> np.array:
        """Global node indices of where the mesh of the cap
        intersects with the mesh of the myocardium. These nodes are
        later used to find similar nodes on the remeshed surface/volume
        """
        return self._global_node_ids

    @global_node_ids_source_mesh.setter
    def global_node_ids_source_mesh(self, value: np.array):
        self._global_node_ids = value

    @property
    def nodes_source_mesh(self) -> np.array:
        """Coordinates of the nodes where the mesh of the cap
        intersects with the mesh of the myocardium. These nodes are
        later used to find similar nodes on the remeshed surface/volume
        """
        return self._nodes

    @nodes_source_mesh.setter
    def nodes_source_mesh(self, value: np.array):
        self._nodes = np.array(value)  # explitely cast to numpy array

    # @property
    # def closing_triangles ( self ) -> np.array:
    #     """These are the triangles that close the cap/valve.
    #     This refers to the global node-ids of the volume mesh
    #     """
    #     return self._closing_triangles

    # @closing_triangles.setter
    # def closing_triangles(self, value: np.array):
    #     self._closing_triangles = value

    @property
    def faces(self) -> np.array:
        return self._faces

    @faces.setter
    def faces(self, value: np.array):
        self._faces = value

    @property
    def centroid(self) -> np.array:
        return self._centroid

    @centroid.setter
    def centroid(self, value: np.array):
        self._centroid = value

    @property
    def normal(self) -> np.array:
        return self._normal

    @normal.setter
    def normal(self, value: np.array):
        self._normal = value

    def dump_nodes_to_file(self, filename):
        """Dumps nodes to a file"""
        np.savetxt(filename, self.nodes_source_mesh, delimiter=",")
        return

    def compute_centroid_source(self):
        """Computes the centroid of the cap"""
        self.centroid = np.mean(self.nodes_source_mesh)
        return self.centroid

    def compute_centroid_edgeloop(self):
        self.centroid = np.mean(self.nodes_cap_edge, axis=0)
        return self.centroid

    def get_nodes_that_close_cap(
        self,
        surface_myocardium: vtk.vtkPolyData,
    ):
        """[OBSOLETE] Uses the surface and original/raw nodes of the initial
        surface mesh to extract a "simplified" closing cap. Requires the surface
        of the connected myocardium as input"""
        # get_valve_node( self.nodes, surface_myocardium)
        return


class Valve(ClosingCap):
    def __init__(
        self,
        name,
        capid,
        global_node_ids=None,
        nodes=None,
        faces=None,
        cap_type=None,
    ):
        super().__init__(name, capid, global_node_ids, nodes, faces, cap_type)


class Cavity:
    """Class that defines the cavity and relevant cavity functions"""

    def __init__(
        self,
        name: str,
        vtk_labels: list,
        vtk_ids: list,
        model_info: ModelInformation,
    ):

        self.name = name
        self.labels = vtk_labels
        self.vtk_ids = vtk_ids
        self.info = model_info

        self.centroid = np.empty(3)

        self.apex_id: int = {"endocardium": 0, "epicardium": 0}
        """Node id of apical points"""

        self.id: int = 0  # cavity id
        self.nodeset_ids = {}
        self.segset_ids = {}

        self._surfaces = (
            dict()
        )  # stores vtk of all surfaces of this cavity: NOTE: original topology

        self._myocardium_surface = (
            vtk.vtkUnstructuredGrid()
        )  # stores vtk representation of the myocardium: NOTE: new topology

        # self.centroid = np.array( (0,0,0) )

        vtklabels_to_vtkids = {}
        for ii, label in enumerate(vtk_labels):
            vtklabels_to_vtkids[label] = vtk_ids[ii]

        # initialize closing caps
        self.closing_caps: List[ClosingCap] = []

        for label in vtk_labels[1:]:
            self.closing_caps.append(ClosingCap(name=label, capid=vtklabels_to_vtkids[label]))

        # store relevant node and segment sets of each cavity
        # list of dictionaries: is this convenient?
        # set contains the list of nodes or segments
        # NOTE: Could replace with namedtuple
        self.node_sets = [
            {"name": "endocardium", "set": np.empty(0), "id": 0},
            {"name": "epicardium", "set": np.empty(0), "id": 0},
        ]

        # same structure as node set: but use append to fill this
        self.segment_sets = []
        self.element_sets = []

    @property
    def labels(self) -> dict:
        return self._labels

    @labels.setter
    def labels(self, value: dict):
        self._labels = value

    def _compute_cap_intersection_with_myocardium(self, nodes, elements, tags):
        """Finds the points where the caps of this cavty intersect with the myocardium
        Needs (global) nodes element definitions and corresponding tags"""

        tags = tags
        elems = elements
        nodes = nodes

        # reference elements - that is, either the Left ventricle myocardium
        # right ventricle myocardium, left atrium myocardium or right atrium myocardium
        ref_elems = elems[tags == self.vtk_ids[0]]

        # loop over each cap of the cavity
        for cap in self.closing_caps:
            LOGGER.debug("Trying to find intersecting nodes for cap {}...".format(cap.name))
            # logger.debug("Looking for intersections on cap {0}...".format(cap.name) )
            # loop over connected parts
            # for ii, label in enumerate( self.labels[1:] ):
            node_indices_intersect = np.intersect1d(
                np.unique(ref_elems.ravel()),
                np.unique(elems[tags == cap.id].ravel()),
            )

            if len(node_indices_intersect) > 0:
                cap.global_node_ids_source_mesh = node_indices_intersect
                cap.nodes_source_mesh = nodes[node_indices_intersect, :].tolist()
            elif len(node_indices_intersect) == 0:
                warnings.warn("Warning: no connecting nodes found for: {0}".format(cap.name))

        return

    def _set_connected_surfaces(self, surface: Union[vtk.vtkPolyData, vtk.vtkUnstructuredGrid]):
        """Adds the surfaces that are involved in this cavity
        as vtk objects to the object"""
        self._surfaces = {}
        for ii, label in enumerate(self.labels):
            id_to_extract = self.vtk_ids[ii]
            surface_threshold, _ = threshold_vtk_data(
                surface, id_to_extract, id_to_extract, data_name="tags"
            )
            surface_threshold = vtk_surface_filter(surface_threshold)
            self._surfaces[label] = surface_threshold

        from ansys.heart.preprocessor.vtk_module import write_vtkdata_to_vtkfile
        from ansys.heart.preprocessor.mesh_module import add_solid_name_to_stl

        for key, polydata in self._surfaces.items():
            if "myocardium" in key:
                write_vtkdata_to_vtkfile(polydata, key + ".vtk")
                vtk_surface_to_stl(polydata, key + ".stl")
                add_solid_name_to_stl(key + ".stl", key, "binary")

        return

    def _find_closing_edge_loop(self):
        """Finds the edgeloops that closes each cavity"""
        # NOTE: Call this from the cap object or cavity object?

        labels = self._surfaces.keys()
        label = [
            label for label in labels if "myocardium" in label
        ]  # get the label that matches myocardium
        label = label[0]

        myocardium_surface = dsa.WrapDataObject(self._myocardium_surface)
        nodes_myocardium = myocardium_surface.Points

        for cap in self.closing_caps:
            LOGGER.debug("Finding nodes for cap: " + cap.name)

            if cap.node_ids_cap_edge is not None:
                LOGGER.debug("Skipping Cap: Node ids already computed and set")
                continue

            cap_normal, node_ids_cap = get_nodes_cap_edge(
                cap.nodes_source_mesh, self._myocardium_surface
            )

            # reorder the nodes on the edge loop again
            node_ids_cap = order_nodes_edgeloop(node_ids_cap, nodes_myocardium)

            # check direction of edge loop - ensure that seen from cavity centroid
            # it is in anti-clockwise direction
            reverse_order = sort_edgeloop_anti_clockwise(
                nodes_myocardium[node_ids_cap, :], self.centroid
            )
            if reverse_order:
                LOGGER.debug("Reversing order of cap nodes: %s" % cap.name)
                node_ids_cap = np.flip(node_ids_cap)

            # get global node ids of volume from "local" surface nodeids
            global_node_ids = myocardium_surface.PointData["GlobalPointIds"][node_ids_cap]

            # form edgeloop from sorted nodes
            nodes1 = global_node_ids
            nodes2 = np.append(global_node_ids[1:], global_node_ids[0])

            # add the edge loop to the cap
            cap.edge_loop = np.vstack((nodes1, nodes2)).transpose()

            # add nodes of edge to object
            cap.nodes_cap_edge = nodes_myocardium[node_ids_cap, :]

            # compute and add ca[ centroid
            cap.compute_centroid_edgeloop()

            cap.normal = cap_normal

            # NOTE: these are the global node ids of the volume mesh in HeartMesh
            cap.node_ids_cap_edge = global_node_ids

        return

    def _compute_centroid_cavity(self):
        """Computes centroid of cavity"""
        # get the myocardium surface
        for surface_name in self._surfaces.keys():
            if "myocardium" in surface_name:
                myocardium_surface = self._surfaces[surface_name]
                break
        self.centroid = np.array(myocardium_surface.GetCenter())
        return

    def _triangulate_caps(self, vtk_volume=vtk.vtkUnstructuredGrid, add_centroid=False):
        """Triangulates the closing caps by using a mid-point and the
        nodes of the edge loop.

        Note
        ------
        Can we write this in a single VTK object/file?

        """
        vtk_volume_obj = dsa.WrapDataObject(vtk_volume)

        all_tris = np.empty((0, 3), dtype=int)
        if add_centroid:
            LOGGER.error("Adding a centroid node is not supported yet")
            for cap in self.closing_caps:
                node_ids_edge = cap.node_ids_cap_edge

                # NOTE: How to add a point to the VTK Object?
                # this is not working yet
                vtk_volume_obj.SetPoints(np.vstack((vtk_volume_obj.Points, cap.centroid)))

                node_idx_centroid = vtk_volume_obj.Points.shape[0] + 1

                num_tris = len(node_ids_edge)
                tris_cap = np.empty((0, 3))
                for ii in range(0, num_tris - 1, 1):
                    tri = [
                        node_ids_edge[ii],
                        node_ids_edge[ii + 1],
                        node_idx_centroid,
                    ]
                    tris_cap = np.vstack((tris_cap, tri))

                # special treatment last triangle that closes the cap
                tri = [node_ids_edge[-1], node_ids_edge[0], node_idx_centroid]
                tris_cap = np.vstack((tris_cap, tri))

                cap.closing_triangles = tris_cap

                all_tris = np.vstack((all_tris, tris_cap))

            all_tris = np.array(all_tris, dtype=int)

        else:
            for cap in self.closing_caps:

                node_ids_edge = cap.node_ids_cap_edge
                num_triangles = len(node_ids_edge) - 2
                tris_cap = np.empty((0, 3))

                for ii in range(0, num_triangles):
                    tri = [
                        node_ids_edge[0],
                        node_ids_edge[ii + 1],
                        node_ids_edge[ii + 2],
                    ]
                    tris_cap = np.vstack((tris_cap, tri))

                cap.closing_triangles = tris_cap

                all_tris = np.vstack((all_tris, tris_cap))

            all_tris = np.array(all_tris, dtype=int)

        # vtk_surface = create_vtk_surface_triangles(vtk_volume_obj.Points, all_tris)
        # write_vtkdata_to_vtkfile(vtk_surface, "caps_of_" + self.name + "vtk")

        return

    def _get_endocardium_epicardium_points(self, surface: vtk.vtkPolyData):
        """Gets node sets on the endocardium and epicardium

        Parameters
        ----------
        surface : vtk.vtkPolyData
            Surface mesh extracted from the volume mesh
        Note
        ----
        This surface mesh should have the global point ids which are used for
        mapping to the volume mesh as a point data field
        """
        # extract tag
        LOGGER.debug("Extracting endo and epicardium for: " + self.labels[0])

        # identify the elements to remove
        nodes, tris, tris_data, point_data = get_tri_info_from_polydata(surface)

        idmap_surface_to_volume = point_data["GlobalPointIds"]

        # duplicate "GlobalPointIds" (this somehow gets overwritten in the geometery filter)
        add_vtk_array(
            surface,
            point_data["GlobalPointIds"],
            name="GlobalPointIds2",
            data_type="point",
            array_type=int,
        )

        # create deep copy of surface
        surface_copy = vtk.vtkPolyData()
        surface_copy.DeepCopy(surface)
        surface_copy.BuildLinks()

        # map triangles to global ids (volume mesh)
        tris_global = idmap_surface_to_volume[tris]

        # check which elements use any of the edgeloop nodes. Mark these for removal
        elements_to_remove = np.empty(0, dtype=int)
        for cap in self.closing_caps:

            # check which elements contain one or more of the edge-loop nodes
            mask = np.isin(tris_global, cap.node_ids_cap_edge)

            elements_to_remove = np.append(elements_to_remove, np.where(np.any(mask, axis=1))[0])

        for element_idx in elements_to_remove:
            surface_copy.DeleteCell(element_idx)
        surface_copy.RemoveDeletedCells()

        # write data
        # write_vtkdata_to_vtkfile(surface_copy, "surface_without_elements.vtk")

        # note vtk_ids[0] is the myocardium of the corresponding cavity
        # NOTE: Do not extract global ids here, they are already contained in the
        # point data
        (surface_tag, _) = threshold_vtk_data(
            surface_copy, self.vtk_ids[0], self.vtk_ids[0], "tags"
        )

        # use PolyData connectivity filter to
        # separate the surface into end/epicardium parts
        # need to create vtk object with ids before filter?

        boundary = vtk.vtkGeometryFilter()
        boundary.SetInputData(surface_tag)
        boundary.Update()

        bnd = dsa.WrapDataObject(boundary.GetOutput())
        if "GlobalPointIds2" not in list(bnd.PointData.keys()):
            raise ValueError("No data field found with global point ids")

        # extract all regions
        connectivity_filter = vtk.vtkPolyDataConnectivityFilter()
        # connectivity_filter = vtk.vtkConnectivityFilter()
        connectivity_filter.SetInputData(boundary.GetOutput())
        connectivity_filter.SetExtractionModeToAllRegions()
        connectivity_filter.SetColorRegions(1)
        connectivity_filter.MarkVisitedPointIdsOn()
        connectivity_filter.Update()
        connectivity_filter.GetNumberOfExtractedRegions()

        # get tri info from polydata?
        nodes, tris, tris_data, point_data = get_tri_info_from_polydata(
            connectivity_filter.GetOutput()
        )

        # get region id which is closest to centroid of cavity. This is
        # 1. region id closest to centroid is endocardium
        # 2. second largest nodeset is the epicardium
        # 3. in case of biventricular model or fourchamber model, the septum is the third largest
        # region

        region_ids = point_data["RegionId"]
        idmap_surface_to_volume = point_data["GlobalPointIds2"]
        # get endocardium:
        myocardium_centroid = np.mean(nodes, axis=0)
        d = np.linalg.norm(nodes - myocardium_centroid, axis=1)
        region_id_endocardium = region_ids[np.argmin(d)]

        # get epicardium:
        unique_region_ids, counts = np.unique(region_ids, return_counts=True)

        # sort by # counts:
        unique_region_ids = np.flip(unique_region_ids[np.argsort(counts)])
        counts = np.flip(counts[np.argsort(counts)])

        # visit other regios
        mask = unique_region_ids != region_id_endocardium
        remaining_region_ids = unique_region_ids[mask]
        remaining_counts = counts[mask]

        # Assumes that the largest remaining region is the epicardium.
        region_id_epicardium = remaining_region_ids[0]

        # Special treatment for left ventricle in case of Bi-Ventricle or FourChamber model
        if self.info.model_type in ["BiVentricle", "FourChamber"] and self.name == "Left ventricle":
            region_id_septum = remaining_region_ids[1]

        else:
            # this results in an empty node set
            region_id_septum = -1

        # map back to original node ids
        node_ids_endocardium = np.argwhere(region_id_endocardium == region_ids).flatten()
        num_nodes_endocardium = len(node_ids_endocardium)

        node_ids_epicardium = np.argwhere(region_id_epicardium == region_ids).flatten()
        num_nodes_epicardium = len(node_ids_epicardium)

        node_ids_septum = np.argwhere(region_id_septum == region_ids).flatten()
        num_nodes_septum = len(node_ids_septum)

        # map back to node ids of the volume mesh.
        global_node_ids_endocardium = idmap_surface_to_volume[node_ids_endocardium]

        global_node_ids_epicardium = idmap_surface_to_volume[node_ids_epicardium]

        global_node_ids_septum = idmap_surface_to_volume[node_ids_septum]

        # store in self
        for nodeset in self.node_sets:
            if nodeset["name"] == "endocardium":
                nodeset["set"] = global_node_ids_endocardium
            elif nodeset["name"] == "epicardium":
                nodeset["set"] = global_node_ids_epicardium

        # Septum nodeset only defined when BiVentricle or Four Chamber models are defined
        if self.info.model_type in ["BiVentricle", "FourChamber"] and self.name == "Left ventricle":
            self.node_sets.append(
                {"name": "epicardium-septum", "set": global_node_ids_septum, "id": 0}
            )

        LOGGER.debug(
            "\tNumber of nodes: [endocardium: {0}, epicardium: {1}]".format(
                num_nodes_endocardium, num_nodes_epicardium
            )
        )
        LOGGER.debug(
            "\tRegion ids: [endocardium: {0}, epicardium: {1}".format(
                region_id_endocardium, region_id_epicardium
            )
        )

        # number of orphan nodes:
        num_orphan_nodes = (
            nodes.shape[0] - num_nodes_endocardium + num_nodes_epicardium + num_nodes_septum
        )

        LOGGER.debug("\tNumber of nodes in all regions: {0}".format(counts))
        LOGGER.debug("\tNumber of orphan nodes: {0}".format(num_orphan_nodes))

        path_to_file = os.path.join(
            self.info.working_directory,
            "extracted_regions_{0}.vtk".format("_".join(self.name.lower().split())),
        )
        write_vtkdata_to_vtkfile(connectivity_filter.GetOutput(), path_to_file)

        return

    def _get_endocardium_epicardium_segments(self, vtk_surface: vtk.vtkPolyData):
        """From the available node sets generate segment sets"""
        surface_obj = dsa.WrapDataObject(vtk_surface)
        num_polygons = vtk_surface.GetNumberOfCells()

        try:
            tris = np.reshape(surface_obj.Polygons, (num_polygons, 4))
            tris = tris[:, 1:]
        except:
            raise Error("Failed to convert cells to triangles. Expecting only triangles")

        # map to global node ids
        tris_global = surface_obj.PointData["GlobalPointIds"][tris]

        # get the triangular elements belonging to each node set
        for node_set in self.node_sets:
            nodes_to_use = node_set["set"]

            # also include the nodes of the edge loop: otherwise those elements will be ignored
            nodes_cap_edges = np.empty(0, dtype=int)
            for cap in self.closing_caps:
                nodes_cap_edges = np.append(nodes_cap_edges, cap.node_ids_cap_edge)

            nodes_to_use = np.append(nodes_to_use, nodes_cap_edges)

            # if all three nodes are used in the triangle, then this is the corresponding surface
            element_indices = np.all(np.isin(tris_global, nodes_to_use), axis=1)

            # mesh = meshio.Mesh(points = surface_obj.Points,
            #                     cells = [("triangle", tris[element_indices, : ] ) ] )
            # mesh.write( "segment" + surface_name + ".vtk" )

            # need to append the segment set:
            segset_to_append = copy.deepcopy(node_set)
            segset_to_append["set"] = tris_global[element_indices, :]
            self.segment_sets.append(segset_to_append)

        return

    def _get_apex(self, nodes: np.array):
        """Extracts an apical point from the cavity

        Parameters
        ----------
        nodes : np.array
            Array with global node coordinates
        """
        LOGGER.debug("Extracting apical points")
        LOGGER.warning("No check whether detected point is at an edge")

        # get reference point (mid-point between two valve centroids)
        # robust enough to extract apical points for either left or right ventricle?
        centroids = np.empty((0, 3))
        for cap in self.closing_caps:
            centroids = np.vstack((centroids, cap.centroid))
        reference_point = np.mean(centroids, axis=0)

        # use the defined node sets to extract the apical point:
        for nodeset in self.node_sets:
            set_name = nodeset["name"]
            # no sense in trying to find the apex for the septum
            if "septum" in nodeset["name"]:
                continue

            points = nodes[nodeset["set"], :]
            id_max = np.argmax(np.linalg.norm(points - reference_point, axis=1))
            global_node_id = nodeset["set"][id_max]

            # validate:
            # checks whether node id is not on edge of segment set
            # if selected point is on edge of segment set select a point
            # which is not on a free edge (free face)
            from ansys.heart.preprocessor.vtk_module import get_free_edges

            for segset in self.segment_sets:
                triangles = segset["set"]
                if set_name == segset["name"]:
                    free_edges, free_triangles = get_free_edges(triangles, True)
                    if np.isin(global_node_id, free_edges):
                        LOGGER.warning("Apical point is on edge of segment set")
                        LOGGER.warning("Selecting point not on edge of segment set")
                        # select triangle that has two points on edge
                        free_edge = free_edges[
                            np.argwhere(np.all(np.isin(free_edges, global_node_id), axis=1)), :
                        ]
                        triangle_to_use = np.argwhere(
                            np.sum(np.isin(free_triangles, free_edge), axis=1) == 2
                        )

                        # select the node that is not used in the free edge
                        global_node_id = triangle_to_use[
                            np.isin(triangle_to_use, free_edge, invert=True)
                        ]
                    break

            self.apex_id[nodeset["name"]] = global_node_id

        return

    def compute_volume(self, stl_path: str):
        """Computes the volume of the enclosed surface"""
        # uses the inner surface of each cavity to compute the volume
        # vtk_surface = vtk.vtkPolyData
        self.volume = compute_volume_stl(stl_path)

        return

    def _dump_cap_nodes_to_file(self, working_dir: str):
        """Dumps the cap nodes to a file"""
        LOGGER.debug("Writing cap raw coordinates to file")
        for cap in self.closing_caps:
            filename = os.path.join(
                working_dir,
                "points_" + cap.name.replace(" ", "_").lower() + ".csv",
            )
            cap.dump_nodes_to_file(filename)

        return

    def _read_cap_nodes_from_file(self, working_dir: str):
        """Reads the nodes that define the caps from a file"""
        for cap in self.closing_caps:
            filename = os.path.join(working_dir, cap.name.lower().replace(" ", "_"))
            try:
                cap.nodes_source_mesh = np.genfromtxt(filename, delimiter=",")
            except:
                raise FileExistsError("Failed to read: " + filename)
        return


if __name__ == "__main__":
    print("Protected")
