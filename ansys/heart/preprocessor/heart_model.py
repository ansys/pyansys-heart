from multiprocessing.sharedctypes import Value
import os

from ansys.heart.preprocessor.heart_mesh import HeartMesh
from ansys.heart.preprocessor.model_information import ModelInformation

# import logger
from ansys.heart.custom_logging import LOGGER, Logger
from ansys.heart.preprocessor.vtk_module import (
    # get_edge_loops,
    get_tri_info_from_polydata,
    read_vtk_polydata_file,
    read_vtk_unstructuredgrid_file,
    threshold_vtk_data,
    threshold_vtk_data_integers,
    vtk_remove_arrays,
    write_vtkdata_to_vtkfile,
)


class HeartBoundaryConditions:
    """Class containing relevant functions for the boundary conditions"""

    def __init__(self):
        """Initializes object"""
        return


class HeartModel:
    """Class which represents the heart model of a specific type"""

    # some properties:

    def __init__(self, model_info: ModelInformation):
        """Initializes the heart model by making use of the specified
        model information object

        Parameters
        ----------
        model_info : ModelInformation
            Class which specifies relevant information
        """

        self.info = model_info

        self._mesh = HeartMesh(model_info)

        return

    def load_model(self, filename: str):
        """Loads all relevant model information for the json file

        Parameters
        ----------
        filename : str
            Path to model information in json format
        """

        # exposed to user
        LOGGER.info("Loading heart model from: %s " % filename)
        model_info = self._mesh.load_mesh(filename)
        self.info = model_info
        return

    def dump_model(self, filename: str = "model_info.json", clean_working_directory: bool = False):
        """Dumps model information for future use. Exports simulation mesh in .vtk format

        Parameters
        ----------
        filename : str, optional
            Path to model info json, by default "model_info.json"
        clean_working_directory : bool, optional
            Flag indicating whether to clean the working directory of any
            temporary files, by default False
        """

        # exposed to user
        if clean_working_directory:
            import glob as glob

            files = glob.glob(os.path.join(self.info.working_directory, "*"))

            if len(files) > 0:
                LOGGER.debug("Files detected: cleaning all files from directory")
            for f in files:
                os.remove(f)

        if not filename:
            filename = os.path.join(self.info.working_directory, "model_info.json")

        self.get_model_characteristics()

        # export the mesh
        self._mesh.export_mesh(filename)

        return

    def extract_simulation_mesh(self, remesh: bool = True):
        """Extracts the simulation mesh based on the model information provided

        Parameters
        ----------
        do_remesh : bool, optional
            Flag indicating whether to perform a remeshing step with the
            specified cell size, by default True
        """
        # exposed to user

        # perform actions on mesh object
        self._mesh.load_raw_mesh()
        self._mesh.extract_parts()
        self._mesh.add_cavities()
        self._mesh.get_cavity_cap_intersections()

        if remesh:
            mesh_size = self.info.mesh_size
            LOGGER.info("Remeshing uniformly with mesh size: %f " % mesh_size)
            self._mesh.remesh_volume(
                mesher="fluent", mesh_size=mesh_size
            )  # need to make this dynamic.
            self.info.is_remeshed = True

            # map data of original mesh to remeshed volume
            # includes (discrete) mapping of the tags
            self._mesh.map_data_to_remeshed_volume()
        else:
            myocardium_volume_vtk = self._mesh._extract_myocardium()
            self._mesh._vtk_volume = myocardium_volume_vtk
            self.info.is_remeshed = False
            self.info._mesh_size = None
            # raise ValueError("Using original input data for simulation " "not yet supported")
            LOGGER.warning("Extracting mesh without remeshing not thoroughly tested yet")

        # extract surface from remeshed volume mesh
        self._mesh.get_surface_from_volume()

        # closes the cavities
        self._mesh.close_cavities()

        # extract endo/epicardium surfaces
        self._mesh.extract_endocardium_epicardium()

        # extract apex
        self._mesh.extract_apical_points()

        # validate cavities
        self._mesh._validate_cavities()

        # extract volumetric region of septum
        if self.info.model_type in ["BiVentricle", "FourChamber"]:
            self._mesh._extract_septum()

        # create element sets for myocardium
        self._mesh._create_myocardium_element_sets()

        return

    def extract_simulation_mesh_improved(self, remesh: bool = True):
        """Extracts the simulation mesh based on the model information provided. Improved version

        Parameters
        ----------
        do_remesh : bool, optional
            Flag indicating whether to perform a remeshing step with the
            specified cell size, by default True
        """
        import copy
        from geomdl import fitting
        import numpy as np

        # model type determines what parts to extract
        from ansys.heart.preprocessor.vtk_module import (
            get_tetra_info_from_unstructgrid,
            write_vtkdata_to_vtkfile,
            vtk_surface_to_stl,
            create_vtk_surface_triangles,
            get_connected_regions,
            get_free_edges,
        )
        from ansys.heart.preprocessor.mesh_connectivity import (
            face_tetra_connectivity,
            get_face_type,
        )
        from ansys.heart.preprocessor.edge_module import edge_connectivity
        from ansys.heart.preprocessor.geodisc_module import project_3d_points

        write_folder = self.info.working_directory

        self._mesh.load_raw_mesh()
        self._mesh._vtk_volume_temp = self._mesh._vtk_volume_raw

        # extract parts
        self._mesh.extract_parts()
        self._mesh.add_cavities()
        self._mesh.get_cavity_cap_intersections()

        vtk_labels = copy.deepcopy(self.info.vtk_labels_to_use)
        # invert label map
        vtk_labels_inverse = dict((v, k) for k, v in vtk_labels.items())
        tags_to_use = list(vtk_labels.values())

        points, tetra, cell_data, point_data = get_tetra_info_from_unstructgrid(
            self._mesh._vtk_volume_raw, deep_copy=True
        )

        original_points = points

        # remove any redundant cells and celldata
        mask = np.isin(cell_data["tags"], tags_to_use)
        tetra = tetra[mask, :]
        for key, value in cell_data.items():
            cell_data[key] = cell_data[key][mask]

        # get part ids
        part_ids = cell_data["tags"]
        part_ids = np.array(part_ids, dtype=int)

        # get face-tetra connectivity
        faces, c0, c1 = face_tetra_connectivity(tetra)

        # get face types
        # face_type == 1: interior face
        # face_type == 2: boundary face
        face_types = get_face_type(faces, np.array([c0, c1]).T)

        mask_interface_faces = np.all(
            np.vstack([part_ids[c0] != part_ids[c1], face_types == 1]), axis=0
        )

        # get number of unique interfaces:
        interface_pairs = np.vstack(
            [part_ids[c0][mask_interface_faces], part_ids[c1][mask_interface_faces]]
        )
        # get all tag pairs
        interface_pairs_sorted = np.sort(interface_pairs, axis=0)
        tag_pairs = np.unique(interface_pairs_sorted, axis=1).transpose()

        # remove myocardium<>myocardium interfaces
        # remove myocardium<>aorta wall and myocardium<>pulmonary artery wall interfaces
        tag_pairs_1 = np.empty((0, 2), dtype=int)
        for tag_pair in tag_pairs:
            labels = [vtk_labels_inverse[tag_pair[0]], vtk_labels_inverse[tag_pair[1]]]
            res1 = ["myocardium" in e for e in labels]
            res2 = ["wall" in e for e in labels]
            if np.all(["myocardium" in e for e in labels]):
                # skip if myocardium-myocardium intersection
                LOGGER.debug(
                    "Skipping {0} | {1} | {2}".format(
                        *tag_pair, "interface myocardium and myocardium"
                    )
                )
                continue
            if (
                np.sum(["myocardium" in e for e in labels]) == 1
                and np.sum(["wall" in e for e in labels]) == 1
            ):
                # skip if myocardium<>pulmonary artery or myocardium<>aorta wall interfaces
                LOGGER.debug(
                    "Skipping {0} | {1} | {2}".format(*tag_pair, "interface myocardium and artery")
                )
                continue
            tag_pairs_1 = np.vstack([tag_pairs_1, tag_pair])
        tag_pairs = tag_pairs_1

        # for labels_pair in label_pairs:
        interfaces = []
        for pair in tag_pairs:
            # skip any unused parts:
            if not np.all(np.isin(pair, tags_to_use)):
                LOGGER.debug("Skipping pair: {0}|{1}".format(*pair))

            labels_pair = [vtk_labels_inverse[pair[0]], vtk_labels_inverse[pair[1]]]
            pair_name = "_".join(labels_pair).replace(" ", "-").lower()
            tags_pair = [vtk_labels[labels_pair[0]], vtk_labels[labels_pair[1]]]

            mask_part1 = (
                np.sum(
                    np.vstack([part_ids[c0] == tags_pair[0], part_ids[c1] == tags_pair[0]]), axis=0
                )
                == 1
            )
            mask_part2 = (
                np.sum(
                    np.vstack([part_ids[c0] == tags_pair[1], part_ids[c1] == tags_pair[1]]), axis=0
                )
                == 1
            )

            # get interface between pairs
            face_ids_interface = np.where(
                np.all(np.vstack([mask_interface_faces, mask_part1, mask_part2]), axis=0)
            )[0]

            faces_interface = faces[face_ids_interface, :]

            # get free edges
            free_edges, _ = get_free_edges(faces_interface, return_free_triangles=True)
            edge_groups, edge_group_types = edge_connectivity(
                free_edges, return_type=True, sort_closed=True
            )

            # print some info
            num_open_loops = np.sum(np.array(edge_group_types) == "open")
            LOGGER.debug("Pair: {0}".format(pair_name))
            LOGGER.debug(
                "Number of edge loops: {0} | open_ended: {1}".format(
                    len(edge_groups), num_open_loops
                )
            )

            # NOTE: do some modifications on edge loops that are not closed?
            if len(edge_groups) > 2:
                LOGGER.debug("Detected {0} edge groups for {1}".format(len(edge_groups), pair_name))
                lengths = []
                mask_closed = np.array(edge_group_types) == "closed"
                if np.sum(mask_closed) == 2:
                    LOGGER.debug("Could improve by selecting closed edge loops for this interface")
                for ii, group in enumerate(edge_groups):
                    lengths.append(len(group))
                    ids_closed = []

            # improve free edges of each interface
            use_geomdl = False
            use_average = True
            if len(edge_groups) == 2:
                for group_idx, edge_group in enumerate(edge_groups):
                    # keep order of unique points:
                    _, idx = np.unique(edge_group.flatten(), return_index=True)
                    node_ids = edge_group.flatten()[np.sort(idx)]

                    projected_points = project_3d_points(points[node_ids, :])[0]

                    # use connected points to compute some kind of walking/windowed mean
                    #  mean( p-1, pn, p+1 )
                    # only do for closed edge loops
                    # perform b-spline interpolation:
                    # https://nurbs-python.readthedocs.io/en/latest/module_fitting.html
                    # from scipy import interpolate
                    # x = projected_points[:,0]
                    # y = projected_points[:,1]
                    # z = projected_points[:,2]
                    # tck = interpolate.splrep(x,y, s=0, k=3)

                    if use_average:
                        num_iters = 2
                        if edge_group_types[group_idx] == "closed":
                            LOGGER.debug("Smoothing points...")
                            tmp = np.append(np.array(node_ids[-1]), node_ids[:])
                            tmp = np.append(tmp, node_ids[0])
                            for iter in range(0, num_iters, 1):
                                for ii in range(1, len(node_ids), 1):
                                    projected_points[ii - 1] = (
                                        points[tmp[ii - 1], :]
                                        + points[tmp[ii], :]
                                        + points[tmp[ii + 1], :]
                                    ) / 3

                    # for open-looped edge loops use different smoothing method
                    if use_geomdl:
                        LOGGER.debug("Using geomdl to fit curve")
                        num_iters = 1
                        iters = 0
                        while iters < num_iters:
                            curve = fitting.interpolate_curve(projected_points.tolist(), degree=2)
                            # curve.delta = 0.01
                            # curve.vis = vis.VisCurve3D()
                            # curve.render()
                            curve.sample_size = len(projected_points)
                            projected_points = np.array(curve.evalpts)
                            iters = iters + 1

                    # assign new points
                    points[node_ids, :] = projected_points

            # collect interface info
            interfaces.append(
                {
                    "name": pair_name,
                    "ids": pair,
                    "faces": faces_interface,
                    "boundary_edges": free_edges,
                    "edgegroups": [edge_groups, edge_group_types],
                }
            )

        # write data for each interface
        for interface in interfaces:
            # Write stl per interface
            vtk_surface = create_vtk_surface_triangles(points, interface["faces"])
            filename = os.path.join(
                write_folder,
                "part_interface_{:02d}_{:02d}.stl".format(interface["ids"][0], interface["ids"][1]),
            )
            vtk_surface_to_stl(vtk_surface, filename, solid_name=interface["name"])

            # write edges of interface in vtk
            tris_edge = np.vstack(
                [interface["boundary_edges"].T, interface["boundary_edges"][:, 1]]
            ).T
            tris_edge_vtk = create_vtk_surface_triangles(points, tris_edge)
            write_vtkdata_to_vtkfile(
                tris_edge_vtk,
                os.path.join(
                    self.info.working_directory,
                    "free_edges_pair_{0}_{1}.vtk".format(*interface["ids"]),
                ),
            )

        # exclude valves/inlets from extraction
        part_ids_to_extract = []
        for label, pid in vtk_labels.items():
            if "inlet" in label or "valve" in label:
                continue
            part_ids_to_extract.append(pid)
        part_ids_to_extract = np.sort(part_ids_to_extract)

        # write stl for each part
        for part_id in part_ids_to_extract:
            mask = np.vstack([face_types == 2, part_ids[c0] == part_id])
            mask = np.all(mask, axis=0)
            boundary_faces = faces[mask, :]

            # extract regions for myocardial parts (LV: 1, RV: 2, LA: 3, RA: 4)
            if "myocardium" in vtk_labels_inverse[part_id]:
                region_ids, vtk_surface = get_connected_regions(points, boundary_faces, True)
                # write_vtkdata_to_vtkfile(
                #     vtk_surface, "..\\tmp\\part_with_region_ids_{:02d}.vtk".format(part_id)
                # )

                # smalles region endocardium, largest region epicardium. In case of Left ventricle
                # smallest region is septum
                name_map = ["endocardium", "epicardium", "septum"]
                unique_region_ids, counts = np.unique(region_ids, return_counts=True)
                for ii, index in enumerate(np.argsort(counts)):
                    unique_region_ids[index]
                    vtk_surface1 = create_vtk_surface_triangles(
                        points, boundary_faces[region_ids == unique_region_ids[index]]
                    )
                    part_name = "-".join(
                        "{0}-{1}".format(vtk_labels_inverse[part_id], name_map[index])
                        .lower()
                        .split()
                    )
                    filename = os.path.join(write_folder, "part_{0}.stl".format(part_name))
                    vtk_surface_to_stl(vtk_surface1, filename, solid_name=part_name)
            else:
                vtk_surface = create_vtk_surface_triangles(points, boundary_faces)
                # write_vtkdata_to_vtkfile(vtk_surface, "..\\tmp\\part_{:02d}.vtk".format(part_id))
                part_name = "-".join("{0}".format(vtk_labels_inverse[part_id]).lower().split())
                filename = os.path.join(write_folder, "part_{0}.stl".format(part_name))
                vtk_surface_to_stl(vtk_surface, filename, solid_name=part_name)

        # start meshing:

        from ansys.heart.preprocessor.mesh_module import mesh_by_fluentmeshing

        LOGGER.debug("Remeshing volume")
        path_to_mesh = os.path.join(self.info.working_directory, "fluent_volume_mesh.msh.h5")
        mesh_by_fluentmeshing(
            filename,
            path_to_output=path_to_mesh,
            mesh_size=1.5,
            journal_type="simplified_geometry",
            show_gui=False,
        )

        # read mesh and extract surfaces
        from ansys.heart.preprocessor.fluenthdf5_module import fluenthdf5_to_vtk

        filename_simulation_mesh_vtk = os.path.join(
            os.path.dirname(path_to_mesh), "simulation_mesh.vtk"
        )

        tetrahedrons, face_zones, points = fluenthdf5_to_vtk(
            path_to_mesh, filename_simulation_mesh_vtk
        )
        # flip all the face zones:
        for face_zone_name, face_zone in face_zones.items():
            LOGGER.warning("Flipping face zone normals of %s " % face_zone_name)
            face_zone["faces"] = face_zone["faces"][:, [0, 2, 1]]

        # read mesh file
        vtk_volume = read_vtk_unstructuredgrid_file(filename_simulation_mesh_vtk)

        # # extract part labels. ignore valves and inlets
        # part_labels = {}
        # for vtk_label in vtk_labels.keys():
        #     if "valve" in vtk_label or "inlet" in vtk_label:
        #         continue
        #     part_labels[vtk_label] = vtk_labels[vtk_label]

        # # extract surface of each part (used to map the volume elements)
        # tetrahedrons_tags = np.zeros(tetrahedrons.shape[0], dtype=int) - 1
        # visited = tetrahedrons_tags != -1
        # tetrahedron_centroids = np.mean(points[tetrahedrons, :], axis=1)

        # tag volume elements with appropiate tag:
        # two methods exist:
        # 1. create high quality surface of individual parts
        # 2. interpolate tags from original mesh onto new mesh

        # map data to remeshed volume:
        # self._mesh._vtk_volume_temp
        self._mesh._vtk_volume = vtk_volume
        self._mesh.map_data_to_remeshed_volume()

        # write_vtkdata_to_vtkfile(
        #     self._mesh._vtk_volume, os.path.join(self.info.working_directory, "mapped_volume.vtk")
        # )

        # assign face zones to cavity segment sets
        # this already includes the endo- and epicardial segments
        for cavity in self._mesh._cavities:
            for name in ["endocardium", "epicardium"]:

                for face_zone_name, face_zone in face_zones.items():
                    if (
                        name in face_zone_name
                        and "-".join(cavity.name.lower().split()) in face_zone_name
                    ):
                        # print("Match: %s " % face_zone_name )
                        cavity.segment_sets.append(
                            {"name": name, "set": face_zone["faces"], "id": face_zone["zone-id"]}
                        )

            if cavity.name == "Left ventricle" and "LeftVentricle" not in self.info.model_type:
                cavity.segment_sets.append(
                    {
                        "name": "epicardium-septum",
                        "set": face_zones["left-ventricle-myocardium-septum"]["faces"],
                        "id": face_zones["left-ventricle-myocardium-septum"]["zone-id"],
                    }
                )

        # get node sets from segment sets
        for cavity in self._mesh._cavities:
            # remove any defined node-sets
            LOGGER.debug("Removing %d nodesets" % len(cavity.node_sets))
            cavity.node_sets = []
            # NOTE: this duplicates the nodes on the intersection of the
            for segset in cavity.segment_sets:
                nodeset = copy.deepcopy(segset)
                nodeset["set"] = np.unique(segset["set"])
                cavity.node_sets.append(nodeset)

        # validate node and segment sets:
        self._mesh._validate_node_sets()
        # self._mesh._validate_segment_sets()  # obsolete

        # find edges of the valve
        for cavity in self._mesh._cavities:
            # find edge loops of endocardial parts
            for segset in cavity.segment_sets:
                if segset["name"] == "endocardium":
                    break
            faces_endocardium = segset["set"]

            # find valves - determined by label pairs (e.g. right ventricle endocardium <> tricuspid valve plane)
            label_myocardium = cavity.labels[0]
            pairs = []
            for label in cavity.labels[1:]:
                pairs.append([label_myocardium, label])

            zone_names = []
            for pair in pairs:
                # match name of face zone from Fluent
                zone_name = "-".join("_".join(pair).lower().split())
                zone_names.append(zone_name)

            # find edge loops for each cap
            for cap in cavity.closing_caps:
                zone_name = "-".join(
                    "_".join([cavity.name + " myocardium", cap.name]).lower().split()
                )
                faces = face_zones[zone_name]["faces"]
                edges = get_free_edges(faces)
                edge_loops, loop_types = edge_connectivity(
                    edges, return_type=True, sort_closed=True
                )
                # find edge loop that is connected to the endocardium
                if np.any(np.array(loop_types) != "closed"):
                    LOGGER.error("Expecting all edge loops to be closed")
                    raise ValueError("Expecting all edge loops to be closed")

                in_endocardium = False
                for edge_loop in edge_loops:
                    nodes_edge = edge_loop[:, 0]  # order is important
                    mask = np.isin(nodes_edge, faces_endocardium)
                    if np.all(mask):
                        in_endocardium = True
                        break
                if not in_endocardium:
                    LOGGER.error("Expecting nodes to be connected to endocardium")
                    raise ValueError("Expecting nodes to be connected to endocardium")

                cap.node_ids_cap_edge = nodes_edge
                cap.nodes_cap_edge = points[nodes_edge, :]
                cap.centroid = np.mean(points[nodes_edge, :], axis=0)

                # triangulate the caps
                # cavity._triangulate_caps(self._mesh._vtk_volume)
        for cavity in self._mesh._cavities:
            for cap in cavity.closing_caps:
                # use three representative points to compute cap normal
                # assumes these are ordered either clockwise or counterclockwise
                # check whether normal is pointing inwards - if not > invert
                offset = int(np.floor(len(cap.nodes_cap_edge) / 3))
                normal = np.cross(
                    cap.nodes_cap_edge[offset] - cap.nodes_cap_edge[0],
                    cap.nodes_cap_edge[offset * 2] - cap.nodes_cap_edge[0],
                )
                cap.normal = normal / np.linalg.norm(normal)
                d1 = np.linalg.norm((cap.centroid + cap.normal) - cavity.centroid)
                d2 = np.linalg.norm((cap.centroid - cap.normal) - cavity.centroid)
                # flip order if d1>d2
                if d1 > d2:
                    LOGGER.debug("Flipping normal of cap: %s" % cap.name)
                    cap.normal = cap.normal * -1
                    cap.nodes_cap_edge = np.flip(cap.nodes_cap_edge)
                    cap.node_ids_cap_edge = np.flip(cap.node_ids_cap_edge)

            cavity._triangulate_caps(self._mesh._vtk_volume)

        # extract apex
        self._mesh.extract_apical_points()

        self._mesh._validate_segment_sets(check_duplicates=False)

        # validate cavities
        self._mesh._validate_cavities()

        # extract volumetric region of septum
        if self.info.model_type in [
            "BiVentricleImproved",
            "FullHeartImproved",
            "FourChamberImproved",
        ]:
            self._mesh._extract_septum()

        # create element sets for myocardium
        self._mesh._create_myocardium_element_sets()

        return

    def extract_simulation_mesh_from_simplified_geometry(self):
        """Extracts a simulation mesh from a modified/simplified geometry
        from Strocchi or Cristobal et al
        """
        from ansys.heart.preprocessor.vtk_module import (
            read_vtk_polydata_file,
        )
        import numpy as np
        import copy

        # node-tag mapping:
        cavity_tag_map = {"Left ventricle": 1, "Right ventricle": 2}
        node_tag_map = {
            "Left ventricle": {
                "epicardium": 0,
                "mitral-valve-edge": 1,
                "aortic-valve-edge": 2,
                "endocardium": 3,
            },
            "Right ventricle": {
                "epicardium": 4,
                "pulmonary-valve-edge": 5,
                "tricuspid-valve-edge": 6,
                "interventricular-edge": 7,
                "endocardium": 8,
            },
        }

        # read surface mesh
        self._mesh._vtk_volume_raw = read_vtk_polydata_file(self.info.path_original_mesh)

        self._mesh._vtk_volume_temp = self._mesh._vtk_volume_raw

        # convert to PolyData
        self._mesh._vtk_surface = (
            self._mesh._vtk_volume_raw
        )  # = convert_to_polydata(self._mesh._vtk_volume_raw)

        # add cavities
        self._mesh.add_cavities()

        self._mesh.get_cavity_cap_intersections_simplified(node_tag_map)

        self._mesh.mesh_volume_from_simplified(node_tag_map, self.info.mesh_size)

        # get node sets from segment sets
        for cavity in self._mesh._cavities:
            # remove any defined node-sets
            LOGGER.debug("Removing %d nodesets" % len(cavity.node_sets))
            cavity.node_sets = []
            # NOTE: this duplicates the nodes on the intersection of the
            for segset in cavity.segment_sets:
                nodeset = copy.deepcopy(segset)
                nodeset["set"] = np.unique(segset["set"])
                cavity.node_sets.append(nodeset)

        self._mesh._validate_node_sets()
        self._mesh._validate_segment_sets()

        # extract the surface from the remeshed volume mesh
        self._mesh.get_surface_from_volume()

        # closes the cavities
        self._mesh.close_cavities()

        # extract apical points
        self._mesh.extract_apical_points()

        # validate cavities
        self._mesh._validate_cavities()

        # extract volumetric region of septum
        if self.info.model_type in ["BiVentricle", "FourChamber"]:
            self._mesh._extract_septum()

        # create element sets for myocardium
        self._mesh._create_myocardium_element_sets()

        return

    def get_model_characteristics(self, write_to_file: bool = True) -> dict:
        """Creates a dictionary of model characteristics

        Parameters
        ----------
        write_to_file : bool, optional
            Flag indicating whether to write to file , by default True

        Returns
        -------
        dict
            Dictionary of model characteristics
        """
        characteristics = {"model": {}, "mesh": {"cavity": []}}

        characteristics["model"]["type"] = self.info.model_type
        characteristics["mesh"]["num-cavities"] = len(self._mesh._cavities)
        num_caps = 0
        cavity_names = []
        for cavity in self._mesh._cavities:
            cavity_names.append(cavity.name)
            cap_names = []
            for cap in cavity.closing_caps:
                cap_names.append(cap.name)
                num_caps += 1
            cavity.volume
            cavity_info = {
                "name": cavity.name,
                "cap-names": cap_names,
                "volume": cavity.volume
                # "caps"      : cavity.closing_caps
            }
            characteristics["mesh"]["cavity"].append(cavity_info)

        characteristics["mesh"]["num-caps"] = num_caps
        characteristics["mesh"]["cavity-names"] = cavity_names

        # write to json
        import json

        json_path = os.path.join(self.info.working_directory, "model_characteristics.json")
        with open(json_path, "w") as outfile:
            json.dump(characteristics, indent=4, fp=outfile)

        return characteristics


if __name__ == "__main__":

    print("Protected")
