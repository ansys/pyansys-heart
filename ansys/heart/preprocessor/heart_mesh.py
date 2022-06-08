"""Contains classes related to the heart mesh class"""
import os
from tkinter import SOLID
import vtk
import copy

import pickle as pickle

from multiprocessing.sharedctypes import Value
from pickletools import read_unicodestring1
from typing import List

from ansys.heart.preprocessor.cavity_module import *
from ansys.heart.preprocessor.model_information import ModelInformation
from ansys.heart.preprocessor.model_information import VALID_MODELS
from ansys.heart.preprocessor.vtk_module import (
    create_vtk_surface_triangles,
    write_vtkdata_to_vtkfile,
)

# import logger
from ansys.heart.custom_logging import logger

from vtk.numpy_interface import (
    dataset_adapter as dsa,
)  # this is an improved numpy integration


class HeartMesh:
    """Class which represents the heart mesh - contains functions to read/write
    the various parts and performs actions such as remeshing, extracting surfaces, etc """

    def __init__(self, model_info: ModelInformation):
        """Initializes heart mesh """
        # TODO What other actions to perform on init?
        # TODO Move some of the methods to seperate public functions

        self.info = model_info

        return

    # public functions
    def load_raw_mesh(self):
        """Reads the raw input mesh
        """
        self._vtk_volume_raw = vtk_read_mesh_file(self.info.path_original_mesh)
        return

    def load_mesh(self, model_info_path: str):
        """Reads the mesh from information provided by model info """

        self.info = ModelInformation()
        self.info.load_from_file(model_info_path)

        # read mesh from vtk file
        self._vtk_volume = vtk_read_mesh_file(self.info.path_mesh)

        # read cavities from pickle
        fid = open(self.info.path_to_pickle, "rb")
        self._cavities = pickle.load(fid)
        fid.close()
        # read cavity._surfaces, and cavity._myocardium?
        return self.info

    def export_mesh(self, model_info_path: str):
        """Exports the mesh and all (relevant) generated data
        for future use
        """
        path_simulation_mesh = os.path.join(
            self.info.working_directory, "simulation_mesh.vtk"
        )

        # save processed volume mesh in designated path
        write_vtkdata_to_vtkfile(self._vtk_volume, path_simulation_mesh)
        # dump
        pickle_path = os.path.join(
            self.info.working_directory, "cavities.pickle"
        )
        fid = open(pickle_path, "wb")

        # remove any vtk objects which cannot be pickled
        # NOTE: should these be saved as vtk instead?
        for cavity in self._cavities:
            del cavity._surfaces
            del cavity._myocardium_surface
        pickle.dump(self._cavities, fid)
        fid.close()

        self.info.path_to_pickle = pickle_path
        self.info.path_mesh = path_simulation_mesh

        self.info.save_to_file(model_info_path)

        return

    def extract_parts(self):
        """Extracts the relevant parts from the mesh given the vtk labels/tags
        to use. Also extracts info such as cell data and point data from the vtk files"""

        write_temporary_files = False
        # Add tags where elements overlap
        # self.separate_overlaps( self._vtk_volume_ugriddata_raw )
        # self.separate_overlaps_surface( self._vtk_volume_ugriddata_raw )

        # extract only the parts that are of interest
        self._vtk_volume_temp = self._extract_relevant_parts(
            self._vtk_volume_raw, self.info.vtk_labels_to_use
        )

        if write_temporary_files:
            write_vtkdata_to_vtkfile(self._vtk_volume_temp, "temp.vtk")
        # self.extract_individual_parts_surface()

        # extract surface data
        self._vtk_surface_polydata = vtk_surface_filter(
            self._vtk_volume_temp, keep_global_ids=True
        )

        if write_temporary_files:
            write_vtkdata_to_vtkfile(
                self._vtk_surface_polydata, "temp_surface.vtk"
            )

        # useful in some calls: NOTE: Just to ensure compatibility: might be removed?
        self._vtk_surface_data = get_surface_info(self._vtk_surface_polydata)

        # replaces previous line eventually?
        (
            nodes_tris,
            tris,
            cell_data_tris,
            point_data_tris,
        ) = get_tri_info_from_polydata(self._vtk_surface_polydata)
        self._nodes_surface = nodes_tris
        self._elements_surface = tris
        self._cell_data_surface = cell_data_tris
        self._point_data_surface = point_data_tris

        # extract nodes and tetra definitions from raw data
        nodes, tetra, cell_data, point_data = get_tetra_info_from_unstructgrid(
            self._vtk_volume_raw
        )
        self._nodes_volume_raw = nodes
        self._elements_volume_raw = tetra
        self._cell_data_volume_raw = cell_data
        self._point_data_volume_raw = point_data

        # extract nodes and tetra definitions from thresholded data

        # NOTE: This is redundant?
        # nodes, tetra, cell_data, point_data = get_tetra_info_from_unstructgrid( self._vtk_volume_temp )

        # self._nodes_volume      = nodes
        # self._elements_volume   = tetra
        # self._cell_data_volume  = cell_data
        # self._point_data_volume = point_data

        # writes each part to seperate stl file
        # self._write_parts_to_stl()

        return

    def add_cavities(self):
        """Adds cavities based on the requested model"""
        dump_to_file = False

        # get list of involved cavities
        self._cavities = self._add_list_of_cavities()

        # add vtk surface object of all relevant parts to each cavity
        for cavity in self._cavities:
            # cavity._set_connected_surfaces( self._vtk_surface_polydata )
            # NOTE: now uses raw surface mesh. This may not be very efficient
            cavity._set_connected_surfaces(self._vtk_volume_temp)
            logger.debug("Processing... " + cavity.name)

            logger.warning("Computing cavity centroid")
            cavity._compute_centroid_cavity()

            # point = create_vtk_polydata_from_points( np.array([cavity.centroid]) )
            # write_vtkdata_to_vtkfile(point, "point.vtk")

            logger.debug("Computing intersection with myocardium")
            cavity._compute_cap_intersection_with_myocardium(
                self._nodes_volume_raw,
                self._elements_volume_raw,
                self._cell_data_volume_raw["tags"],
            )

        # validate cavities:
        # for cavity in self._cavities:
        #     print(cavity.name)
        #     for cap in cavity.closing_caps:
        #         logger.debug(cap.name)

        # write nodes that close the cap to file
        if dump_to_file:
            for cavity in self._cavities:
                cavity._dump_cap_nodes_to_file(self.info.working_directory)

        # NOTE: this uses another object - node numbering should be consistent with previous
        # can lead to issues?

        return

    def close_cavities(self):
        """Finds the edge loops that close the cavities
        e.g. at the aortic, mitral, tricuspid, and pulmoanry valves. """
        # vtk_volume = dsa.WrapDataObject( self._vtk_volume )
        write_caps_to_file = True

        for cavity in self._cavities:
            logger.debug("Closing cavity... %s" % cavity.name)
            cavity._find_closing_edge_loop()

            # TODO create shell elements by using the closing edge loop
            cavity._triangulate_caps(self._vtk_volume)

            # write caps of each cavity to stl file:
            if write_caps_to_file:
                all_tris = np.empty((0, 3))
                volume = dsa.WrapDataObject(self._vtk_volume)

                for cap in cavity.closing_caps:
                    all_tris = np.vstack((all_tris, cap.closing_triangles))

                logger.debug("Writing caps of %s" % cavity.name)

                mesh = meshio.Mesh(
                    points=volume.Points,
                    cells=[("triangle", np.array(all_tris, dtype=int))],
                )

                filename = "caps_" + "_".join(cavity.name.split()) + ".stl"
                file_path = os.path.join(self.info.working_directory, filename)
                mesh.write(file_path)

            # for cap in cavity.closing_caps:
            #     vtk_pointcloud = create_vtk_polydata_from_points( cap.nodes_cap_edge )
            #     write_vtkdata_to_vtkfile( vtk_pointcloud,
            #                                 "points_" + "_".join( cavity.name.split() ) + "_"
            #                                 + "_".join(cap.name.split()).lower() + ".vtk" )

        return

    def remesh_volume(self, mesher: str = "gmesh", mesh_size: float = 2.0):
        """Remeshes the volume. Uses the surface as input, uses SpaceClaim to smooth the 
        geometry and gmesh to remesh the bulk"""

        logger.debug("Remeshing volume...")

        if mesher == "gmesh":
            use_gmesh = True
            use_fluent = False
        elif mesher == "fluent":
            use_fluent = True
            use_gmesh = False
        else:
            raise ValueError("Specified mesher not recognized")

        working_dir = self.info.working_directory

        # only extract myocardium parts
        vtk_tags_to_use = []
        vtk_tags_to_use1 = {}
        for cavity in self._cavities:
            for ii, label in enumerate(cavity.labels):
                if "myocardium" in label:
                    vtk_tags_to_use.append(cavity.vtk_ids[ii])
                    vtk_tags_to_use1[label] = cavity.vtk_ids[ii]

        spaceclaim_input = self._extract_relevant_parts(
            self._vtk_volume_temp, vtk_tags_to_use1
        )

        input_stl_spaceclaim = os.path.join(
            working_dir, "input_surface_spaceclaim.stl"
        )
        output_stl_spaceclaim = os.path.join(
            working_dir, "output_surface_spaceclaim.stl"
        )

        input_stl_fluent = os.path.abspath(
            os.path.join(working_dir, "input_surface_spaceclaim.stl")
        )
        output_msh_fluent = os.path.abspath(
            os.path.join(working_dir, "fluent_volume_mesh.msh.h5")
        )

        output_vtk = os.path.join(working_dir, "output_meshing.vtk")

        vtk_surface_to_stl(spaceclaim_input, input_stl_spaceclaim)

        if use_gmesh:
            # launch spaceclaim to wrap the surface
            logger.debug("\tLaunching SpaceClaim...")
            shrink_by_spaceclaim(input_stl_spaceclaim, output_stl_spaceclaim)

            logger.debug("\tLaunching GMESH...")
            run_gmsh(output_stl_spaceclaim, output_vtk, mesh_size)

            # extract tetrahedrons
            convert_vtk_into_tetra_only(output_vtk)

        elif use_fluent:
            logger.debug("\tLaunching Fluent meshing...")
            mesh_by_fluentmeshing(
                input_stl_fluent, output_msh_fluent, mesh_size
            )

            fluenthdf5_to_vtk(output_msh_fluent, output_vtk)

        self.info.mesh_size = mesh_size
        # store path in info
        self.info.path_mesh = output_vtk

        # cleanup        return
        return

    def map_data_to_remeshed_volume(self):
        """Maps the data from the original mesh file to 
        the remeshed vtk file 
        """
        logger.debug("Mapping data from original mesh...")

        # NOTE Source contains all point data and cell data fields
        # NOTE Map the tags of only the myocardium

        source = self._vtk_volume_temp
        target = self._vtk_volume

        target = vtk_map_continuous_data(source, target)

        # for debugging purposes
        # write_vtkdata_to_vtkfile( source, "source_interp.vtk" )
        # write_vtkdata_to_vtkfile( target, "target_interp_1.vtk" )

        # map the tags (consider this to be discrete data)
        # store in as _vtk_volume
        self._vtk_volume = self._map_tags_to_remeshed_volume(source, target)

        return

    def extract_endocardium_epicardium(self):
        """Extracts the endo and epicardium from each of myocardial parts.
        """
        # extract surface
        surface_mesh = vtk_surface_filter(
            self._vtk_volume, keep_global_ids=True
        )

        # get the endo and epicardium points for each cavity
        for cavity in self._cavities:
            logger.debug(
                "Extracting endocardium and epicardium node sets of... "
                + cavity.name
            )
            cavity._get_endocardium_epicardium_points(surface_mesh)

        # extract endo/epicardium segments for each cavity
        for cavity in self._cavities:
            logger.debug(
                "Extracting endocardium and epicardium segments of... "
                + cavity.name
            )
            cavity._get_endocardium_epicardium_segments(surface_mesh)

            # write to file to check
            volume = dsa.WrapDataObject(self._vtk_volume)
            points = volume.Points
            cells = []
            celldata = {"id": []}
            ii = 0
            for seg_set in cavity.segment_sets:
                cells.append(("triangle", seg_set["set"]))
                celldata["id"].append(np.ones(seg_set["set"].shape[0]) * ii)
                ii = ii + 1

            mesh = meshio.Mesh(points=points, cells=cells, cell_data=celldata)
            mesh.write(
                os.path.join(
                    self.info.working_directory,
                    "segset_" + "_".join(cavity.name.lower().split()) + ".vtk",
                )
            )

        # validate the nodesets and segment sets to remove any duplicate nodes/segments
        self._validate_node_sets()
        self._validate_segment_sets()

        # write segment sets to file
        for cavity in self._cavities:
            # write to file to check
            volume = dsa.WrapDataObject(self._vtk_volume)
            points = volume.Points
            cells = []
            celldata = {"id": []}
            ii = 0
            for seg_set in cavity.segment_sets:
                cells.append(("triangle", seg_set["set"]))
                celldata["id"].append(np.ones(seg_set["set"].shape[0]) * ii)
                ii = ii + 1

            mesh = meshio.Mesh(points=points, cells=cells, cell_data=celldata)
            mesh.write(
                os.path.join(
                    self.info.working_directory,
                    "segset_" + "_".join(cavity.name.lower().split()) + ".vtk",
                )
            )

        return

    def extract_apical_points(self):
        """Extracts apical points from each ventricle"""

        # get nodes
        volume_dsa = dsa.WrapDataObject(self._vtk_volume)
        nodes = volume_dsa.Points
        
        for cavity in self._cavities:
            # skip if atrium
            if "atrium" in cavity.name:
                continue

            cavity._get_apex(nodes)

        return

    def _validate_node_sets(self):
        """Validates the node sets and makes sure there are no references to duplicate nodes. 
        Endocardium of any cavity takes precedence over septum for instance.
        """
        logger.debug("Validating node sets...")
        nodes_used = []
        part_id = []
        for cavity in self._cavities:
            for nodeset in cavity.node_sets:
                for ii, nodes in enumerate(nodes_used):
                    # check if any of the nodes are already in use
                    use_idx = np.argwhere(
                        np.isin(nodeset["set"], nodes, invert=True)
                    ).flatten()
                    if len(use_idx) != len(nodeset["set"]):
                        logger.warning(
                            "Found duplicate nodes in: {0} {1}. Already used in {2} {3}".format(
                                cavity.name,
                                nodeset["name"],
                                part_id[ii][0],
                                part_id[ii][1],
                            )
                        )

                    nodeset["set"] = nodeset["set"][use_idx]

                part_id.append([cavity.name, nodeset["name"]])

                nodes_used.append(nodeset["set"])

            # write points of cavity to working directory
        for cavity in self._cavities:
            volume = dsa.WrapDataObject(self._vtk_volume)
            for nodeset in cavity.node_sets:

                vtk_poly_points = create_vtk_polydata_from_points(
                    volume.Points[nodeset["set"]]
                )

                path_to_file = os.path.join(
                    self.info.working_directory,
                    "points_{0}_{1}.vtk".format(
                        nodeset["name"], "_".join(cavity.name.lower().split())
                    ),
                )

                write_vtkdata_to_vtkfile(vtk_poly_points, path_to_file)

        return

    def _validate_segment_sets(self):
        """Validates the segment sets. 1. Checks whether to add the epicardium-septum of
        the left ventricle to the right-ventricle cavity. 2. Checks if normal 
        of endocardium part is pointing inwards
        """
        logger.debug("Validating segment sets...")
        # add epicardium-septum to right ventricle segment set
        if self.info.model_type in ["BiVentricle", "FourChamber"]:
            segset_septum = None
            for cavity in self._cavities:
                if cavity.name == "Left ventricle":
                    for ii, segset in enumerate(cavity.segment_sets):
                        if segset["name"] == "epicardium-septum":
                            segset_septum = copy.deepcopy(segset)
                    # remove key from segment set list
                    cavity.segment_sets.remove(segset)

            if segset_septum is None:
                raise Error(
                    "Did not find segment set in Left ventricle cavity "
                    "with name epicardium-septum"
                )

            for cavity in self._cavities:
                if cavity.name == "Right ventricle":
                    segset_septum["name"] = "endocardium-septum"
                    cavity.segment_sets.append(segset_septum)
                    logger.debug(
                        "Assigning septum of left ventricle to "
                        "segment set to Right ventricle"
                    )

        # remove any duplicate segments within each cavity
        logger.debug("Detecting duplicate segments...")
        all_tris = np.empty((0, 3), dtype=int)
        for cavity in self._cavities:

            for seg_set in cavity.segment_sets:
                np.sort(all_tris, axis=1)
                seg_set1 = np.sort(seg_set["set"], axis=1)

                # find triangles to remove
                idx_remove = []
                for idx, tri in enumerate(seg_set1):
                    if np.any(np.all(tri == all_tris, axis=1)):
                        idx_remove.append(idx)
                logger.debug(
                    "Found {0} duplicate segments in {1}".format(
                        len(idx_remove), seg_set["name"]
                    )
                )

                # remove the duplicate faces
                seg_set["set"] = np.delete(seg_set["set"], idx_remove, axis=0)

                all_tris = np.vstack((all_tris, seg_set1))

        ## TODO 2
        ## check if normals of the endocardium segments
        ## are pointing towards the centroid
        # volume = dsa.WrapDataObject( self._vtk_volume )
        # for cavity in self._cavities:
        #     for seg_name, seg_set in cavity.segment_sets.items():
        #         if "endocardium" in seg_name:
        #             seg_set

        return

    def _map_tags_to_remeshed_volume(
        self, source: vtk.vtkUnstructuredGrid, target: vtk.vtkUnstructuredGrid
    ):
        """Maps the tags of the myocardium to the remeshed volume.

        Notes
        -----
        This uses the original reference file as a source for interpolation. 
        """
        # _, _, cell_data_source, _ = get_tetra_info_from_unstructgrid(source)
        source_obj = dsa.WrapDataObject(source)
        target_obj = dsa.WrapDataObject(target)

        # get unique tags
        tags = source_obj.CellData["tags"]

        unique_tags = np.unique(tags)

        label_dict = copy.deepcopy(self.info.vtk_labels_to_use)
        tag_map = []

        # modify tags such that the valve/caps id point to the
        # attached myocardial part. E.g. Aortic valve tag -> Left ventricle tag
        # NOTE: should this come from CAVITY_DEFINITIONS? Can we extend
        # the vtk_labels_to_use dictionary?
        visited_parts = []
        # for cavity, parts in CAVITY_DEFINITIONS.items():
        for cavity in self._cavities:
            ref_part = cavity.labels[0]
            if "myocardium" not in ref_part:
                raise ValueError("Expecting myocardium substring in part")

            for part in cavity.labels[1:]:
                # skip if already visited
                if [s for s in visited_parts if part in s]:
                    # logger.debug("Already visited")
                    continue
                tag_map.append([label_dict[ref_part], label_dict[part]])
                visited_parts.append(part)

        # re-assigning tags
        tag_map = np.array(tag_map, dtype=int)
        for tag_pair in tag_map:
            tags[tags == tag_pair[1]] = tag_pair[0]

        # overwrite source array with new tags
        add_vtk_array(
            source, tags, name="tags", data_type="cell", array_type=int
        )

        # map discrete cell data
        target = vtk_map_discrete_cell_data(source, target, data_name="tags")

        return target

    def set_volume_mesh_vtk(self, filename):
        """Reads volume mesh (vtk) and sets the data"""
        if self._vtk_volume_temp:
            logger.warning("Overwriting existing data...")

        self._vtk_volume = vtk_read_mesh_file(filename)
        return

    def get_surface_from_volume(self):
        """Extracts the surface from the (remeshed) VTK object (volume)
        and assigns this mesh to each cavity for reference
        """

        self._vtk_surface = vtk_surface_filter(
            self._vtk_volume, keep_global_ids=True
        )

        # write_vtkdata_to_vtkfile(self._vtk_surface, "surface_of_remeshed_volume.vtk")

        for cavity in self._cavities:
            if "myocardium" in cavity.labels[0]:
                # surface_tag = cavity.vtk_ids[0]
                # surface, globalids = threshold_vtk_data( self._vtk_surface, surface_tag, surface_tag, "tags" )
                cavity._myocardium_surface = self._vtk_surface
            else:
                logger.error(
                    "Myocardium not in list: " + ", ".join(cavity.labels)
                )

        return

    def add_surface_to_cavities(self):
        """Adds the surface mesh to the cavity objects
        """
        for cavity in self._cavities:
            cavity.myocardium

    # private functions
    def _add_list_of_cavities(self):
        """ Initializes list of cavities given the model type and cavity definitions"""
        cavity_definitions = VALID_MODELS[self.info.model_type][
            "CavityDefinition"
        ]
        vtk_label_to_tag = self.info.vtk_labels_to_use

        cavities: List[Cavity] = []

        # instantiate cavity classes
        for part_names in cavity_definitions:
            # logger.debug(cavity)
            cavity_name = part_names[0].replace(" myocardium", "")

            vtk_ids = []
            for part_name in part_names:
                vtk_ids.append(vtk_label_to_tag[part_name])

            cavity_object = Cavity(
                name=cavity_name,
                vtk_labels=part_names,
                vtk_ids=vtk_ids,
                model_info=self.info,
            )

            cavity_object.id = vtk_ids[0]

            cavities.append(cavity_object)

        return cavities

    def _extract_relevant_parts(self, vtk_ugrid, vtk_labels_to_use):
        """Extracts the relevant parts from the input mesh file """

        # get labels to use
        logger.debug("Extracting tags...")
        labels_to_use = []
        for key, value in vtk_labels_to_use.items():
            labels_to_use.append(value)
            logger.debug("\tName: " + key + " | id: " + str(value))

        labels_to_use.sort()
        tags_to_use = labels_to_use
        # labels_to_use = np.array( labels_to_use )

        # extract relevant part of the mesh.
        vtk_thresholded = threshold_vtk_data_integers(
            vtk_ugrid, tags_to_use, data_name="tags"
        )

        return vtk_thresholded

    def _add_nodal_areas(self):
        """Computes and adds the nodal areas as point data to the
        vtk object
        """
        # self._vtk_volume

        return

    def _validate_cavities(self):
        """Validates list of defined cavities."""

        # NOTE: Any other validation steps that need to be done?
        # Now:
        # 1. check if cavity id is specified
        # 2. check if normal of cap points inward
        # 3. computes the volume (this is not a real validation step - so can be moved )

        # checks whether cavity ids are specified
        for cavity in self._cavities:
            if cavity.id == 0:
                raise ValueError("Cavity ID not set")

        # ensure that the normal of each cap is pointing towards the cavity center
        visited_caps = []
        for cavity in self._cavities:
            for cap in cavity.closing_caps:
                if cap.id in visited_caps:
                    # skip if already visited
                    continue

                d1 = np.linalg.norm(
                    cap.centroid + cap.normal - cavity.centroid
                )
                d2 = np.linalg.norm(
                    cap.centroid - cap.normal - cavity.centroid
                )
                if d2 > d1:
                    logger.debug(
                        "Flip cap normal such that it points inwards: %s"
                        % " ".join([cavity.name, cap.name])
                    )
                    cap.normal = cap.normal * -1
                visited_caps.append(cap.id)

        # For each cavity compute the volume of the enclosed surface
        # that is, the volume enclosed by the endocardium and cavity caps
        # NOTE: all normals should point inwards
        for cavity in self._cavities:
            # each cavity. Skip atrium cavities for now
            if "atrium" in cavity.name:
                logger.warning("Skipping volume computation %s" % cavity.name)
                continue

            cavity_triangles = np.empty((0, 3), dtype=int)
            for segsets in cavity.segment_sets:
                if "endocardium" in segsets["name"]:
                    cavity_triangles = np.vstack(
                        (cavity_triangles, segsets["set"])
                    )

            for cap in cavity.closing_caps:
                cavity_triangles = np.vstack(
                    (cavity_triangles, cap.closing_triangles)
                )

            # prepare stl file which is used to compute the volume
            # NOTE: could parse vtk.vtkPolyData of surface instead
            # but limited advantage over writing a (binary) stl file
            volume_mesh = dsa.WrapDataObject(self._vtk_volume)

            cavity_triangles = np.array(cavity_triangles, dtype=int)
            # write stl and compute its volume and store in cavity
            mesh = meshio.Mesh(
                points=volume_mesh.Points,
                cells=[("triangle", cavity_triangles)],
            )

            stl_path = os.path.join(
                self.info.working_directory,
                "closed_volume_{0}.stl".format(
                    "_".join(cavity.name.lower().split())
                ),
            )
            mesh.write(stl_path)

            cavity.compute_volume(stl_path)

            # os.remove( stl_path )

        return

    def _write_parts_to_stl(self):
        """Writes each unique tag to seperate .stl file (can be used in meshing step)
        """
        logger.debug("Writing parts to stl...[For Fluent use in future]")
        vtk_labels_to_use = self.info.vtk_labels_to_use

        nodes = self._nodes_surface
        faces = self._elements_surface
        tags = self._cell_data_surface["tags"]

        points = nodes
        for label, id in vtk_labels_to_use.items():
            cells = [("triangle", faces[tags == id, :])]
            # cell_data = { "tags": [tags] }

            mesh = meshio.Mesh(points, cells)

            path_to_file = os.path.join(
                self.info.working_directory, "part_{:0>3.0f}.stl".format(id)
            )

            mesh.write(path_to_file)

            add_solid_name_to_stl(path_to_file, label.replace(" ", "-"))

        return

    def extract_individual_parts_surface(self):
        """Extracts the surfaces of each individual part"""
        logger.debug(
            "Extracting individual parts as stl surfaces... [for Fluent?]"
        )
        for key, value in self.info.vtk_labels_to_use.items():
            vtk_part_volume, _ = threshold_vtk_data(
                self._vtk_volume_temp, value, value, "tags"
            )
            filename = "part_surface_{:0>3.0f}.stl".format(value)
            filepath = os.path.join(self.info.working_directory, filename)
            vtk_surface_to_stl(vtk_part_volume, filepath)

        return

    def _extract_septum(self):
        """Extract ids of elements that make up the septum"""
        # use septum segment set, convert to vtk polydata and
        # extrude this in (negative?) normal direction. This forms an
        # enclosed surface which can be used to tag the elements of the septum

        volume_obj = dsa.WrapDataObject(self._vtk_volume)
        points_volume = volume_obj.Points

        # select septum and convert this to vtk object
        for cavity in self._cavities:
            find = False
            for segset in cavity.segment_sets:
                if segset["name"] == "endocardium-septum":

                    # remove two layers of triangles connected to the boundary edge
                    endocardium_septum_segset = segset["set"]
                    endocardium_septum_segset_reduced = remove_triangle_layers_from_trimesh(
                        triangles=endocardium_septum_segset, iters=2
                    )
                    # convert to vtk object
                    endocardium_septum_vtk = create_vtk_surface_triangles(
                        points_volume, endocardium_septum_segset_reduced
                    )
                    # smooth surface
                    endocardium_septum_vtk = smooth_polydata(
                        endocardium_septum_vtk
                    )

                    # extrude in normal direction
                    endocardium_septum_vtk_extruded = extrude_polydata(
                        vtk_surface=endocardium_septum_vtk, extrude_by=-20
                    )

                    # get cell ids that are inside the extruded surface
                    cell_ids_septum = cell_ids_inside_enclosed_surface(
                        self._vtk_volume, endocardium_septum_vtk_extruded
                    )

                    cavity.element_sets.append(
                        {"name": "septum", "set": cell_ids_septum, "id": 1}
                    )

                    find = True
                    break
            if find:
                break

        return

    def _deprecated_export_mesh(self, directory: str):
        """Exports the mesh and all necessary data to load the 
        model again at a later time
        """
        # write mesh for simulation to file
        file_path = os.path.join(
            self.info.working_directory, "simulation_volume_mesh.vtk"
        )
        self.info.path_mesh = file_path
        write_vtkdata_to_vtkfile(self._vtk_volume, file_path)

        return

    def _OBSOLETE_separate_overlaps(self, volume: vtk.vtkUnstructuredGrid):
        """[OBSOLETE]Finds overlapping parts and separate these. These are then useful for 
        generating a consistent volume mesh

        Parameters
        ----------
        surface : vtk.vtkPolyData
            Surface representation as vtkPolyData object
        """
        nodes, tetra, cell_data, point_data = get_tetra_info_from_unstructgrid(
            volume
        )

        # extract surface
        surface = vtk_surface_filter(volume)
        (
            nodes_tris,
            tris,
            cell_data_tris,
            point_data_tris,
        ) = get_tri_info_from_polydata(surface)

        # find duplicate nodes:
        unique_nodes = np.unique(nodes.round(decimals=2), axis=0)

        tags = np.array(cell_data["tags"], dtype=int)
        # tags_nodes = np.zeros( nodes.shape[0], dtype = float )

        mesh = meshio.Mesh(
            points=nodes,
            cells=[("tetra", tetra)],
            cell_data={"tags": [np.array(tags, dtype=float)]},
        )

        mesh.write("volume_mesh_original.vtk")

        cavities_and_parts = copy.deepcopy(
            VALID_MODELS[self.info.model_type]["CavityDefinition"]
        )
        vtk_labels_to_use = self.info.vtk_labels_to_use

        # create new tags for rings
        used_tags = list(vtk_labels_to_use.values())
        if np.max(used_tags) < 100:
            new_tag = 100
        else:
            new_tag = np.max(used_tags) + 1

        for part_list in cavities_and_parts:
            part_name_myocardium = [
                part_name
                for part_name in part_list
                if "myocardium" in part_name
            ]
            if len(part_name_myocardium) != 1:
                raise ValueError(
                    "Expecting only one myocardium but found {0} ".format(
                        len(part_name_myocardium)
                    )
                )
            else:
                part_name_myocardium = part_name_myocardium[0]

            vtk_tag_myocardium = vtk_labels_to_use[part_name_myocardium]

            part_list.remove(part_name_myocardium)

            tetra_myocardium = tetra[tags == vtk_tag_myocardium, :]

            # mesh = meshio.Mesh( points = nodes,
            #                         cells = [("tetra", tetra_myocardium )] )
            # mesh.write("myocardium_{:0>3.0f}.vtk".format( vtk_tag_myocardium ) )

            # loop over remaining parts to find overlaps with reference part ( myocardium )
            # these are the "valve rings" - and should be added as such
            for part_name in part_list:
                part_id = vtk_labels_to_use[part_name]
                tetra_part = tetra[tags == part_id, :]
                node_indices_intersect = np.intersect1d(
                    np.unique(tetra_myocardium.ravel()),
                    np.unique(tetra_part.ravel()),
                )
                # find elements that refer these nodes and change tag
                # finds the overlapping faces
                mask = np.vstack(
                    (
                        np.isin(tetra_part[:, 0], node_indices_intersect),
                        np.isin(tetra_part[:, 1], node_indices_intersect),
                        np.isin(tetra_part[:, 2], node_indices_intersect),
                        np.isin(tetra_part[:, 3], node_indices_intersect),
                    )
                )

                tet_indices_to_tag = np.where(np.sum(mask, axis=0) == 3)[
                    0
                ]  # NOTE: THIS IS NOT CORRECT
                tetra_tagged = tetra_part[tet_indices_to_tag, :]

                global_indices_tags = np.where(tags == part_id)[0]

                tags[global_indices_tags[tet_indices_to_tag]] = new_tag
                # tags[ tags == part_id ][ tet_indices_to_tag ] = new_tag

                vtk_labels_to_use[part_name.replace("plane", "ring")] = new_tag

                # TODO: get triangles that make up shared boundary
                mask_tagged = mask[:, tet_indices_to_tag].transpose()
                tetra_tagged_1d = np.reshape(tetra_tagged, (tetra_tagged.size))
                masked_tagged_1d = np.reshape(mask_tagged, (tetra_tagged.size))
                triangles_tagged = np.reshape(
                    tetra_tagged_1d[np.where(masked_tagged_1d)],
                    (tetra_tagged.shape[0], 3),
                )

                mesh = meshio.Mesh(
                    points=nodes, cells=[("triangle", triangles_tagged)]
                )
                mesh.write(
                    "triangles_tagged_pid_{:0>3.0f}.stl".format(part_id)
                )

                new_tag = new_tag + 1

        # write file with specified tags
        # get list of ids to write
        tags_to_write = []
        for key, value in vtk_labels_to_use.items():
            if "plane" not in key:
                tags_to_write.append(value)
        tags_to_write = np.array(tags_to_write, dtype=int)
        tet_indices_to_write = np.isin(tags, tags_to_write)
        tet_indices_to_write = np.where(tet_indices_to_write)[0]

        tags = np.array(tags, dtype=float)

        tets_to_write = vtk_labels_to_use
        mesh = meshio.Mesh(
            points=nodes,
            cells=[("tetra", tetra[tet_indices_to_write, :])],
            cell_data={"tags": [tags[tet_indices_to_write]]},
        )

        mesh.write("volume_mesh_retagged.vtk")

        vtk_retagged = read_vtk_unstructuredgrid_file(
            "volume_mesh_retagged.vtk"
        )

        # extract surface
        surface = vtk_surface_filter(vtk_retagged)

        (
            nodes_tris,
            tris,
            cell_data_tris,
            point_data,
        ) = get_tri_info_from_polydata(surface)

        tags_tris = np.array(cell_data_tris["tags"], dtype=int)
        # export each tag as separate stl file
        unique_tags = np.unique(tags_tris)
        for unique_tag in unique_tags:
            tris_to_write = tris[tags_tris == unique_tag, :]
            mesh = meshio.Mesh(
                points=nodes_tris, cells=[("triangle", tris_to_write)]
            )
            mesh.write("triangles_pid_{:0>3.0f}.stl".format(unique_tag))

        vtk_surface_to_stl(surface, "surface_mesh.vtk")

        return

    def _OBSOLETE_separate_overlaps_surface(
        self, vtk_volume: vtk.vtkUnstructuredGrid
    ):
        """[OBSOLETE]
        """

        tags_to_use = list(self.info.vtk_labels_to_use.values())

        nodes = np.empty((0, 3), dtype=float)
        tris = np.empty((0, 3), dtype=int)
        tags = np.array([], dtype=int)

        # extract surface for each tag to use
        offset = 0
        for tag in tags_to_use:
            logger.debug("Extrating part id: {0}".format(tag))
            volume_thresholded, global_ids = threshold_vtk_data(
                vtk_volume, tag, tag, data_name="tags"
            )

            surface_threshold = vtk_surface_filter(volume_thresholded)

            (
                nodes_to_add,
                tris_to_add,
                cell_data,
                point_data_tris,
            ) = get_tri_info_from_polydata(surface_threshold)
            tags_to_add = np.array(cell_data["tags"], dtype=int)

            tris = np.vstack((tris, tris_to_add + offset))
            nodes = np.vstack((nodes, nodes_to_add))
            tags = np.append(tags, tags_to_add)
            offset = nodes.shape[0]

        mesh = meshio.Mesh(
            points=nodes,
            cells=[("triangle", tris)],
            cell_data={"tags": [np.array(tags, dtype=float)]},
        )

        mesh.write("surface_mesh_after_filter.vtk")

        # find duplicate nodes and remap triangles / tags accordingly
        nodes, tris = remove_duplicate_nodes(nodes, tris)

        mesh = meshio.Mesh(
            points=nodes,
            cells=[("triangle", tris)],
            cell_data={"tags": [np.array(tags, dtype=float)]},
        )

        mesh.write("surface_mesh_after_unique.vtk")

        # extract rings of valve parts

        cavities_and_parts = copy.deepcopy(
            VALID_MODELS[self.info.model_type]["CavityDefinition"]
        )
        vtk_labels_to_use = self.info.vtk_labels_to_use

        # create new tags for rings
        used_tags = list(vtk_labels_to_use.values())
        if np.max(used_tags) < 100:
            new_tag = 100
        else:
            new_tag = np.max(used_tags) + 1

        for part_list in cavities_and_parts:
            part_name_myocardium = [
                part_name
                for part_name in part_list
                if "myocardium" in part_name
            ]
            if len(part_name_myocardium) != 1:
                raise ValueError(
                    "Expecting only one myocardium but found {0} ".format(
                        len(part_name_myocardium)
                    )
                )
            else:
                part_name_myocardium = part_name_myocardium[0]

            vtk_tag_myocardium = vtk_labels_to_use[part_name_myocardium]

            part_list.remove(part_name_myocardium)

            tris_myocardium = tris[tags == vtk_tag_myocardium, :]

            mesh = meshio.Mesh(
                points=nodes, cells=[("triangle", tris_myocardium)]
            )
            mesh.write("myocardium_{:0>3.0f}.vtk".format(vtk_tag_myocardium))

            # loop over remaining parts to find overlaps with reference part ( myocardium )
            # these are the "valve rings" - and should be added as such
            for part_name in part_list:
                part_id = vtk_labels_to_use[part_name]
                tris_part = tris[tags == part_id, :]
                node_indices_intersect = np.intersect1d(
                    np.unique(tris_myocardium.ravel()),
                    np.unique(tris_part.ravel()),
                )

                np.savetxt(
                    "intersecting_nodes.csv",
                    nodes[node_indices_intersect, :],
                    delimiter=",",
                )
                # find elements that refer these nodes and change tag
                # finds the overlapping faces
                mask = np.vstack(
                    (
                        np.isin(tris_part[:, 0], node_indices_intersect),
                        np.isin(tris_part[:, 1], node_indices_intersect),
                        np.isin(tris_part[:, 2], node_indices_intersect),
                    )
                )

                tri_indices_to_tag = np.where(np.sum(mask, axis=0) == 3)[
                    0
                ]  # NOTE: THIS IS NOT CORRECT
                tris_tagged = tris_part[tri_indices_to_tag, :]

                # # remove the elements that use the intersecting nodes from the surface
                indices_duplicate = find_duplicate_elements(
                    tris_tagged, tris_myocardium
                )
                tris_myocardium[indices_duplicate, :] = -1

                global_indices_tags = np.where(tags == part_id)[0]

                tags[global_indices_tags[tri_indices_to_tag]] = new_tag
                # tags[ tags == part_id ][ tet_indices_to_tag ] = new_tag

                vtk_labels_to_use[part_name.replace("plane", "ring")] = new_tag

                # TODO: get triangles that make up shared boundary
                mask_tagged = mask[:, tri_indices_to_tag].transpose()
                tris_tagged_1d = np.reshape(tris_tagged, (tris_tagged.size))
                masked_tagged_1d = np.reshape(mask_tagged, (tris_tagged.size))
                triangles_tagged = np.reshape(
                    tris_tagged_1d[np.where(masked_tagged_1d)],
                    (tris_tagged.shape[0], 3),
                )

                # identify overlapping elements

                mesh = meshio.Mesh(
                    points=nodes, cells=[("triangle", triangles_tagged)]
                )
                stl_name = "triangles_pid_{:0>3.0f}.stl".format(part_id)
                mesh.write(stl_name)
                add_solid_name_to_stl(stl_name, "-".join(part_name.split()))

                new_tag = new_tag + 1

            # write myocardium without the overlapping elements
            tris_myo_to_write = tris_myocardium[:, 0] > -1
            mesh = meshio.Mesh(
                points=nodes,
                cells=[("triangle", tris_myocardium[tris_myo_to_write, :])],
            )
            stl_name = "triangles_pid_{:0>3.0f}.stl".format(vtk_tag_myocardium)
            mesh.write(stl_name)
            add_solid_name_to_stl(
                stl_name, "-".join(part_name_myocardium.split())
            )

        # find intersection between the myocardial parts

        return


if __name__ == "__main__":
    logger.info("Protected")
