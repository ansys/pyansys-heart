import os

from ansys.heart.preprocessor.heart_mesh import HeartMesh
from ansys.heart.preprocessor.model_information import ModelInformation

# import logger
from ansys.heart.custom_logging import LOGGER
from ansys.heart.preprocessor.vtk_module import (
    read_vtk_polydata_file,
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

        write_folder = self.info.working_directory

        self._mesh.load_raw_mesh()
        vtk_labels = copy.deepcopy(self.info.vtk_labels_to_use)
        # invert label map
        vtk_labels_inverse = dict((v, k) for k, v in vtk_labels.items())
        tags_to_use = list(vtk_labels.values())

        # model type determines what parts to extract
        from ansys.heart.preprocessor.vtk_module import (
            get_tetra_info_from_unstructgrid,
            write_vtkdata_to_vtkfile,
            vtk_surface_to_stl,
            create_vtk_surface_triangles,
            get_connected_regions,
        )
        from ansys.heart.preprocessor.mesh_module import add_solid_name_to_stl
        from ansys.heart.preprocessor.extract_regions_manual import (
            face_tetra_connectivity,
            get_face_type,
        )
        import numpy as np

        points, tetra, cell_data, point_data = get_tetra_info_from_unstructgrid(
            self._mesh._vtk_volume_raw
        )

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

        # # NOTE: delete any tag pairs?
        # tag_pairs[np.any(tag_pairs == -1, axis=1), :]

        # # assign -1 to any cell which is not used
        # part_ids[np.isin(part_ids, tags_to_use, invert=True)] = -1

        # remove myocardium<>myocardium interfaces
        # remove myocardium<>aorta wall and myocardiun<>pulmonary artery wall interfaces
        tag_pairs_1 = np.empty((0, 2), dtype=int)
        for tag_pair in tag_pairs:
            labels = [vtk_labels_inverse[tag_pair[0]], vtk_labels_inverse[tag_pair[1]]]
            res1 = ["myocardium" in e for e in labels]
            res2 = ["wall" in e for e in labels]
            if np.all(["myocardium" in e for e in labels]):
                # skip if myocardium-myocardium intersection
                print("Skipping {0} | {1}".format(*tag_pair))
                continue
            if (
                np.sum(["myocardium" in e for e in labels]) == 1
                and np.sum(["wall" in e for e in labels]) == 1
            ):
                # skip if myocardium<>pulmonary artery or myocardium<>aorta wall interfaces
                print("Skipping {0} | {1}".format(*tag_pair))
                continue
            tag_pairs_1 = np.vstack([tag_pairs_1, tag_pair])
        tag_pairs = tag_pairs_1

        # for labels_pair in label_pairs:
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
            vtk_surface = create_vtk_surface_triangles(points, faces_interface)
            filename = os.path.join(
                write_folder, "part_interface_{:02d}_{:02d}.stl".format(tags_pair[0], tags_pair[1])
            )
            vtk_surface_to_stl(vtk_surface, filename, solid_name=pair_name)

        # exclude valves/inlets from extraction
        part_ids_to_extract = []
        for label, pid in vtk_labels.items():
            if "inlet" in label or "valve" in label:
                continue
            part_ids_to_extract.append(pid)
        part_ids_to_extract = np.sort(part_ids_to_extract)

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
