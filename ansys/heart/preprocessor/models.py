"""Module containing classes for the various heart models which
can be created with the preprocessor
"""
# import json
import pathlib
import typing
import os
import numpy as np
import pickle, json

from ansys.heart.preprocessor.model_definitions import HEART_PARTS, LABELS_TO_ID
from ansys.heart.preprocessor.mesh.objects import Part, Mesh, SurfaceMesh, Cap, Cavity, Point
import ansys.heart.preprocessor.mesh.mesher as mesher
import ansys.heart.preprocessor.mesh.connectivity as connectivity
from ansys.heart.custom_logging import LOGGER
import ansys.heart.preprocessor.mesh.vtkmethods as vtkmethods


class ModelInfo:
    """Contains model information"""

    @property
    def database(self) -> str:
        return self._database

    @database.setter
    def database(self, value: str):
        valid_databases = ["Strocchi2020", "Cristobal2021"]
        if value not in valid_databases:
            raise ValueError(
                "{0} not a valid database name. Please specify one of the following database names: {1}".format(
                    value, valid_databases
                )
            )
        self._database = value

    def __init__(
        self,
        database: str,
        work_directory: pathlib.Path,
        path_to_case: pathlib.Path,
        path_to_simulation_mesh: pathlib.Path = None,
    ) -> None:
        self.database = database
        """Name of the database to use"""
        self.workdir = work_directory
        """Path to the working directory"""
        self.path_to_original_mesh = path_to_case
        """Path to the original mesh file"""
        self.path_to_simulation_mesh = path_to_simulation_mesh
        """Path to simulation mesh (in vtk format)"""
        self.labels_to_ids = LABELS_TO_ID[database]
        """Dict that maps labels > part/tag id"""
        self.ids_to_labels = dict((v, k) for k, v in LABELS_TO_ID[database].items())
        """Inverted dict that maps part/tag id > labels"""

        pass

    def clean_workdir(
        self,
        extensions_to_remove: typing.List[str] = [".stl", ".vtk", ".msh.h5"],
        remove_all: bool = False,
    ):
        """Removes files with extension present in the working directory

        Parameters
        ----------
        extensions_to_remove : typing.List[str], optional
            List of extensions to remove, by default [".stl", ".vtk", ".msh.h5"]
        remove_all: bool, optional
            Flag indicating whether to remove files with any extension.
            Keeps files/folder without extension
        """
        import glob as glob

        if not remove_all:
            for ext in extensions_to_remove:
                files = glob.glob(os.path.join(self.workdir, "*" + ext))
                for file in files:
                    os.remove(file)
        elif remove_all:
            files = glob.glob(os.path.join(self.workdir, "*.*"))
            for file in files:
                os.remove(file)

        return

    def create_workdir(self):
        """Creates the working directory if it doesn't exist"""
        if not os.path.isdir(self.workdir):
            os.makedirs(self.workdir)
        return

    def dump_info(self, filename: pathlib.Path = None):
        """Dumps model information to file"""
        if not filename:
            filename = os.path.join(self.workdir, "model_info.json")
        with open(filename, "w") as file:
            file.write(json.dumps(self.__dict__, indent=4))
        return


class HeartModel:
    """Parent class for heart models"""

    @property
    def parts(self) -> typing.List[Part]:
        """Returns list of parts"""
        parts = []
        for key, value in self.__dict__.items():
            attribute = getattr(self, key)
            if type(attribute) == Part:
                parts.append(attribute)
        return parts

    @property
    def part_names(self) -> typing.List[str]:
        """Returns list of part names"""
        part_names = []
        for part in self.parts:
            part_names.append(part.name)
        return part_names

    def __init__(self, info: ModelInfo) -> None:
        self.info = info
        """Model meta information"""
        self.mesh = Mesh()
        """Modified mesh"""
        self.mesh_raw = Mesh()
        """Raw input mesh"""

        self._add_subparts()
        """Adds any subparts"""
        self._add_labels_to_parts()
        """Adds appropiate vtk labels to the parts"""
        pass

    def extract_simulation_mesh(self, clean_up: bool = False):
        """Updates the model"""
        self.read_input_mesh()
        self._remove_unused_tags()
        self._prepare_for_meshing()
        self._remesh()
        self._update_parts()

        if clean_up:
            self.info.clean_workdir(["part*.stl", "cavity*.stl"])

        return

    def get_part(self, name: str, by_substring: bool = False):
        """Gets specific part based on part name"""
        found = False
        for part in self.parts:
            if part.name == name:
                return part
            if by_substring:
                if name in part.name:
                    return part
        if not found:
            return None

    def add_part(self, part_name: str):
        """Dynamically adds a part as an attribute to the object"""
        setattr(self, "_".join(part_name.lower().split()), Part(name=part_name))
        return

    def print_info(self):
        """Prints information about the model"""
        LOGGER.info("*****************************************")
        LOGGER.info("*****************************************")
        LOGGER.info("Mesh info:")
        LOGGER.info("Number of tetra: {:d}".format(self.mesh.tetrahedrons.shape[0]))
        LOGGER.info("Number of nodes: {:d}".format(self.mesh.nodes.shape[0]))
        LOGGER.info("-----------------------------------------")
        for ii, part in enumerate(self.parts):
            LOGGER.info("{:d}. part name: {:}".format(ii + 1, part.name))
            LOGGER.info("\tnumber of tetrahedrons: {:d}\n".format(len(part.element_ids)))

            for surface in part.surfaces:
                LOGGER.info(
                    "\tsurface: {:} | # faces: {:d}".format(surface.name, surface.faces.shape[0])
                )
            for cap in part.caps:
                LOGGER.info("\tcap: {:} | # nodes {:d}".format(cap.name, len(cap.node_ids)))
            if part.cavity:
                LOGGER.info(
                    "\tcavity: {:} | volume: {:.1f} [mm3]".format(
                        part.cavity.name, part.cavity.volume
                    )
                )
            LOGGER.info("-----------------------------------------")
        LOGGER.info("*****************************************")
        LOGGER.info("*****************************************")
        return

    def read_input_mesh(self):
        """Reads the input mesh defined in ModelInfo"""
        if not os.path.isfile(self.info.path_to_original_mesh):
            raise ValueError("Please specify a valid path to the input file")

        if self.info.database == "Cristobal2021":
            LOGGER.warning("Changing data fields of vtk file")
            self.mesh_raw.read_mesh_file_cristobal2021(self.info.path_to_original_mesh)
        else:
            self.mesh_raw.read_mesh_file(self.info.path_to_original_mesh)

        try:
            self.mesh_raw.cell_data["tags"]
        except:
            raise KeyError("Expecting a field 'tags' in mesh_raw.cell_data")
        return

    def dump_model(self, filename: pathlib.Path = None):
        """Saves model to file"""
        LOGGER.debug("Writing model to disk")
        if not filename:
            filename = os.path.join(self.info.workdir, "heart_model.pickle")
        # cleanup model object for more efficient storage
        # NOTE deleting faces, nodes of surfaces does not affect size
        del self.mesh_raw
        with open(filename, "wb") as file:
            pickle.dump(self, file)
        self.info.dump_info()
        return

    def load_model(self, filename: pathlib.Path):
        """Loads the model from file"""
        with open(filename, "rb") as file:
            self = pickle.load(file)
        return

    def _add_labels_to_parts(self):
        """Uses model definitions to add corresponding vtk labels to the part"""
        for part in self.parts:
            if part.name == "Septum":
                continue
            part.tag_labels = HEART_PARTS[part.name]["VTKLabels"]
            # add tag ids
            part.tag_ids = []
            for tag_label in part.tag_labels:
                part.tag_ids.append(LABELS_TO_ID[self.info.database][tag_label])
        return

    def _add_subparts(self):
        """Adds subparts to parts of type ventricle"""
        for part in self.parts:
            if part.part_type in ["ventricle"]:
                part._add_myocardium_part()
                if "Left ventricle" in part.name:
                    part._add_septum_part()

        return

    def _remove_unused_tags(self):
        """Extract only the tags of interest"""
        # collect all used tags
        tag_ids = []
        for part in self.parts:
            if not part.tag_ids:
                continue
            tag_ids.extend(part.tag_ids)

        self.mesh_raw.keep_elements_with_value(tag_ids, "tags")

        return

    def _get_used_element_ids(self):
        """Returns array of used element ids"""
        element_ids = np.empty(0, dtype=int)
        for part in self.parts:
            element_ids = np.append(element_ids, part.element_ids)

        return element_ids

    def _get_endo_epicardial_surfaces(self):
        """get endo and epicardial surfaces

        Note
        ----
        Also obtains the septum in case of a BiVentricle, FourChamber or FullHeart model

        """
        surfaces_to_add = []
        for boundary in self.mesh_raw.boundaries:
            if not "ventricle" in boundary.name and not "atrium" in boundary.name:
                continue

            if "left-ventricle-myocardium" in boundary.name and isinstance(
                self, (BiVentricle, FourChamber, FullHeart)
            ):
                boundary_name_suffix = ["epicardium", "endocardium", "septum"]
            else:
                boundary_name_suffix = ["epicardium", "endocardium"]
            LOGGER.debug("Extracting : {} from {}".format(boundary_name_suffix, boundary.name))
            region_ids = boundary.separate_connected_regions()

            unique_regions, counts = np.unique(region_ids, return_counts=True)

            assert len(unique_regions) == len(
                boundary_name_suffix
            ), "Number of extracted regions does not match number of expected regions"

            # NOTE: number of cells do not a guarantee to distinguish between endo and epicardium
            # sort by bounding box volume instead.
            volumes = []
            surfaces: typing.List[SurfaceMesh] = []
            for ii, region_id in enumerate(unique_regions):
                mask = region_ids == region_id
                surface = SurfaceMesh(faces=boundary.faces[mask, :], nodes=self.mesh_raw.nodes)
                volumes.append(surface.compute_bounding_box()[1])
                surfaces.append(surface)
            # sort by volume of bounding box
            order = np.flip(np.argsort(volumes))
            surfaces = [surfaces[idx] for idx in order]
            # update names:
            for ii, surface in enumerate(surfaces):
                surface.name = boundary.name.replace("myocardium", boundary_name_suffix[ii])

            surfaces_to_add += surfaces

        self.mesh_raw.boundaries = self.mesh_raw.boundaries + surfaces_to_add

        # rename septum boundary from left to right ventricle
        if isinstance(self, (BiVentricle, FourChamber, FullHeart)):
            septum_boundary = [
                b for b in self.mesh_raw.boundaries if b.name == "left-ventricle-septum"
            ]
            if septum_boundary and len(septum_boundary) == 1:
                septum_boundary[0].name = "right-ventricle-septum"
            else:
                raise ValueError("Expecting a single boundary named 'left-ventricle-septum'")

        return

    def _prepare_for_meshing(self):
        """Prepares the input for volumetric meshing with Fluent meshing

        Note
        ----
        Extracts boundary surfaces of the part and the interface surfaces between parts.
        These surfaces are written in .stl format which can be imported in Fluent Meshing

        """

        # establish connectivity of the raw mesh
        self.mesh_raw.establish_connectivity()
        # mark interface pairs - and get all unique interface-pairs
        mask, interface_pairs = self.mesh_raw.get_mask_interface_faces(return_pairs=True)

        # loop over each interface pair and write an stl file for each relevant interface
        # skip any interfaces that do not include valve/inlets
        # NOTE: Collect this into a separate method.
        interface_pairs1 = []
        interface_names = []
        for pair in interface_pairs:
            labels = [self.info.ids_to_labels[pair[0]], self.info.ids_to_labels[pair[1]]]
            # skip any interfaces that do not involve a inlet or valve
            if not (
                np.any(["inlet" in e for e in labels]) or np.any(["valve" in e for e in labels])
            ):
                LOGGER.debug("Skipping interface pair: {0} | {1}".format(labels[0], labels[1]))
                continue
            interface_pairs1.append(pair)
            interface_name = "-".join("_".join(labels).lower().split())
            interface_names.append(interface_name)
        interface_pairs = interface_pairs1

        # add interfaces and perform smoothing
        self.mesh_raw.add_interfaces(interface_pairs, interface_names)
        self.mesh_raw.smooth_interfaces()

        # extract boundaries of relevant parts (ignore valves and inlets)
        for part in self.parts:
            if not part.tag_labels:
                continue
            labels_to_use = [l for l in part.tag_labels if "inlet" not in l and "valve" not in l]
            part_ids_use = [self.info.labels_to_ids[l] for l in labels_to_use]
            part_names = ["-".join(l.split()).lower() for l in labels_to_use]
            self.mesh_raw.add_boundaries(part_ids_use, part_names)

        self._get_endo_epicardial_surfaces()

        # write interfaces and boundaries to .stl file (input for mesher)
        surfaces_to_skip = [
            "left-ventricle-myocardium",
            "right-ventricle-myocardium",
            "left-atrium-myocardium",
            "right-atrium-myocardium",
        ]
        # boundaries_to_write = [
        #     b for b in self.mesh_raw.boundaries if "valve" not in b.name and "inlet" not in b.name
        # ]
        for surface in self.mesh_raw.interfaces + self.mesh_raw.boundaries:
            if surface.name in surfaces_to_skip:
                continue
            surface.write_to_stl(os.path.join(self.info.workdir, "part_" + surface.name + ".stl"))

        return

    def _remesh(self):
        """Uses the generated files to remesh the surfaces and volume"""
        LOGGER.info("Remeshing volume...")
        path_mesh_file = os.path.join(self.info.workdir, "fluent_volume_mesh.msh.h5")
        mesher.fluentmeshing(
            self.info.workdir,
            path_mesh_file,
            mesh_size=1.5,
            journal_type="improved",
            show_gui=True,
        )
        path_mesh_file_vtk = path_mesh_file.replace(".msh.h5", ".vtk")
        tetra, face_zones, nodes = mesher.hdf5.fluenthdf5_to_vtk(path_mesh_file, path_mesh_file_vtk)

        # update mesh object
        self.mesh.tetrahedrons = tetra
        self.mesh.nodes = nodes
        for face_zone_name, face_zone in face_zones.items():
            self.mesh.boundaries.append(
                SurfaceMesh(
                    name=face_zone_name,
                    faces=face_zone["faces"][
                        :, [0, 2, 1]
                    ],  # ensures normals pointing away from the volume mesh
                    nodes=self.mesh.nodes,
                    sid=face_zone["zone-id"],
                )
            )
        self._map_data_to_remeshed_volume()

        return

    def _map_data_to_remeshed_volume(self):
        """Maps the data from the original (volume) mesh to the remeshed (volume) mesh
        including part-ids
        """

        # get list of tag ids to keep for mapping
        mapper = self.info.labels_to_ids
        labels = [l for p in self.parts if p.tag_labels for l in p.tag_labels]
        tag_ids = [mapper[l] for l in labels if "valve" not in l and "inlet" not in l]

        self.mesh_raw.keep_elements_with_value(tag_ids, "tags")

        # get list of all arrays in original mesh
        array_names = list(self.mesh_raw.cell_data.keys()) + list(self.mesh_raw.point_data.keys())

        # write (original) raw mesh and and new mesh to disk
        filename_original = os.path.join(self.info.workdir, "mesh_raw.vtk")
        filename_remeshed = os.path.join(self.info.workdir, "mesh.vtk")

        self.mesh_raw.write_to_vtk(filename_original)
        self.mesh.write_to_vtk(filename_remeshed)

        source = vtkmethods.read_vtk_unstructuredgrid_file(filename_original)
        target = vtkmethods.read_vtk_unstructuredgrid_file(filename_remeshed)

        # map uvc arrays
        uvc_array_names = [k for k in self.mesh_raw.point_data.keys() if "uvc" in k]
        target = vtkmethods.vtk_map_continuous_data(
            source=source,
            target=target,
            array_names_to_include=uvc_array_names,
        )
        # map tags
        target = vtkmethods.vtk_map_discrete_cell_data(source, target, data_name="tags")

        # interpolate remaining fields
        remaining_arrays = set(array_names) - set(uvc_array_names + ["tags"])

        target = vtkmethods.vtk_map_continuous_data(
            source=source,
            target=target,
            array_names_to_include=remaining_arrays,
        )

        # assign cell and point data to mesh object
        (
            _,
            _,
            self.mesh.cell_data,
            self.mesh.point_data,
        ) = vtkmethods.get_tetra_info_from_unstructgrid(target)
        # self.mesh.part_ids = self.mesh.cell_data["tags"].astype(int)

        # For any non-ventricular points assign -100 to uvc coordinates
        mapper = self.info.ids_to_labels
        ventricular_tags = [tid for tid in tag_ids if "ventricle" in mapper[tid]]
        mask = np.isin(self.mesh.part_ids, ventricular_tags, invert=True)
        node_ids_to_modify = np.unique(self.mesh.tetrahedrons[mask, :])
        for key, value in self.mesh.point_data.items():
            if "uvc_" in key:
                self.mesh.point_data[key][node_ids_to_modify] = -100

        # mesh with interpolated data is the simulation mesh
        path_to_simulation_mesh = os.path.join(self.info.workdir, "simulation_mesh.vtk")
        self.mesh.write_to_vtk(path_to_simulation_mesh)
        self.info.path_to_simulation_mesh = path_to_simulation_mesh

        # cleanup
        os.remove(filename_original)
        os.remove(filename_remeshed)
        return

    def _update_parts(self):
        """Updates the parts using the (re)meshed volume

        Notes
        -----
        1. Extracts septum
        2. Updates Parts to include element ids of the respective part
        3. Assign surfaces to each part
        4. Extracts the closing caps
        5. Creates cavities

        """

        self._extract_septum()
        self._assign_elements_to_parts()
        self._assign_surfaces_to_parts()
        self._assign_caps_to_parts()
        self._assign_cavities_to_parts()
        self._extract_apex()

        return

    def _extract_septum(self):
        """Separates the septum elements from the left ventricle

        Note
        ----
        Uses the septum surface of the right ventricle
        """
        if not isinstance(self, (BiVentricle, FourChamber, FullHeart)):
            LOGGER.warning("Model type: {0} Not extracting septum elements".format(type(self)))
            return None

        surface_septum = [s for s in self.mesh.boundaries if "septum" in s.name]
        if len(surface_septum) != 1:
            raise ValueError("Expecting only one surface that contains string: 'septum'")
        surface_septum = surface_septum[0]

        # extrude septum surface
        faces_septum = connectivity.remove_triangle_layers_from_trimesh(
            surface_septum.faces, iters=1
        )

        septum_surface_vtk = vtkmethods.create_vtk_surface_triangles(self.mesh.nodes, faces_septum)

        septum_surface_vtk = vtkmethods.smooth_polydata(septum_surface_vtk)

        septum_surface_vtk_extruded = vtkmethods.extrude_polydata(septum_surface_vtk, -20)

        filename_vtk = os.path.join(self.info.workdir, "volume_mesh.vtk")
        self.mesh.write_to_vtk(filename_vtk)
        volume_vtk = vtkmethods.read_vtk_unstructuredgrid_file(filename_vtk)

        element_ids_septum = vtkmethods.cell_ids_inside_enclosed_surface(
            volume_vtk, septum_surface_vtk_extruded
        )

        # assign to septum
        part = next(part for part in self.parts if part.name == "Septum")
        part.element_ids = element_ids_septum
        # remove these element ids from the left-ventricle
        part = next(part for part in self.parts if part.name == "Left ventricle")
        mask = np.isin(part.element_ids, element_ids_septum, invert=True)
        part.element_ids = part.element_ids[mask]

        os.remove(filename_vtk)

        return

    def _extract_apex(self):
        """Extracts the apex for both the endocardium and epicardium of each ventricle

        Note
        ----
        Apex defined as the point furthest from the mid-point between caps/valves

        """
        ventricles = [p for p in self.parts if "ventricle" in p.name]
        surface_substrings = ["endocardium", "epicardium"]
        for ventricle in ventricles:
            # get reference point (center point between two caps)
            cap_centroids = [c.centroid for c in ventricle.caps]
            ref_point = np.mean(np.array(cap_centroids), axis=0)
            for surface_substring in surface_substrings:
                surface = next(s for s in ventricle.surfaces if surface_substring in s.name)
                apical_node_id = surface.node_ids[
                    np.argmax(
                        np.linalg.norm(surface.nodes[surface.node_ids, :] - ref_point, axis=1)
                    )
                ]
                ventricle.apex_points.append(
                    Point(
                        name="apex " + surface_substring,
                        node_id=apical_node_id,
                        xyz=surface.nodes[apical_node_id, :],
                    )
                )

        return

    def _assign_elements_to_parts(self):
        """Gets the element ids of each part and assign these to the Part objects"""
        # get element ids of each part
        used_element_ids = self._get_used_element_ids()
        for part in self.parts:
            if len(part.element_ids) > 0:
                LOGGER.warning(
                    "Part {0} seems to already have elements assigned: skipping".format(part.name)
                )
                continue

            element_ids = np.where(np.isin(self.mesh.part_ids, part.tag_ids))[0]
            element_ids = element_ids[np.isin(element_ids, used_element_ids, invert=True)]
            part.element_ids = element_ids

        summ = 0
        for part in self.parts:
            LOGGER.debug("Num elements in {0}: {1}".format(part.name, part.element_ids.shape[0]))
            summ = summ + part.element_ids.shape[0]
        LOGGER.debug("Total num elements: {}".format(summ))

        if summ != self.mesh.tetrahedrons.shape[0]:
            raise ValueError(
                "{0} elements assigned to parts - but {1} exist in mesh".format(
                    summ, self.mesh.tetrahedrons.shape[0]
                )
            )

    def _assign_surfaces_to_parts(self):
        """Assigns surfaces generated during remeshing to model parts"""

        for part in self.parts:
            for surface in part.surfaces:
                boundary_name = "-".join(surface.name.lower().split())
                boundary_surface = self.mesh.get_surface_from_name(boundary_name)
                if boundary_surface:
                    surface.faces = boundary_surface.faces
                    surface.nodes = boundary_surface.nodes
                else:
                    LOGGER.warning("Could not find matching surface for: {0}".format(surface.name))

        return

    def _assign_caps_to_parts(self):
        """Uses connectivity to obtain cap boundaries and adds these to their respective parts"""

        used_boundary_surface_names = [s.name for p in self.parts for s in p.surfaces]
        remaining_surfaces = list(set(self.mesh.boundary_names) - set(used_boundary_surface_names))
        remaining_surfaces1: typing.List[SurfaceMesh] = []
        for surface in self.mesh.boundaries:
            if surface.name not in remaining_surfaces:
                continue
            surface.get_boundary_edges()
            remaining_surfaces1.append(surface)

        # find intersection between remaining surfaces and part surfaces
        # This will find the valve/cap nodes
        for part in self.parts:
            for surface in part.surfaces:
                surface.get_boundary_edges()
                if not "endocardium" in surface.name:
                    continue
                for edge_group in surface.edge_groups:
                    if edge_group.type != "closed":
                        raise ValueError("Expecting closed group of edges")

                    for surf in remaining_surfaces1:
                        if "valve" not in surf.name and "inlet" not in surf.name:
                            continue
                        if "myocardium" not in surf.name:
                            continue

                        if np.any(np.isin(edge_group.edges, surf.boundary_edges)):
                            name_valve = next(
                                n for n in surf.name.split("_") if "valve" in n or "inlet" in n
                            )
                            name_valve = name_valve.replace("-plane", "").replace("-inlet", "")

                            cap = Cap(name=name_valve, node_ids=edge_group.edges[:, 0])

                            # get approximate cavity centroid to check normal of cap
                            cavity_centroid = surface.compute_centroid()

                            cap.tesselate()
                            p1 = (
                                surf.nodes[
                                    cap.triangles[:, 1],
                                ]
                                - surf.nodes[
                                    cap.triangles[:, 0],
                                ]
                            )
                            p2 = (
                                surf.nodes[
                                    cap.triangles[:, 2],
                                ]
                                - surf.nodes[
                                    cap.triangles[:, 0],
                                ]
                            )
                            normals = np.cross(p1, p2)
                            cap_normal = np.mean(normals, axis=0)
                            cap_normal = cap_normal / np.linalg.norm(cap_normal)
                            cap_centroid = np.mean(surf.nodes[cap.node_ids, :], axis=0)
                            d1 = np.linalg.norm(cap_centroid + cap_normal - cavity_centroid)
                            d2 = np.linalg.norm(cap_centroid - cap_normal - cavity_centroid)
                            if d1 > d2:
                                LOGGER.debug(
                                    "Flipping order of nodes on cap to ensure normal pointing inward"
                                )
                                cap.node_ids = np.flip(cap.node_ids)
                                cap.tesselate()

                            cap.centroid = np.mean(surf.nodes[cap.node_ids, :], axis=0)

                            part.caps.append(cap)
                            LOGGER.debug("Cap: {0} closes {1}".format(name_valve, surface.name))
                            break

        # replace caps of atria by caps of ventricle
        for part in self.parts:
            if not "atrium" in part.name:
                continue
            for cap in part.caps:
                # replace with cap in ventricle
                cap_ref = [
                    c
                    for p in self.parts
                    if "ventricle" in p.name
                    for c in p.caps
                    if c.name == cap.name
                ]
                if len(cap_ref) == 1:
                    LOGGER.debug(
                        "Replacing cap {0} of part{1}: with that of the ventricle".format(
                            cap.name, part.name
                        )
                    )
                    # note: flip order to make sure normal is pointing inwards
                    cap.node_ids = np.flip(cap_ref[0].node_ids)
                    cap.tesselate()

        # As a consequence we need to add interface region to endocardium of atria or ventricle
        # current approach is to add these to the atria
        for part in self.parts:
            if "Left atrium" in part.name:
                interface_name = "mitral-valve-plane"
            elif "Right atrium" in part.name:
                interface_name = "tricuspid-valve-plane"
            else:
                continue
            interfaces = [s for s in remaining_surfaces1 if interface_name in s.name]
            endocardium = next(s for s in part.surfaces if "endocardium" in s.name)
            # append interface faces to endocardium
            for interface in interfaces:
                endocardium.faces = np.vstack([endocardium.faces, interface.faces])

        return

    def _assign_cavities_to_parts(self):
        """Creates cavities based on endocardium surfaces and cap definitions
        And assigns these to the parts of the model"""

        # rename septum to right ventricle endocardium septum
        part = self.get_part("Right ventricle", True)
        for surface in part.surfaces:
            if "Right ventricle septum" in surface.name:
                surface.name = surface.name.replace("septum", "endocardium septum")

        # construct cavities with endocardium and caps
        for part in self.parts:
            if "atrium" not in part.name and "ventricle" not in part.name:
                continue
            cavity_faces = np.empty((0, 3), dtype=int)

            surfaces = [s for s in part.surfaces if "endocardium" in s.name]
            for surface in surfaces:
                cavity_faces = np.vstack([cavity_faces, surface.faces])

            for cap in part.caps:
                cavity_faces = np.vstack([cavity_faces, cap.triangles])

            surface = SurfaceMesh(
                name=part.name + " cavity", faces=cavity_faces, nodes=self.mesh.nodes
            )
            part.cavity = Cavity(surface=surface, name=part.name)
            part.cavity.compute_centroid()

            volume = part.cavity.compute_volume()
            LOGGER.debug("Volume of cavity: {0} = {1}".format(part.cavity.name, volume))

            part.cavity.surface.write_to_stl(
                os.path.join(self.info.workdir, "-".join(part.cavity.surface.name.lower().split()))
            )

        return


class LeftVentricle(HeartModel):
    """Model of just the left ventricle"""

    def __init__(self, info: ModelInfo) -> None:

        self.left_ventricle = Part(name="Left ventricle", part_type="ventricle")
        """Left ventricle part"""
        # remove septum - not used in left ventricle only model
        del self.left_ventricle.septum

        super().__init__(info)
        pass


class BiVentricle(HeartModel):
    """Model of the left and right ventricle"""

    def __init__(self, info: ModelInfo) -> None:

        self.left_ventricle = Part(name="Left ventricle", part_type="ventricle")
        """Left ventricle part"""
        self.right_ventricle = Part(name="Right ventricle", part_type="ventricle")
        """Right ventricle part"""
        self.septum = Part(name="Septum", part_type="septum")
        """Septum"""

        super().__init__(info)
        pass


class FourChamber(HeartModel):
    """Model of the left/right ventricle and left/right atrium"""

    def __init__(self, info: ModelInfo) -> None:

        self.left_ventricle = Part(name="Left ventricle", part_type="ventricle")
        """Left ventricle part"""
        self.right_ventricle = Part(name="Right ventricle", part_type="ventricle")
        """Right ventricle part"""
        self.septum = Part(name="Septum", part_type="septum")
        """Septum"""

        self.left_atrium = Part(name="Left atrium", part_type="atrium")
        """Left atrium part"""
        self.right_atrium = Part(name="Right atrium", part_type="atrium")
        """Right atrium part"""

        super().__init__(info)

        pass


class FullHeart(HeartModel):
    """Model of the left/right ventricle,  left/right atrium, aorta
    and pulmonary artery
    """

    def __init__(self, info: ModelInfo) -> None:

        self.left_ventricle = Part(name="Left ventricle", part_type="ventricle")
        """Left ventricle part"""
        self.right_ventricle = Part(name="Right ventricle", part_type="ventricle")
        """Right ventricle part"""
        self.septum = Part(name="Septum", part_type="septum")
        """Septum"""
        self.left_atrium = Part(name="Left atrium", part_type="atrium")
        """Left atrium part"""
        self.right_atrium = Part(name="Right atrium", part_type="atrium")
        """Right atrium part"""
        self.aorta = Part(name="Aorta", part_type="artery")
        """Aorta part"""
        self.pulmonary_artery = Part(name="Pulmonary artery", part_type="artery")
        """Pulmonary artery part"""

        super().__init__(info)

        pass


if __name__ == "__main__":
    info = ModelInfo(database="Strocchi2020", work_directory="tmp", path_to_case="test.case")

    model = LeftVentricle(info)
    print("LeftVentricle:")
    model.print_parts()

    model = BiVentricle(info)
    print("BiVentricle:")
    model.print_parts()

    model = FourChamber(info)
    print("FourChamber:")
    model.print_parts()

    model = FullHeart(info)
    print("FullHeart:")
    model.print_parts()

    print(model.part_names)
