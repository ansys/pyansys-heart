"""Module containing classes for the various heart models."""
import json

# from ansys.heart import LOG as LOGGER
import logging
import os

# import json
import pathlib
import pickle
from typing import List, Union

LOGGER = logging.getLogger("pyheart_global.preprocessor")
import ansys.heart.preprocessor.mesh.connectivity as connectivity
import ansys.heart.preprocessor.mesh.mesher as mesher
from ansys.heart.preprocessor.mesh.objects import (
    BeamMesh,
    Cap,
    Cavity,
    Mesh,
    Part,
    Point,
    SurfaceMesh,
)
import ansys.heart.preprocessor.mesh.vtkmethods as vtkmethods
from ansys.heart.preprocessor.model_definitions import HEART_PARTS, LABELS_TO_ID
import meshio
import numpy as np
import pyvista as pv
from scipy.spatial.transform import Rotation as R


class ModelInfo:
    """Contains model information."""

    @property
    def database(self) -> str:
        """Name of the database to use."""
        return self._database

    @database.setter
    def database(self, value: str):
        valid_databases = ["Strocchi2020", "Cristobal2021", "LabeledSurface"]
        if value not in valid_databases:
            raise ValueError(
                "{0} not a valid database name. Please specify one of the"
                " following database names: {1}".format(value, valid_databases)
            )
        self._database = value

    def __init__(
        self,
        database: str,
        work_directory: pathlib.Path,
        path_to_case: pathlib.Path,
        path_to_simulation_mesh: pathlib.Path = None,
        path_to_model: pathlib.Path = None,
        mesh_size: float = 1.5,
        add_blood_pool: bool = False,
    ) -> None:
        self.database = database
        """Name of the database to use."""
        self.workdir = work_directory
        """Path to the working directory."""
        self.path_to_original_mesh = path_to_case
        """Path to the original mesh file."""
        self.path_to_simulation_mesh = path_to_simulation_mesh
        """Path to simulation(in .vtk format)."""
        self.path_to_model = path_to_model
        """Path to model (in .pickle format)."""
        self.labels_to_ids = LABELS_TO_ID[database]
        """Dict that maps labels > part/tag id."""
        self.ids_to_labels = dict((v, k) for k, v in LABELS_TO_ID[database].items())
        """Inverted dict that maps part/tag id > labels."""
        self.model_type: str = None
        """Model (geometric) type."""
        self.mesh_size: float = mesh_size
        """Mesh size used for remeshing."""
        self.add_blood_pool: bool = add_blood_pool
        """Flag indicating whether to add blood to the cavities."""

        pass

    def clean_workdir(
        self,
        extensions_to_remove: List[str] = [".stl", ".vtk", ".msh.h5"],
        remove_all: bool = False,
    ) -> None:
        """Remove files with extension present in the working directory.

        Parameters
        ----------
        extensions_to_remove : List[str], optional
            List of extensions to remove, by default [".stl", ".vtk", ".msh.h5"]
        remove_all: bool, optional
            Flag indicating whether to remove files with any extension.
            Keeps files/folder without extension
        """
        import glob as glob

        files = []
        if not remove_all:
            for ext in extensions_to_remove:
                files += glob.glob(os.path.join(self.workdir, "*" + ext))
        elif remove_all:
            files = glob.glob(os.path.join(self.workdir, "*.*"))

            for file in files:
                try:
                    os.remove(file)
                except:
                    LOGGER.debug(f"Unable to delete: {file}")
        return

    def create_workdir(self) -> None:
        """Create the working directory if it doesn't exist."""
        if not os.path.isdir(self.workdir):
            os.makedirs(self.workdir)
        return

    def dump_info(self, filename: pathlib.Path = None) -> None:
        """Dump model information to file."""
        if not filename:
            filename = os.path.join(self.workdir, "model_info.json")
        with open(filename, "w") as file:
            file.write(json.dumps(self.__dict__, indent=4))
        return


class HeartModel:
    """Parent class for heart models."""

    @property
    def parts(self) -> List[Part]:
        """Return list of parts."""
        parts = []
        for key, value in self.__dict__.items():
            attribute = getattr(self, key)
            if isinstance(attribute, Part):
                parts.append(attribute)
        return parts

    @property
    def part_names(self) -> List[str]:
        """Return list of part names."""
        part_names = []
        for part in self.parts:
            part_names.append(part.name)
        return part_names

    @property
    def part_ids(self) -> List[int]:
        """Return list of used part ids."""
        return [part.pid for part in self.parts]

    @property
    def cavities(self) -> List[Cavity]:
        """Return list of cavities in the model."""
        return [part.cavity for part in self.parts if part.cavity]

    def __init__(self, info: ModelInfo) -> None:
        self.info = info
        """Model meta information."""
        self.mesh = Mesh()
        """Modified mesh of the tissue."""
        self.mesh_raw = Mesh()
        """Raw input mesh."""

        self.fluid_mesh = Mesh()
        """Generated fluid mesh."""

        self._add_subparts()
        """Add any subparts."""
        self._add_labels_to_parts()
        """Add appropriate vtk labels to the parts."""

        if not self.info.mesh_size:
            self._set_default_mesh_size()
            """Set default mesh size."""

        self.model_type = self.__class__.__name__
        """Model type."""
        self.info.model_type = self.__class__.__name__

        """
        Elemental labels of AHA17 segments.
        NaN for elements not belong to left ventricle.
        """
        self.aha_ids = None

        self.cap_centroids: List[Point] = []
        """Centroid point to create cap shell."""

        self.electrodes: List[Point] = []
        """Electrodes positions for ECG computing."""

    def extract_simulation_mesh(self, clean_up: bool = False) -> None:
        """Update the model.

        Examples
        --------
        Processes a model from the public database and generates a simulation mesh

        1. Instantiate object that stores relevant information for the preprocessor

        >>> info = ModelInfo(
                       database="Strocchi2020",
                       path_to_case = "01.case",
                       work_directory = "workdir",
                       mesh_size = 1.5
                       )

        2. Instantiate the model object

        >>> model = FullHeart(info)

        3. Run the extract simulation mesh method

        >>> model.extract_simulation_mesh()

        4. Save model to disk

        >>> model.dump_model()

        5. Print info about the model

        >>> model.print_info()


        """
        # TODO factor in several methods
        if self.info.database == "LabeledSurface":
            input_surface = pv.read(self.info.path_to_original_mesh)

            points = input_surface.points
            try:  # Unstructure grid
                tris = input_surface.cells.reshape(-1, 4)[:, 1:]
            except:  # Polydata
                tris = input_surface.faces.reshape(-1, 4)[:, 1:]
            if isinstance(self, LeftVentricle):
                label_to_id = {
                    "Left ventricle endocardium": 1,
                    "Left ventricle epicardium": 2,
                    "Left ventricle myocardium mitral valve": 3,
                    "Left ventricle myocardium aortic valve": 4,
                }
                self.parts[0].tag_ids = [1]

            elif isinstance(self, BiVentricle):
                label_to_id = {
                    "Left ventricle endocardium": 1,
                    "Left ventricle epicardium": 2,
                    "Left ventricle myocardium mitral valve": 3,
                    "Left ventricle myocardium aortic valve": 4,
                    "Right ventricle endocardium": 5,
                    "Right ventricle epicardium": 6,
                    "Right ventricle myocardium pulmonary valve": 7,
                    "Right ventricle myocardium tricuspid valve": 8,
                    "Right ventricle septum": 9,
                }
                for part in self.parts:
                    if "Left" in part.name:
                        part.tag_ids = [1]
                    elif "Right" in part.name:
                        part.tag_ids = [2]
            else:
                print("Under development")
                exit()

            # create lists of relevant surfaces for remeshing
            # write stl input for volume meshing
            for key, tagid in label_to_id.items():
                surface_name = "-".join(key.lower().split())
                mask = input_surface["cell-tags"] == tagid
                surface = SurfaceMesh(surface_name, tris[mask, :], nodes=points)
                surface.write_to_stl(
                    os.path.join(self.info.workdir, "part_" + surface_name + ".stl")
                )
                if "valve" in key or "septum" in key:
                    self.mesh_raw.boundaries.append(surface)
                else:
                    self.mesh_raw.interfaces.append(surface)

            # Start remeshing
            self._remesh(post_mapping=False)

            if isinstance(self, BiVentricle):
                # split epicardium boundary into 'real' epicardium and septum
                for b in self.mesh.boundaries:
                    if "left-ventricle-epicardium" in b.name:
                        triangles = b.faces.reshape(-1, 4)[:, 1:]
                        region_ids = vtkmethods.get_connected_regions(b.nodes, triangles)
                        unique_region_ids, counts = np.unique(region_ids, return_counts=True)
                        sort_idx = np.argsort(counts)
                        assert len(unique_region_ids) == 2, "Expecting two regions"
                        unique_region_ids = unique_region_ids[sort_idx]
                        mask_septum = region_ids == unique_region_ids[0]
                        mask_epicardium = region_ids == unique_region_ids[1]
                        self.mesh.boundaries.append(
                            SurfaceMesh(
                                "right-ventricle-septum",
                                triangles=triangles[mask_septum, :],
                                nodes=b.nodes,
                            )
                        )
                        vertex = np.ones((len(triangles[mask_epicardium]), 1), dtype=int) * 3
                        bb = np.hstack((vertex, triangles[mask_epicardium, :]))
                        b.faces = bb.ravel()

                # write boundaries for debugging purposes
                # for b in model.mesh.boundaries:
                #     b.save(os.path.join(model.info.workdir, b.name + ".stl"))

                # update lv and rv parts (assign tet ids to these parts)
                vertex = np.linspace(0, len(points) - 1, num=len(points), dtype=int).reshape(-1, 1)
                meshio.write_points_cells(
                    "tags.vtk",
                    points,
                    [("vertex", vertex)],
                    cell_data={"tags": [input_surface["tags"]]},
                )
                vertex_tags = vtkmethods.read_vtk_unstructuredgrid_file("tags.vtk")
                os.remove("tags.vtk")
                self.mesh = vtkmethods.vtk_map_discrete_cell_data(vertex_tags, self.mesh, "tags")

            elif isinstance(self, LeftVentricle):
                self.mesh.cell_data["tags"] = np.ones(self.mesh.tetrahedrons.shape[0], dtype=float)
            # extract septum
            self._extract_septum()

            # call some additional methods to assign parts, surfaces, caps, etc
            self._assign_elements_to_parts()
            self._assign_surfaces_to_parts()
            self._assign_caps_to_parts()
            self._assign_cavities_to_parts()
            self._extract_apex()
            self._add_nodal_areas()

            def _define_uvc():
                """Temporal method use to compute longitudinal coordinates."""
                # add uvc coordinates to the self
                # 0. Compute longitudinal axis with: average cap centroid average apical point
                # 1. use longitudal axis to compute rotate to
                # xy plane and compute uvc_l by z-coordinate
                lv_apex = self.left_ventricle.apex_points[1].xyz
                mv_centroid = [
                    c.centroid for p in self.parts for c in p.caps if "mitral" in c.name
                ][0]
                longitudinal_axis = lv_apex - mv_centroid
                # np.savetxt(os.path.join(self.info.workdir, "points.txt"),
                #            np.vstack([lv_apex,mv_centroid]),delimiter=",")
                from ansys.heart.preprocessor.mesh.geodisc import rodrigues_rot

                Prot = rodrigues_rot(self.mesh.nodes - lv_apex, longitudinal_axis, [0, 0, -1])
                Prot[:, 2] = Prot[:, 2] - np.min(Prot, axis=0)[2]
                scaling = Prot[:, 2] / np.max(Prot[:, 2])
                # self.mesh.point_data = {}
                self.mesh.point_data["uvc_longitudinal"] = scaling
                self.mesh.cell_data["fiber"] = np.tile(
                    [0.0, 0.0, 1.0], (self.mesh.tetrahedrons.shape[0], 1)
                )
                self.mesh.cell_data["sheet"] = np.tile(
                    [0.0, 1.0, 0.0], (self.mesh.tetrahedrons.shape[0], 1)
                )

            _define_uvc()

            # rename cap
            for part in self.parts:
                for cap in part.caps:
                    cap.name = cap.name.replace("left-ventricle-myocardium-", "")
                    cap.name = cap.name.replace("right-ventricle-myocardium-", "")

        else:
            self.read_input_mesh()
            self._remove_unused_tags()
            self._prepare_for_meshing()
            self._remesh()
            self._update_parts()

        if clean_up:
            self.info.clean_workdir(["part*.stl", "cavity*.stl"])

        return

    def get_part(self, name: str, by_substring: bool = False) -> Union[Part, None]:
        """Get specific part based on part name."""
        found = False
        for part in self.parts:
            if part.name == name:
                return part
            if by_substring:
                if name in part.name:
                    return part
        if not found:
            return None

    def add_part(self, part_name: str) -> None:
        """Dynamically add a part as an attribute to the object."""
        setattr(self, "_".join(part_name.lower().split()), Part(name=part_name))
        return

    def remove_part(self, part_name: str) -> None:
        """Remove a part with a specific name from the model."""
        keys = self.__dict__.keys()
        for key in keys:
            attribute = getattr(self, key)
            if isinstance(attribute, Part):
                if part_name == attribute.name:
                    delattr(self, key)
                    return
        return

    def print_info(self) -> None:
        """Print information about the model."""
        if not isinstance(self.mesh.tetrahedrons, np.ndarray):
            LOGGER.info("Nothing to print")
            return

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
                    "\tsurface: {:} | # faces: {:d}".format(
                        surface.name, surface.triangles.shape[0]
                    )
                )
            for cap in part.caps:
                LOGGER.info("\tcap: {:} | # nodes {:d}".format(cap.name, len(cap.node_ids)))
            if part.cavity:
                LOGGER.info(
                    "\tcavity: {:} | volume: {:.1f} [mm3]".format(
                        part.cavity.name, part.cavity.surface.volume
                    )
                )
            LOGGER.info("-----------------------------------------")
        LOGGER.info("*****************************************")
        LOGGER.info("*****************************************")
        return

    def read_input_mesh(self) -> None:
        """Read the input mesh defined in ModelInfo."""
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

    def dump_model(self, filename: pathlib.Path = None, remove_raw_mesh: bool = True) -> None:
        """Save model to file.

        Examples
        --------
        >>> model.dump_model("my_heart_model.pickle")

        """
        LOGGER.debug("Writing model to disk")
        if not filename and not self.info.path_to_model:
            filename = os.path.join(self.info.workdir, "heart_model.pickle")
        elif not filename:
            filename = self.info.path_to_model
        else:
            self.info.path_to_model = filename

        # cleanup model object for more efficient storage
        # NOTE deleting faces, nodes of surfaces does not affect size
        # due to copy-by-reference
        if remove_raw_mesh:
            self.mesh_raw = None

        with open(filename, "wb") as file:
            pickle.dump(self, file)
        self.info.dump_info()
        return

    def plot_mesh(
        self, plot_raw_mesh: bool = False, show_edges: bool = True, color_by: str = "tags"
    ):
        """Plot the volume mesh of the heart model.

        Parameters
        ----------
        plot_raw_mesh : bool, optional
            Whether to plot the raw mesh, by default False
        show_edges : bool, optional
            Whether to plot the edges, by default True
        color_by : str, optional
            Color by cell/point data, by default "tags"

        Examples
        --------
        >>> import ansys.heart.preprocessor.models as models
        >>> model = models.HeartModel.load_model("heart_model.pickle")
        >>> model.plot_mesh(show_edges=True)
        """
        try:
            import pyvista
        except ImportError:
            LOGGER.warning("pyvista not found. Install with: pip install pyvista")
            return

        plotter = pyvista.Plotter()
        if plot_raw_mesh:
            plotter.add_mesh(self.mesh_raw, show_edges=show_edges, scalars="tags")
        else:
            plotter.add_mesh(self.mesh, show_edges=show_edges, scalars="tags")

        plotter.show()
        return

    def plot_fibers(self, n_seed_points: int = 1000, plot_raw_mesh: bool = False):
        """Plot the mesh and fibers as streamlines.

        Parameters
        ----------
        plot_raw_mesh : bool, optional
            Flag indicating whether to plot the streamlines on the raw mesh, by default False
        n_seed_points : int, optional
            Number of seed points. Recommended to use 5000, by default 1000

        Examples
        --------
        >>> import ansys.heart.preprocessor.models as models
        >>> model = models.HeartModel.load_model("my_model.pickle")
        >>> model.plot_fibers(n_seed_points=5000)
        """
        try:
            import pyvista
        except ImportError:
            LOGGER.warning("pyvista not found. Install with: pip install pyvista")
            return
        plotter = pyvista.Plotter()

        if plot_raw_mesh:
            if not isinstance(self.mesh_raw, Mesh):
                LOGGER.info("Raw mesh not available.")
                return
            mesh = self.mesh_raw
        else:
            mesh = self.mesh

        mesh = mesh.ctp()
        streamlines = mesh.streamlines(vectors="fiber", source_radius=75, n_points=n_seed_points)
        tubes = streamlines.tube()
        plotter.add_mesh(mesh, opacity=0.5, color="white")
        plotter.add_mesh(tubes, color="white")
        plotter.show()
        return

    def plot_surfaces(self, show_edges: bool = True):
        """Plot all the surfaces in the model.

        Examples
        --------
        Import modules and load model.
        >>> import ansys.heart.preprocessor.models as models
        >>> model = models.HeartModel.load_model("my_model.pickle")
        Plot the model
        >>> model.plot(show_edges=True)
        """
        try:
            import pyvista as pv
        except ImportError:
            LOGGER.warning(
                "PyVista not found: visualization not supported."
                "Install pyvista with: pip install pyvista"
            )
            return
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            LOGGER.warning("matplotlib not found. Install matplotlib with: pip install matplotlib")
            return

        surfaces_to_plot = [s for p in self.parts for s in p.surfaces]
        valves = [b for b in self.mesh.boundaries if "valve" in b.name or "border" in b.name]
        surfaces_to_plot = surfaces_to_plot + valves

        color_map = plt.cm.get_cmap("tab20", len(surfaces_to_plot))
        colors = color_map.colors[:, 0:3]
        plotter = pv.Plotter()
        ii = 0
        for surface in surfaces_to_plot:
            plotter.add_mesh(
                surface,
                color=colors[ii, :],
                show_edges=show_edges,
                label=surface.name,
            )
            plotter.add_legend(face=None, size=(0.25, 0.25), loc="lower left")
            ii += 1

        plotter.show()
        return

    def plot_purkinje(self):
        """Plot the mesh and Purkinje network."""
        if not len(self.mesh.beam_network) > 0:
            LOGGER.info("No Purkinje network to plot.")
            return

        try:
            import pyvista as pv
        except ImportError:
            LOGGER.warning(
                "PyVista not found: visualization not supported."
                "Install pyvista with: pip install pyvista"
            )
            return

        try:
            plotter = pv.Plotter()
            plotter.add_mesh(self.mesh, color="w", opacity=0.3)
            for beams in self.mesh.beam_network:
                plotter.add_mesh(beams, color=np.random.uniform(size=3), line_width=3)
            plotter.show()
        except:
            LOGGER.warning("Failed to plot mesh.")
        return

    @staticmethod
    def load_model(filename: pathlib.Path):
        """Load a preprocessed model from file.

        Examples
        --------
        >>> model = HeartModel.load_model("my_model.pickle")

        """
        # NOTE: need to suppress some vtk errors in pickled pyvista objects.
        # change the verbosity in the vtk loggerger and suppress the python logger.
        import copy
        import logging

        import vtk

        logger = copy.deepcopy(logging.getLogger("pyheart_global"))
        # setting propagate to False is workaround for VTK changing log behavior
        logger.propagate = False

        # to suppress vtk errors
        vtk_logger = vtk.vtkLogger
        vtk_logger.SetStderrVerbosity(vtk.vtkLogger.VERBOSITY_OFF)
        with open(filename, "rb") as file:
            model = pickle.load(file)

        # logger = logging.getLogger("pyheart_global.preprocessor")

        # logger.disabled = False
        # vtk_logger.SetStderrVerbosity(vtk.vtkLogger.VERBOSITY_1)
        return model

    def _set_default_mesh_size(self) -> None:
        """Set the default mesh size."""
        LOGGER.warning("No mesh size given setting default mesh size")
        self.info.mesh_size = 1.5
        return

    def _add_labels_to_parts(self) -> None:
        """Use model definitions to add corresponding vtk labels to the part."""
        for part in self.parts:
            if part.name == "Septum":
                continue
            part.tag_labels = HEART_PARTS[part.name]["VTKLabels"]
            # add tag ids
            part.tag_ids = []
            for tag_label in part.tag_labels:
                part.tag_ids.append(LABELS_TO_ID[self.info.database][tag_label])
        return

    def _add_subparts(self) -> None:
        """Add subparts to parts of type ventricle."""
        for part in self.parts:
            if part.part_type in ["ventricle"]:
                part._add_myocardium_part()
                if "Left ventricle" in part.name:
                    part._add_septum_part()
        return

    def _remove_unused_tags(self) -> None:
        """Extract only the tags of interest."""
        # collect all used tags
        tag_ids = []
        for part in self.parts:
            if not part.tag_ids:
                continue
            tag_ids.extend(part.tag_ids)
        self.mesh_raw.keep_elements_with_value(tag_ids, "tags")
        return

    def _get_used_element_ids(self) -> np.ndarray:
        """Return array of used element ids."""
        element_ids = np.empty(0, dtype=int)
        for part in self.parts:
            element_ids = np.append(element_ids, part.element_ids)

        return element_ids

    def _get_endo_epicardial_surfaces(self) -> None:
        """Get endo- and epicardial surfaces.

        Note
        ----
        Also obtains the septum in case of a BiVentricle, FourChamber or FullHeart model

        """
        surfaces_to_add: List[SurfaceMesh] = []
        orphan_surfaces: List[SurfaceMesh] = []
        separated_regions = []
        for boundary in self.mesh_raw.boundaries:
            if not "ventricle" in boundary.name and not "atrium" in boundary.name:
                continue

            if "left-ventricle-myocardium" in boundary.name and isinstance(
                self, (BiVentricle, FourChamber, FullHeart)
            ):
                boundary_name_suffix = ["epicardium", "endocardium", "septum"]
            else:
                boundary_name_suffix = ["epicardium", "endocardium"]

            num_expected_regions = len(boundary_name_suffix)

            LOGGER.debug("Extracting : {} from {}".format(boundary_name_suffix, boundary.name))
            boundary_conn = boundary.clean().connectivity()
            region_ids = boundary_conn.cell_data["RegionId"] + 1

            unique_regions, counts = np.unique(region_ids, return_counts=True)

            num_regions = len(unique_regions)

            # NOTE: number of cells do not a guarantee to distinguish between endo and epicardium
            # sort by bounding box volume instead.
            volumes = []
            surfaces: List[SurfaceMesh] = []
            for ii, region_id in enumerate(unique_regions):
                mask = region_ids == region_id
                surface = SurfaceMesh(
                    triangles=boundary.triangles[mask, :], nodes=self.mesh_raw.nodes
                )
                volumes.append(surface.compute_bounding_box()[1])
                surfaces.append(surface)

            # sort by volume of bounding box
            order = np.flip(np.argsort(volumes))
            surfaces = [surfaces[idx] for idx in order]

            # get list of orphan surfaces if any
            if num_expected_regions < num_regions:
                LOGGER.warning(
                    "Expecting {0} regions but getting {1} regions for boundary: {2}".format(
                        num_expected_regions, num_regions, boundary.name
                    )
                )
                for jj, surface in enumerate(surfaces):
                    LOGGER.warning(
                        "\t{0}: num_faces: {1}".format(jj + 1, surface.triangles.shape[0])
                    )
                LOGGER.warning(
                    "First {0} surfaces have largest bounding box".format(num_expected_regions)
                )
                num_orphans = num_regions - num_expected_regions
                orphan_surfaces.extend(surfaces[-num_orphans:])

            surfaces_to_use = surfaces[:num_expected_regions]
            for ii, surface in enumerate(surfaces_to_use):
                surface.name = boundary.name.replace("myocardium", boundary_name_suffix[ii])

            separated_regions.append(boundary.name)
            surfaces_to_add += surfaces_to_use

        # try to assign the orphan faces to any other connected boundary
        for orphan_surface in orphan_surfaces:
            orphan_surface.get_boundary_edges()
            num_connected_nodes = []
            surfaces_unseparated = [
                b for b in self.mesh_raw.boundaries if b.name not in separated_regions
            ]
            for surface in surfaces_unseparated:
                if surface.boundary_edges.shape[0] == 0:
                    get_boundary_edges = True
                else:
                    get_boundary_edges = False

                if get_boundary_edges:
                    surface.get_boundary_edges()

                num_connected_nodes.append(
                    np.sum(np.isin(orphan_surface.boundary_nodes, surface.boundary_nodes))
                )
            if np.any(num_connected_nodes) > 0:
                surface_to_copy_to = surfaces_unseparated[np.argmax(num_connected_nodes)]
                LOGGER.warning(
                    "Copying {0} orphan faces to {1}".format(
                        orphan_surface.triangles.shape[0], surface_to_copy_to.name
                    )
                )
                surface_to_copy_to.triangles = np.vstack(
                    [surface_to_copy_to.triangles, orphan_surface.triangles]
                )

            else:
                LOGGER.warning("Could not find suitable candidate surface - proceed with caution")
                orphan_surface.write_to_stl(os.path.join(self.info.workdir, "orphan_surface.stl"))
                for surface in surfaces_to_add:
                    surface.write_to_stl(os.path.join(self.info.workdir, surface.name + ".stl"))

                if orphan_surface.triangles.shape[0] < 5:
                    LOGGER.warning("Deleting orphan surface: consists of less than 5 faces")
                else:
                    LOGGER.warning(
                        "Could not find suitable candidate surface to merge "
                        "orphan faces into - proceed with caution"
                    )

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

    def _prepare_for_meshing(self) -> None:
        """Prepare the input for volumetric meshing with Fluent meshing.

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

    def _remesh(self, post_mapping: bool = True) -> None:
        """
        Use the generated files to remesh the surfaces and volume.

        Parameters
        ----------
        post_mapping: bool,True
            Mapping raw mesh data to new mesh
        """
        LOGGER.info("Remeshing volume...")
        path_mesh_file = os.path.join(self.info.workdir, "fluent_volume_mesh.msh.h5")

        if not self.info.mesh_size:
            LOGGER.warning("No mesh size set: setting to uniform size of 1.5 mm")
            self.info.mesh_size = 1.5

        mesher.mesh_heart_model_by_fluent(
            self.info.workdir,
            path_mesh_file,
            mesh_size=self.info.mesh_size,
            add_blood_pool=self.info.add_blood_pool,
            show_gui=False,
        )

        fluent_mesh = mesher.hdf5.FluentMesh()
        fluent_mesh.load_mesh(path_mesh_file)

        tissue_cell_zone = next(cz for cz in fluent_mesh.cell_zones if "heart-tet-cells" in cz.name)
        tetra_tissue = tissue_cell_zone.cells

        # update mesh object of the tissue
        self.mesh.tetrahedrons = tetra_tissue
        self.mesh.nodes = fluent_mesh.nodes

        # ensures normals pointing into the cavity
        # NOTE: not sure what determines the ordering when adding the blood pool
        # E.g. the normals of the endo AND epicardium are now pointing inwards with
        # this switch deactivated when adding add blood pool
        if self.info.add_blood_pool:
            flip_face_order = False
        else:
            flip_face_order = True

        for face_zone in fluent_mesh.face_zones:
            # don't create surfaces for interior face zones
            if "interior" in face_zone.name or face_zone.zone_type != 3:
                continue

            if flip_face_order:
                faces = face_zone.faces[:, [0, 2, 1]]
            else:
                faces = face_zone.faces

            face_zone_surface_mesh = SurfaceMesh(
                name=face_zone.name,
                triangles=faces,
                nodes=self.mesh.nodes,
                id=face_zone.id,
            )
            # face_zone_surface_mesh.write_to_stl(
            #     os.path.join(self.info.workdir, face_zone.name + ".stl")
            # )
            self.mesh.boundaries.append(face_zone_surface_mesh)

        # update mesh object of the fluid
        if self.info.add_blood_pool:
            fluid_cell_zones = [
                cz for cz in fluent_mesh.cell_zones if "heart-tet-cells" not in cz.name
            ]
            tetras_fluid = np.empty((0, 4), dtype=int)
            for cell_zone in fluid_cell_zones:
                tetras_fluid = np.vstack([tetras_fluid, cell_zone.cells])
            self.fluid_mesh.tetrahedrons = tetras_fluid
            self.fluid_mesh.nodes = fluent_mesh.nodes

        if post_mapping:
            self._map_data_to_remeshed_volume()

        return

    def _map_data_to_remeshed_volume(self) -> None:
        """Map the data from the original mesh to the remeshed mesh."""
        # get list of tag ids to keep for mapping
        mapper = self.info.labels_to_ids
        labels = [l for p in self.parts if p.tag_labels for l in p.tag_labels]
        tag_ids = [mapper[l] for l in labels if "valve" not in l and "inlet" not in l]

        self.mesh_raw.keep_elements_with_value(tag_ids, "tags")

        # get list of all arrays in original mesh
        array_names = list(self.mesh_raw.cell_data.keys()) + list(self.mesh_raw.point_data.keys())
        array_names1 = self.mesh_raw.array_names

        source = self.mesh_raw  # vtkmethods.read_vtk_unstructuredgrid_file(filename_original)
        target = self.mesh  # vtkmethods.read_vtk_unstructuredgrid_file(filename_remeshed)

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
            mapped_cell_data,
            mapped_point_data,
        ) = vtkmethods.get_tetra_info_from_unstructgrid(target)
        for key, value in mapped_cell_data.items():
            self.mesh.cell_data[key] = value
        for key, value in mapped_point_data.items():
            self.mesh.point_data[key] = value

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

        return

    def _add_nodal_areas(self):
        """Compute and add nodal areas to surface nodes."""
        LOGGER.debug("Adding nodal areas")
        for surface in self.mesh.boundaries:
            nodal_areas = vtkmethods.compute_surface_nodal_area_pyvista(surface)
            surface.point_data["nodal_areas"] = nodal_areas

        # compute nodal areas for explicitly named surfaces
        for part in self.parts:
            for surface in part.surfaces:
                nodal_areas = vtkmethods.compute_surface_nodal_area_pyvista(surface)
                surface.point_data["nodal_areas"] = nodal_areas

        # add nodal areas to volume mesh. Note that nodes can be part of
        # multiple surfaces - so we need to perform a summation.
        # interior nodes will have an area of 0.
        self.mesh.point_data["nodal_areas"] = np.zeros(self.mesh.nodes.shape[0])

        for surface in self.mesh.boundaries:
            self.mesh.point_data["nodal_areas"] += surface.point_data["nodal_areas"]
        self.mesh.write_to_vtk(os.path.join(self.info.workdir, "volume_nodal_areas.vtk"))

        return

    def _add_surface_normals(self):
        """Add surface normal as point data and cell data to all 'named' surfaces in the model."""
        LOGGER.debug("Adding normals to all 'named' surfaces")
        for part in self.parts:
            for surface in part.surfaces:
                surface_with_normals = surface.compute_normals(
                    cell_normals=True, point_normals=True
                )
                surface.cell_data["normals"] = surface_with_normals.cell_data["Normals"]
                surface.point_data["normals"] = surface_with_normals.point_data["Normals"]

        return

    def _update_parts(self):
        """Update the parts using the (re)meshed volume.

        Notes
        -----
        1. Extracts septum
        2. Updates Parts to include element ids of the respective part
        3. Assign surfaces to each part
        4. Extracts the closing caps
        5. Creates cavities
        6. Extracts apical points
        7. Computes left-ventricle axis
        8. Computes left-ventricle 17 segments
        9. Adds nodal areas
        10. Adds surface normals to boundaries
        """
        self._extract_septum()
        self._assign_elements_to_parts()
        self._assign_surfaces_to_parts()
        self._assign_caps_to_parts()
        self._assign_cavities_to_parts()
        self._extract_apex()
        #
        self.compute_left_ventricle_anatomy_axis()
        self.compute_left_ventricle_aha17()

        self._add_nodal_areas()
        self._add_surface_normals()

        return

    def _extract_septum(self) -> None:
        """Separate the septum elements from the left ventricle.

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
            surface_septum.triangles, iters=1
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

    def _extract_apex(self) -> None:
        """Extract the apex for both the endocardium and epicardium of each ventricle.

        Note
        ----
        Apex defined as the point furthest from the mid-point between caps/valves
        If this point is on the edge, another point of the same element will be picked.

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

                surface.get_boundary_edges()
                if np.any(surface.boundary_edges == apical_node_id):
                    # Apical node is on the edge, need to adjust
                    element_id = np.argwhere(np.any(surface.triangles == apical_node_id, axis=1))[
                        0
                    ][0]
                    triangle = surface.triangles[element_id, :]

                    # get another point on the same element
                    apical_node_id = triangle[
                        np.argwhere(
                            np.isin(
                                triangle,
                                surface.boundary_edges,
                                invert=True,
                            )
                        )[
                            0
                        ][0]
                    ]
                    LOGGER.warning(
                        "Initial apical point is on edge of {0}, a close point is picked".format(
                            surface.name,
                        )
                    )

                #     assign apex point
                ventricle.apex_points.append(
                    Point(
                        name="apex " + surface_substring,
                        node_id=apical_node_id,
                        xyz=surface.nodes[apical_node_id, :],
                    )
                )

        return

    def _assign_elements_to_parts(self) -> None:
        """Get the element ids of each part and assign these to the Part objects."""
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

        return

    def _assign_surfaces_to_parts(self) -> None:
        """Assign surfaces generated during remeshing to model parts."""
        for part in self.parts:
            for surface in part.surfaces:
                boundary_name = "-".join(surface.name.lower().split())
                boundary_surface = self.mesh.get_surface_from_name(boundary_name)
                if boundary_surface:
                    surface.triangles = boundary_surface.triangles
                    surface.nodes = boundary_surface.nodes
                else:
                    LOGGER.warning("Could not find matching surface for: {0}".format(surface.name))

        return

    def _assign_caps_to_parts(self, unique_mitral_tricuspid_valve=True) -> None:
        """
        Use connectivity to obtain cap boundaries and adds these to their respective parts.

        Parameters
        ----------
        unique_mitral_tricuspid_valve
        If True, mitral/tricuspid valves defined for ventricles are also used for atrium.
        """
        used_boundary_surface_names = [s.name for p in self.parts for s in p.surfaces]
        remaining_surfaces = list(set(self.mesh.boundary_names) - set(used_boundary_surface_names))
        remaining_surfaces1: List[SurfaceMesh] = []
        for surface in self.mesh.boundaries:
            if surface.name not in remaining_surfaces:
                continue
            surface.get_boundary_edges()
            remaining_surfaces1.append(surface)

        # find intersection between remaining surfaces and part surfaces
        # This will find the valve/cap nodes
        for part in self.parts:
            for surface in part.surfaces:
                # special treatment since a part of surface is defined in septum
                if surface.name == "Right ventricle endocardium":
                    surface.get_boundary_edges(append_triangles=part.surfaces[2].triangles)
                else:
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

                            if unique_mitral_tricuspid_valve and "atrium" in part.name:
                                if "mitral" in name_valve or "tricuspid" in name_valve:
                                    LOGGER.debug(
                                        f"{name_valve} has been create in ventricular parts."
                                    )
                                    # Create dummy cap (only name) and will be filled later
                                    part.caps.append(Cap(name=name_valve))
                                    break

                            cap = Cap(name=name_valve, node_ids=edge_group.edges[:, 0])
                            cap.centroid = np.mean(surf.nodes[cap.node_ids, :], axis=0)

                            # # tessellation 0 : pick a node and create segments
                            # cap.tessellate()
                            # p1 = surf.nodes[cap.triangles[:, 1]] - surf.nodes[cap.triangles[:, 0]]
                            # p2 = surf.nodes[cap.triangles[:, 2]] - surf.nodes[cap.triangles[:, 0]]

                            # tessellation 1 : add a center node
                            cap.centroid_id = len(self.mesh.nodes) + len(
                                self.cap_centroids
                            )  # center node ID, 0 based
                            self.cap_centroids.append(
                                Point(
                                    name=name_valve + "_center",
                                    xyz=cap.centroid,
                                    node_id=cap.centroid_id,
                                )
                            )

                            cap.tessellate(center_point_id=[cap.centroid_id])
                            p1 = surf.nodes[cap.triangles[:, 1],] - cap.centroid
                            p2 = surf.nodes[cap.triangles[:, 2],] - cap.centroid

                            # get approximate cavity centroid to check normal of cap
                            cavity_centroid = surface.compute_centroid()
                            normals = np.cross(p1, p2)
                            cap_normal = np.mean(normals, axis=0)
                            cap_normal = cap_normal / np.linalg.norm(cap_normal)
                            cap.normal = cap_normal
                            cap_centroid = np.mean(surf.nodes[cap.node_ids, :], axis=0)
                            d1 = np.linalg.norm(cap_centroid + cap_normal - cavity_centroid)
                            d2 = np.linalg.norm(cap_centroid - cap_normal - cavity_centroid)
                            if d1 > d2:
                                LOGGER.debug(
                                    "Flipping order of nodes on cap to ensure normal "
                                    "pointing inward"
                                )
                                cap.node_ids = np.flip(cap.node_ids)
                                # flip segments
                                cap.triangles[:, [1, 2]] = cap.triangles[:, [2, 1]]
                                cap.normal = cap.normal * -1

                            part.caps.append(cap)
                            LOGGER.debug("Cap: {0} closes {1}".format(name_valve, surface.name))
                            break
        if unique_mitral_tricuspid_valve:
            # replace caps of atria by caps of ventricle
            for part in self.parts:
                if not "atrium" in part.name:
                    continue
                for cap in part.caps:
                    # replace with cap in ventricle (mitral and tricuspid valve)
                    cap_ref = [
                        c
                        for p in self.parts
                        if "ventricle" in p.name
                        for c in p.caps
                        if c.name == cap.name
                    ]

                    if len(cap_ref) == 1:
                        cap.centroid_id = cap_ref[0].centroid_id
                        # note: flip order to make sure normal is pointing inwards
                        cap.node_ids = np.flip(cap_ref[0].node_ids)
                        # flip segments
                        cap.triangles = cap_ref[0].triangles[:, [0, 2, 1]]

                        LOGGER.debug(
                            "Replacing cap {0} of part{1}: with that of the ventricle".format(
                                cap.name, part.name
                            )
                        )

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
                    endocardium.triangles = np.vstack([endocardium.triangles, interface.triangles])

        return

    def _assign_cavities_to_parts(self) -> None:
        """Create cavities based on endocardium surfaces and cap definitions."""
        # rename septum to right ventricle endocardium septum
        if isinstance(self, (BiVentricle, FourChamber, FullHeart)):
            part = self.get_part("Right ventricle", True)
            for surface in part.surfaces:
                if "Right ventricle septum" in surface.name:
                    surface.name = surface.name.replace("septum", "endocardium septum")

        # construct cavities with endocardium and caps

        if len(self.cap_centroids) == 0:
            nodes = self.mesh.nodes
        else:
            # a center node for each cap has been created, add them into create the cavity
            nodes = np.vstack((self.mesh.nodes, np.zeros((len(self.cap_centroids), 3))))
            for cap_center in self.cap_centroids:
                nodes[cap_center.node_id] = cap_center.xyz

        for part in self.parts:
            if "atrium" not in part.name and "ventricle" not in part.name:
                continue
            cavity_faces = np.empty((0, 3), dtype=int)

            surfaces = [s for s in part.surfaces if "endocardium" in s.name]
            for surface in surfaces:
                cavity_faces = np.vstack([cavity_faces, surface.triangles])

            for cap in part.caps:
                cavity_faces = np.vstack([cavity_faces, cap.triangles])

            surface = SurfaceMesh(name=part.name + " cavity", triangles=cavity_faces, nodes=nodes)
            part.cavity = Cavity(surface=surface, name=part.name)
            part.cavity.compute_centroid()

            LOGGER.debug("Volume of cavity: {0} = {1}".format(part.cavity.name, part.cavity.volume))

            part.cavity.surface.write_to_stl(
                os.path.join(self.info.workdir, "-".join(part.cavity.surface.name.lower().split()))
            )

        return

    def compute_left_ventricle_anatomy_axis(self, first_cut_short_axis=0.2):
        """
        Compute the long and short axes of the left ventricle.

        Parameters
        ----------
        first_cut_short_axis: default=0.2, use to avoid cut on aortic valve
        """
        for cap in self.left_ventricle.caps:
            if cap.name == "mitral-valve":
                mv_center = cap.centroid
            elif cap.name == "aortic-valve":
                av_center = cap.centroid

        # apex is defined on epicardium
        for ap in self.left_ventricle.apex_points:
            if ap.name == "apex epicardium":
                apex = ap.xyz

        # 4CAV long axis across apex, mitral and aortic valve centers
        center = np.mean(np.array([av_center, mv_center, apex]), axis=0)
        normal = np.cross(av_center - apex, mv_center - apex)
        self.l4cv_axis = {"center": center, "normal": normal / np.linalg.norm(normal)}

        # short axis: from mitral valve center to apex
        sh_axis = apex - mv_center
        # the highest possible but avoid to cut aortic valve
        center = mv_center + first_cut_short_axis * sh_axis
        self.short_axis = {"center": center, "normal": sh_axis / np.linalg.norm(sh_axis)}
        # LOGGER.info("Plane of short axis:", self.short_axis)

        # 2CAV long axis: normal to 4cav axe and pass mv center and apex
        center = np.mean(np.array([mv_center, apex]), axis=0)
        p1 = center + 10 * self.l4cv_axis["normal"]
        p2 = mv_center
        p3 = apex
        normal = np.cross(p1 - p2, p1 - p3)
        self.l2cv_axis = {"center": center, "normal": normal / np.linalg.norm(normal)}

        return

    def compute_left_ventricle_aha17(self, seg=17, p_junction=None) -> None:
        """
        Compute AHA17 label for left ventricle elements.

        Parameters
        ----------
        seg ::  default 17, or 16 segments
        p_junction: use CASIS definition for the first cut
        """
        self.aha_ids = np.empty(len(self.mesh.tetrahedrons))
        self.aha_ids[:] = np.nan

        # get lv elements
        try:
            ele_ids = np.hstack((self.left_ventricle.element_ids, self.septum.element_ids))
        except:  # if no septum exists?
            ele_ids = np.hstack(self.left_ventricle.element_ids)

        # left ventricle elements center
        elem_center = np.mean(self.mesh.nodes[self.mesh.tetrahedrons[ele_ids]], axis=1)
        label = np.empty(len(elem_center))
        label[:] = np.nan

        # anatomical points
        for cap in self.left_ventricle.caps:
            if cap.name == "mitral-valve":
                mv_center = cap.centroid
        for apex in self.left_ventricle.apex_points:
            if "endocardium" in apex.name:
                apex_ed = apex.xyz
            elif "epicardium" in apex.name:
                apex_ep = apex.xyz

        # short axis
        short_axis = self.short_axis["normal"]
        p_highest = self.short_axis["center"]

        # define reference cut plane
        if p_junction is not None:
            # CASIS definition: LV and RV junction point
            vec = (p_junction - p_highest) / np.linalg.norm(p_junction - p_highest)
            axe_60 = R.from_rotvec(np.radians(90) * short_axis).apply(vec)
        else:
            # default: rotate 60 from long axis
            axe_60 = R.from_rotvec(np.radians(60) * short_axis).apply(self.l4cv_axis["normal"])

        axe_120 = R.from_rotvec(np.radians(60) * short_axis).apply(axe_60)
        axe_180 = -R.from_rotvec(np.radians(60) * short_axis).apply(axe_120)
        axe_45 = R.from_rotvec(np.radians(-15) * short_axis).apply(axe_60)
        axe_135 = R.from_rotvec(np.radians(90) * short_axis).apply(axe_45)

        p1_3 = 1 / 3 * (apex_ep - p_highest) + p_highest
        p2_3 = 2 / 3 * (apex_ep - p_highest) + p_highest

        for i, n in enumerate(elem_center):
            # This part contains valves, do not considered by AHA17
            if np.dot(n - p_highest, mv_center - p_highest) >= 0:
                continue
            # Basal: segment 1 2 3 4 5 6
            elif np.dot(n - p1_3, mv_center - p1_3) >= 0:
                if np.dot(n - p1_3, axe_60) >= 0:
                    if np.dot(n - p1_3, axe_120) >= 0:
                        if np.dot(n - p1_3, axe_180) >= 0:
                            label[i] = 6
                        else:
                            label[i] = 5
                    else:
                        label[i] = 1
                else:
                    if np.dot(n - p1_3, axe_180) <= 0:
                        if np.dot(n - p1_3, axe_120) >= 0:
                            label[i] = 4
                        else:
                            label[i] = 3
                    else:
                        label[i] = 2
            # Mid cavity: segment 7 8 9 10 11 12
            elif np.dot(n - p2_3, mv_center - p2_3) >= 0:
                if np.dot(n - p1_3, axe_60) >= 0:
                    if np.dot(n - p1_3, axe_120) >= 0:
                        if np.dot(n - p1_3, axe_180) >= 0:
                            label[i] = 12
                        else:
                            label[i] = 11
                    else:
                        label[i] = 7
                else:
                    if np.dot(n - p1_3, axe_180) <= 0:
                        if np.dot(n - p1_3, axe_120) >= 0:
                            label[i] = 10
                        else:
                            label[i] = 9
                    else:
                        label[i] = 8
            # Apical
            else:
                if seg == 17:
                    if np.dot(n - apex_ed, apex_ep - apex_ed) >= 0:
                        label[i] = 17
                    else:
                        if np.dot(n - p1_3, axe_45) >= 0:
                            if np.dot(n - p1_3, axe_135) >= 0:
                                label[i] = 16
                            else:
                                label[i] = 13
                        else:
                            if np.dot(n - p1_3, axe_135) >= 0:
                                label[i] = 15
                            else:
                                label[i] = 14

                else:
                    if np.dot(n - p1_3, axe_45) >= 0:
                        if np.dot(n - p1_3, axe_135) >= 0:
                            label[i] = 16
                        else:
                            label[i] = 13
                    else:
                        if np.dot(n - p1_3, axe_135) >= 0:
                            label[i] = 15
                        else:
                            label[i] = 14

        self.aha_ids[ele_ids] = label

        return

    def compute_left_ventricle_element_cs(self):
        """Compute elemental coordinate system for aha17 elements."""
        ele_ids = np.where(~np.isnan(self.aha_ids))[0]
        elems = self.mesh.tetrahedrons[ele_ids]
        elem_center = np.mean(self.mesh.nodes[elems], axis=1)

        # compute longitudinal direction, i.e. short axis
        e_l = np.tile(self.short_axis["normal"], (len(ele_ids), 1))

        # compute radial direction
        center_offset = elem_center - self.left_ventricle.apex_points[1].xyz
        e_r = center_offset - (np.sum(e_l * center_offset, axis=1) * e_l.T).T
        # normalize each row
        e_r /= np.linalg.norm(e_r, axis=1)[:, np.newaxis]

        # compute circumferential direction
        e_c = np.cross(e_l, e_r)

        return e_l, e_r, e_c

    def _compute_uvc_rotation_bc(self, mesh: pv.UnstructuredGrid):
        """Select node set on long axis plane."""
        mesh["cell_ids"] = np.arange(0, mesh.n_cells, dtype=int)
        mesh["point_ids"] = np.arange(0, mesh.n_points, dtype=int)
        slice = mesh.slice(origin=self.l4cv_axis["center"], normal=self.l4cv_axis["normal"])
        crinkled = mesh.extract_cells(np.unique(slice["cell_ids"]))
        free_wall_center, septum_center = crinkled.clip(
            origin=self.l2cv_axis["center"],
            normal=-self.l2cv_axis["normal"],
            crinkle=True,
            return_clipped=True,
        )

        rotation_mesh = mesh.remove_cells(free_wall_center["cell_ids"])
        print(f"{mesh.n_points - rotation_mesh.n_points} nodes are removed from clip.")

        vn = mesh.points[free_wall_center["point_ids"]] - self.l4cv_axis["center"]
        v0 = np.tile(self.l4cv_axis["normal"], (len(free_wall_center["point_ids"]), 1))

        dot = np.einsum("ij,ij->i", v0, vn)  # dot product row by row
        set1 = free_wall_center["point_ids"][dot >= 0]  # -pi
        set2 = free_wall_center["point_ids"][dot < 0]  # pi
        set3 = np.setdiff1d(septum_center["point_ids"], free_wall_center["point_ids"])  # 0

        # visu
        # mesh["bc"] = np.zeros(mesh.n_points)
        # mesh["bc"][set1] = 1
        # mesh["bc"][set2] = -1
        # mesh["bc"][set3] = 2
        # mesh.set_active_scalars("bc")
        # mesh.plot()

        return set1, set2, set3

    def _compute_uvc_apex_set(self, radius=3):
        """Todo Only for LeftVentricle."""
        import scipy.spatial as spatial

        point_tree = spatial.cKDTree(self.mesh.points)
        set = point_tree.query_ball_point(self.parts[0].apex_points[1].xyz, radius)

        # print(set)
        return np.array(set)


class LeftVentricle(HeartModel):
    """Model of just the left ventricle."""

    def __init__(self, info: ModelInfo) -> None:
        self.left_ventricle: Part = Part(name="Left ventricle", part_type="ventricle")
        """Left ventricle part."""
        # remove septum - not used in left ventricle only model
        del self.left_ventricle.septum

        super().__init__(info)
        pass


class BiVentricle(HeartModel):
    """Model of the left and right ventricle."""

    def __init__(self, info: ModelInfo) -> None:
        self.left_ventricle: Part = Part(name="Left ventricle", part_type="ventricle")
        """Left ventricle part."""
        self.right_ventricle: Part = Part(name="Right ventricle", part_type="ventricle")
        """Right ventricle part."""
        self.septum: Part = Part(name="Septum", part_type="septum")
        """Septum."""

        super().__init__(info)
        pass


class FourChamber(HeartModel):
    """Model of the left/right ventricle and left/right atrium."""

    def __init__(self, info: ModelInfo) -> None:
        self.left_ventricle: Part = Part(name="Left ventricle", part_type="ventricle")
        """Left ventricle part."""
        self.right_ventricle: Part = Part(name="Right ventricle", part_type="ventricle")
        """Right ventricle part."""
        self.septum: Part = Part(name="Septum", part_type="septum")
        """Septum."""

        self.left_atrium: Part = Part(name="Left atrium", part_type="atrium")
        """Left atrium part."""
        self.right_atrium: Part = Part(name="Right atrium", part_type="atrium")
        """Right atrium part."""

        super().__init__(info)

    def compute_SA_node(self) -> Point:
        """
        Compute SinoAtrial node.

        SinoAtrial node is defined on endocardium surface and
        between sup vena cave and inf vena cave.
        """
        right_atrium_endo = self.right_atrium.endocardium

        for cap in self.right_atrium.caps:
            if "superior" in cap.name:
                sup_vcava_centroid = cap.centroid
            elif "inferior" in cap.name:
                inf_vcava_centroid = cap.centroid

        # define SinoAtrial node at here:
        target_coord = sup_vcava_centroid - (inf_vcava_centroid - sup_vcava_centroid) / 2

        target_id = pv.PolyData(self.mesh.nodes[right_atrium_endo.node_ids, :]).find_closest_point(
            target_coord
        )
        SA_node_id = right_atrium_endo.node_ids[target_id]
        SA_point = Point(name="SA_node", xyz=self.mesh.nodes[SA_node_id, :], node_id=SA_node_id)
        self.right_atrium.points.append(SA_point)

        return SA_point

    def compute_AV_node(self) -> Point:
        """
        Compute Atrio-Ventricular node.

        AtrioVentricular node is defined on endocardium surface and closest to septum.

        Returns
        -------
        Point
            returns the AV node.
        """
        right_atrium_endo = self.right_atrium.endocardium

        for surface in self.right_ventricle.surfaces:
            if "endocardium" in surface.name and "septum" in surface.name:
                right_septum = surface

        # define AtrioVentricular as the closest point to septum
        target_coord = pv.PolyData(
            self.mesh.points[right_atrium_endo.node_ids, :]
        ).find_closest_point(pv.PolyData(self.mesh.points[right_septum.node_ids, :]).center)

        # assign a point
        target_id = right_atrium_endo.node_ids[target_coord]
        AV_point = Point(name="AV_node", xyz=self.mesh.nodes[target_id, :], node_id=target_id)
        self.right_atrium.points.append(AV_point)

        return AV_point

        #

        # Notes
        # -----
        # 1. with create_new_nodes=True, node ID of AV point will be modified
        # 2. todo: multiple paths

        # Parameters
        # ----------
        # create_new_nodes: if duplicate news from solid elements.

    def compute_av_conduction(self, create_new_nodes=True) -> BeamMesh:
        """Compute Atrio-Ventricular conduction by means of beams following a geodesic path.

        Parameters
        ----------
        create_new_nodes : bool, optional
            Duplicate nodes found of the computed geodesic path, by default True

        Returns
        -------
        BeamMesh
            Beam mesh.

        Raises
        ------
        NotImplementedError
            Not implemented error.
        """
        if not create_new_nodes:
            raise NotImplementedError

        right_atrium_endo = self.right_atrium.endocardium

        try:
            SA_node = self.right_atrium.get_point("SA_node").node_id
        except AttributeError:
            LOGGER.info("SA node is not defined, creating with default option.")
            SA_node = self.compute_SA_node().node_id

        try:
            AV_node = self.right_atrium.get_point("AV_node").node_id
        except AttributeError:
            LOGGER.info("AV node is not defined, creating with default option.")
            AV_node = self.compute_AV_node().node_id

        path_SAN_AVN = right_atrium_endo.geodesic(SA_node, AV_node)
        edges = path_SAN_AVN["vtkOriginalPointIds"]

        if create_new_nodes:
            # duplicate nodes inside the line, connect only SA node (the first) with 3D
            edges[1:] = len(self.mesh.nodes) + np.linspace(
                0, len(edges) - 2, len(edges) - 1, dtype=int
            )

            # modify ID of AV point
            self.right_atrium.get_point("AV_node").node_id = edges[-1]

            # build connectivity table
            edges = np.vstack((edges[:-1], edges[1:])).T

            # add to mesh
            beam = self.mesh.add_beam_network(
                new_nodes=path_SAN_AVN.points[1:, :], edges=edges, name="SAN_to_AVN"
            )

            return beam

    def compute_His_conduction(self, beam_number=4) -> BeamMesh:
        """Compute His bundle conduction.

        Parameters
        ----------
        beam_number : int, optional
            beam number, by default 4

        Returns
        -------
        BeamMesh
            Beam mesh
        """
        start_point, end_point = self._define_hisbundle_start_end_point(beam_number)

        # https://www.researchgate.net/publication/353154291_Morphometric_analysis_of_the_His_bundle_atrioven
        # tricular_fascicle_in_humans_and_other_animal_species_Histological_and_immunohistochemical_study
        # # (1.06 ± 0.6 mm)

        # create nodes from start to end
        new_nodes = np.array(
            [
                start_point.xyz,
                end_point.xyz,
            ]
        )

        step = (end_point.xyz - start_point.xyz) / (beam_number + 1)
        for i in range(1, beam_number):
            new_nodes = np.insert(new_nodes, i, start_point.xyz + i * step, axis=0)

        # path start from AV point, to septum start point, then to septum end point
        AV_node = self.right_atrium.get_point("AV_node")
        edges = np.concatenate(
            (
                np.array([AV_node.node_id]),
                len(self.mesh.nodes)
                + np.linspace(0, len(new_nodes) - 1, len(new_nodes), dtype=int),
            )
        )
        # create beam
        connect_table = np.vstack((edges[:-1], edges[1:])).T
        beam = self.mesh.add_beam_network(new_nodes=new_nodes, edges=connect_table, name="His")

        # TODO:
        # Find upper right septal point "URSP"
        # connect this point with AV node
        # get total length of bounding box
        # define HIS as fraction of line (few edge lengths) connecting URSP and septum centroid

        # plotter = pv.Plotter()
        # heartsurface = self.mesh.extract_surface()
        # plotter.add_mesh(heartsurface, opacity=0.01, color="w")
        # plotter.add_mesh(self.mesh.beam_network[0], opacity=0.1, line_width=4)
        # plotter.add_mesh(self.mesh.beam_network[1], opacity=0.1, color="k", line_width=4)
        # plotter.add_mesh(new_nodes[0, :], opacity=0.1, color="r")
        # plotter.add_mesh(new_nodes[0, :], opacity=0.1, color="r")
        # plotter.add_mesh(new_nodes[1, :], opacity=0.1, color="r")
        # plotter.add_mesh(self.mesh.nodes[AV_node.node_id, :], opacity=0.1, color="r")
        # plotter.add_mesh(self.mesh.nodes[His_septum_start_id, :], opacity=0.1, color="r")

        # plotter.show()

        return beam

    def _define_hisbundle_start_end_point(self, beam_number) -> (Point, Point):
        """
        Define start and end points of the bundle of His.

        Start point: create a point in septum, it's closest to AV node.
        End point: create a point in septum, it's close to AV node but randomly chosen.
        """
        AV_node = self.right_atrium.get_point("AV_node")

        septum_point_ids = np.unique(np.ravel(self.mesh.tetrahedrons[self.septum.element_ids,]))
        septum_pointcloud = pv.PolyData(self.mesh.nodes[septum_point_ids, :])

        # Define start point: closest to artria
        pointcloud_id = septum_pointcloud.find_closest_point(AV_node.xyz)

        His_septum_start_id = septum_point_ids[pointcloud_id]
        His_septum_start_xyz = self.mesh.nodes[His_septum_start_id, :]

        # Define end point:  a random close point
        n = 50
        pointcloud_id = septum_pointcloud.find_closest_point(AV_node.xyz, n=n)[n - 1]

        His_septum_end_id = septum_point_ids[pointcloud_id]
        His_septum_end_xyz = self.mesh.nodes[His_septum_end_id, :]

        # create Points
        His_septum_start = Point(
            name="His septum start", xyz=His_septum_start_xyz, node_id=len(self.mesh.nodes)
        )
        self.septum.points.append(His_septum_start)

        His_septum_end = Point(
            name="His septum end",
            xyz=His_septum_end_xyz,
            node_id=beam_number + len(self.mesh.nodes),
        )
        self.septum.points.append(His_septum_end)

        return His_septum_start, His_septum_end

    def compute_bundle_branches(self) -> (BeamMesh, BeamMesh):
        """Compute Buncle branches conduction system."""
        His_end = self.septum.get_point("His septum end")

        left_bundle = self._compute_bundle_oneside(
            His_end,
            self.left_ventricle.endocardium,
            self.left_ventricle.endocardium.node_ids,
            side="Left",
        )

        face = np.hstack(
            (self.right_ventricle.endocardium.faces, self.right_ventricle.septum.faces)
        )
        right_endo = pv.PolyData(self.mesh.points, face)

        right_endo_nids = np.unique(
            np.concatenate(
                (self.right_ventricle.septum.node_ids, self.right_ventricle.endocardium.node_ids)
            )
        )
        right_bundle = self._compute_bundle_oneside(
            His_end, right_endo, right_endo_nids, side="Right"
        )

        return left_bundle, right_bundle

    def _compute_bundle_oneside(self, His_end: Point, endo_surface, endo_nids, side: str):
        """Bundle brunch."""
        start_xyz = pv.PolyData(self.mesh.points[endo_nids, :]).find_closest_point(His_end.xyz)
        start_id = endo_nids[start_xyz]

        if side == "Left":
            ventricle = self.left_ventricle
        elif side == "Right":
            ventricle = self.right_ventricle

        bundle_branch = endo_surface.geodesic(
            endo_surface.find_closest_point(self.mesh.points[start_id, :]),
            endo_surface.find_closest_point(ventricle.apex_points[0].xyz),
        )

        # build branch net
        # exclude last (apex) node which belongs to purkinje beam
        new_nodes = bundle_branch.points[0:-1, :]

        # first node "his septum end" and last node "apex"
        edges = np.concatenate(([His_end.node_id], bundle_branch["vtkOriginalPointIds"]))

        # duplicate nodes except the two ends
        edges[1:-1] = len(self.mesh.nodes) + np.linspace(
            0, len(edges) - 3, len(edges) - 2, dtype=int
        )

        edges = np.vstack((edges[:-1], edges[1:])).T
        bundle_beam = self.mesh.add_beam_network(
            new_nodes=new_nodes, edges=edges, name=side + " bundle branch"
        )
        return bundle_beam

    def compute_Bachman_bundle(self):
        """Compute Bachman bundle conduction system."""
        raise NotImplementedError
        return

    def compute_cavity_interfaces(self):
        """Compute AtrioVentricular conduction system."""
        self.mesh.establish_connectivity()
        left_ventricle_left_atrium = []
        right_ventricle_right_atrium = []
        left_ventricle_right_atrium = []
        right_ventricle_left_atrium = []
        left_ventricle_left_atrium_name = "left-ventricle_left-atrium"
        right_ventricle_right_atrium_name = "right-ventricle_right-atrium"
        left_ventricle_right_atrium_name = "left-ventricle_right-atrium"
        right_ventricle_left_atrium_name = "right-ventricle_left-atrium"

        # build atrio-ventricular tag-id pairs
        # labels_to_ids stores the mapping between tag-ids and the corresponding label.
        labels_to_tag_ids = self.info.labels_to_ids
        left_ventricle_left_atrium = [
            labels_to_tag_ids["Left ventricle myocardium"],
            labels_to_tag_ids["Left atrium myocardium"],
        ]
        right_ventricle_right_atrium = [
            labels_to_tag_ids["Right ventricle myocardium"],
            labels_to_tag_ids["Right atrium myocardium"],
        ]
        left_ventricle_right_atrium = [
            labels_to_tag_ids["Left ventricle myocardium"],
            labels_to_tag_ids["Right atrium myocardium"],
        ]
        right_ventricle_left_atrium = [
            labels_to_tag_ids["Right ventricle myocardium"],
            labels_to_tag_ids["Left atrium myocardium"],
        ]

        # build atrioventricular tag_id pairs
        left_ventricle_left_atrium = np.unique(left_ventricle_left_atrium)
        right_ventricle_right_atrium = np.unique(right_ventricle_right_atrium)
        left_ventricle_right_atrium = np.unique(left_ventricle_right_atrium)
        right_ventricle_left_atrium = np.unique(right_ventricle_left_atrium)
        # find atrioventricular shared nodes/interfaces
        self.mesh.add_interfaces(
            [
                left_ventricle_left_atrium,
                right_ventricle_right_atrium,
                left_ventricle_right_atrium,
                right_ventricle_left_atrium,
            ],
            [
                left_ventricle_left_atrium_name,
                right_ventricle_right_atrium_name,
                left_ventricle_right_atrium_name,
                right_ventricle_left_atrium_name,
            ],
        )


class FullHeart(FourChamber):
    """Model of both ventricles, both atria, aorta and pulmonary artery."""

    def __init__(self, info: ModelInfo) -> None:
        self.left_ventricle: Part = Part(name="Left ventricle", part_type="ventricle")
        """Left ventricle part."""
        self.right_ventricle: Part = Part(name="Right ventricle", part_type="ventricle")
        """Right ventricle part."""
        self.septum: Part = Part(name="Septum", part_type="septum")
        """Septum."""
        self.left_atrium: Part = Part(name="Left atrium", part_type="atrium")
        """Left atrium part."""
        self.right_atrium: Part = Part(name="Right atrium", part_type="atrium")
        """Right atrium part."""

        self.aorta: Part = Part(name="Aorta", part_type="artery")
        """Aorta part."""
        self.pulmonary_artery: Part = Part(name="Pulmonary artery", part_type="artery")
        """Pulmonary artery part."""

        HeartModel.__init__(self, info)
        # super().__init__(info)

        pass


if __name__ == "__main__":
    info = ModelInfo(database="Strocchi2020", work_directory="tmp", path_to_case="test.case")

    model = LeftVentricle(info)
    print("LeftVentricle:")
    model.print_info()
    model.remove_part("Left ventricle")

    model = BiVentricle(info)
    print("BiVentricle:")
    model.print_info()

    model = FourChamber(info)
    print("FourChamber:")
    model.print_info()

    model = FullHeart(info)
    print("FullHeart:")
    model.print_info()

    print(model.part_names)
