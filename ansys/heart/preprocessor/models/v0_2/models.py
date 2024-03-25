"""Module containing classes for the various heart models."""

import copy
import json
import logging
import os

# import json
import pathlib
import pickle
from typing import List, Union

LOGGER = logging.getLogger("pyheart_global.preprocessor")
# from ansys.heart.preprocessor.input import HEART_MODELS
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
    _create_line,
)
import ansys.heart.preprocessor.mesh.vtkmethods as vtkmethods
from ansys.heart.preprocessor.models.v0_2.input import _InputModel
import numpy as np
import pyvista as pv
from scipy.spatial.transform import Rotation as R
import yaml


class ModelInfo:
    """Contains model information."""

    def __init__(
        self,
        input: Union[pathlib.Path, str, pv.PolyData, pv.UnstructuredGrid] = None,
        scalar: str = "part-id",
        part_definitions: dict = None,
        work_directory: pathlib.Path = ".",
        path_to_simulation_mesh: pathlib.Path = None,
        path_to_model: pathlib.Path = None,
        mesh_size: float = 1.5,
        add_blood_pool: bool = False,
    ) -> None:
        self.input = input
        """Input to the workflow."""
        self.scalar = scalar
        """Scalar field name with part/surface ids."""
        self.part_definitions = part_definitions
        """Part definitions."""

        self.workdir = work_directory
        """Path to the working directory."""
        self.path_to_simulation_mesh = path_to_simulation_mesh
        """Path to simulation(in .vtk format)."""
        self.path_to_model = path_to_model
        """Path to model (in .pickle format)."""

        self.mesh_size: float = mesh_size
        """Mesh size used for remeshing."""
        self.add_blood_pool: bool = add_blood_pool
        """Flag indicating whether to add blood to the cavities."""

        self.model_type: str = ""
        """Deprecated dummy value of model type."""

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
        if not isinstance(self.input, (str, pathlib.Path)):
            self.input = None

        if not filename:
            filename = os.path.join(self.workdir, "model_info.json")

        extension = os.path.splitext(filename)[1]
        with open(filename, "w") as file:
            if extension == ".json":
                formatted_string = json.dumps(self.__dict__, indent=4)
            elif extension == ".yml":
                formatted_string = yaml.dump(self.__dict__, indent=4, sort_keys=False)

            file.write(formatted_string)

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
    def surfaces(self) -> List[SurfaceMesh]:
        """Return list of all defined surfaces."""
        return [s for p in self.parts for s in p.surfaces]

    @property
    def surface_names(self) -> List[str]:
        """Return list of all defined surface names."""
        return [s.name for s in self.surfaces]

    @property
    def surface_ids(self) -> List[str]:
        """Return list of all defined surface names."""
        return [s.id for s in self.surfaces]

    @property
    def cavities(self) -> List[Cavity]:
        """Return list of cavities in the model."""
        return [part.cavity for part in self.parts if part.cavity]

    @property
    def part_name_to_part_id(self) -> dict:
        """Dictionary that maps the part name to the part id."""
        return {p.name: p.pid for p in self.parts}

    @property
    def part_id_to_part_name(self) -> dict:
        """Dictionary that maps part id to part name."""
        return {p.pid: p.name for p in self.parts}

    @property
    def surface_name_to_surface_id(self) -> dict:
        """Dictionary that maps surface name to surface id."""
        return {s.name: s.id for p in self.parts for s in p.surfaces}

    @property
    def surface_id_to_surface_name(self) -> dict:
        """Dictionary that maps surface name to surface id."""
        return {s.id: s.name for p in self.parts for s in p.surfaces}

    @property
    def cap_centroids(self):
        """Return list of cap centroids."""
        return [
            Point(name=c.name + "_center", xyz=c.centroid, node_id=c.centroid_id)
            for p in self.parts
            for c in p.caps
        ]

    def __init__(self, info: ModelInfo) -> None:
        self.info = info
        """Model meta information."""
        self.mesh = Mesh()
        """Computational mesh."""

        self.fluid_mesh = Mesh()
        """Generated fluid mesh."""

        self._input = _InputModel()
        """Input model."""

        self._add_subparts()
        """Add any subparts."""

        self._set_part_ids()
        """Set incremenetal part ids."""

        self.model_type = self.__class__.__name__
        """Model type."""

        self.aha_ids = None
        """American Heart Association ID's."""

        self.electrodes: List[Point] = []
        """Electrodes positions for ECG computing."""

        self.beam_network: List[BeamMesh] = []
        """List of beam networks in the mesh."""

        self.electrodes: List[Point] = []
        """Electrodes positions for ECG computing."""

        return

    def create_part_by_ids(self, eids: List[int], name: str) -> Union[None, Part]:
        """Create a new part by element ids.

        Parameters
        ----------
        eids : List[int]
            element id list
        name : str
            part name

        Returns
        -------
        Union[None, Part]
           return the part if succeed
        """
        if len(eids) == 0:
            LOGGER.error(f"Element list is empty to create {name}")
            return None

        for part in self.parts:
            try:
                part.element_ids = np.setdiff1d(part.element_ids, eids)
            except:
                LOGGER.error(f"Failed to create part {name}")
                return None

        self.add_part(name)
        new_part = self.get_part(name)
        new_part.part_type = "myocardium"
        new_part.element_ids = eids

        return new_part

    def add_purkinje_from_kfile(self, filename: pathlib.Path, name: str) -> None:
        """Read an LS-DYNA file containing purkinje beams and nodes.

        Parameters
        ----------
        filename : pathlib.Path

        name : str
            beamnet name
        """
        # Open file and import beams and created nodes
        with open(filename, "r") as file:
            start_nodes = 0
            lines = file.readlines()
        # find line ids delimiting node data and edge data
        start_nodes = np.array(np.where(["*NODE" in line for line in lines]))[0][0]
        end_nodes = np.array(np.where(["*" in line for line in lines]))
        end_nodes = end_nodes[end_nodes > start_nodes][0]
        start_beams = np.array(np.where(["*ELEMENT_BEAM" in line for line in lines]))[0][0]
        end_beams = np.array(np.where(["*" in line for line in lines]))
        end_beams = end_beams[end_beams > start_beams][0]

        # load node data
        node_data = np.loadtxt(
            filename, skiprows=start_nodes + 1, max_rows=end_nodes - start_nodes - 1
        )
        new_ids = node_data[:, 0].astype(int) - 1
        beam_nodes = node_data[:, 1:4]

        # load beam data
        beam_data = np.loadtxt(
            filename, skiprows=start_beams + 1, max_rows=end_beams - start_beams - 1, dtype=int
        )
        edges = beam_data[:, 2:4] - 1
        pid = beam_data[0, 1]

        # TODO: physically, this is not fully understood: Merging the end of bundle branch, the
        # origin of Purkinje and the apex of myiocardium seems logical, but it has more chance
        # the EP wave will not be triggered.
        # so I remove it, it means: end of bundle branch connect to apex, origin of Purkinje
        # is another point on the same location.

        # # replace origin (new created beam mesh) of purkinje by apex point (on solid mesh)
        # if "left" in name.lower():
        #     apex = self.left_ventricle.apex_points[0]
        # elif "right" in name.lower():
        #     apex = self.right_ventricle.apex_points[0]

        # to_be_replaced_id = new_ids[0]
        # new_ids = new_ids[1:]  # remove this point
        # beam_nodes = beam_nodes[1:]  # remove this point
        # edges[edges == to_be_replaced_id] = apex.node_id  # replace by apex Id

        #
        mask = np.isin(edges, new_ids)  # True for new created nodes
        edges[mask] -= new_ids[0]  # beam nodes id start from 0

        beam = self.add_beam_net(beam_nodes, edges, mask, pid=pid, name=name)

        return beam

    def add_beam_net(
        self, beam_nodes: np.ndarray, edges: np.ndarray, mask: np.ndarray, pid=0, name: str = None
    ) -> BeamMesh:
        """Add a BeamMesh object on the model.

        Parameters
        ----------
        beam_nodes : np.ndarray
            new nodes coordinates.
        edges : np.ndarray
            beam elements connectivity
            If `mask` is true, it's Id of `beam_nodes` (start by 0),
            it will be offset when creating BeamMesh object.
            If `mask` is false, it's Id of existed nodes, it will not be offset.
        mask : np.ndarray
            with the same shape of `edges`
        pid : int, optional
            part Id, will be reassigned when writing, by default 0
        name : str, optional
            name, by default None

        Returns
        -------
        BeamMesh
            BeamMesh object
        """
        edges[mask] += len(self.mesh.nodes) + len(BeamMesh.all_beam_nodes)

        if len(BeamMesh.all_beam_nodes) == 0:
            BeamMesh.all_beam_nodes = beam_nodes
        else:
            BeamMesh.all_beam_nodes = np.vstack((BeamMesh.all_beam_nodes, beam_nodes))

        # nodes is just for pyvista plot, edges used in writer will be offsetted
        # TODO only save necessary nodes, cells, and with a 'global id' array
        beam_net = BeamMesh(
            nodes=np.vstack((self.mesh.nodes, BeamMesh.all_beam_nodes)),
            edges=edges,
            beam_nodes_mask=mask,
        )
        beam_net.pid = pid
        beam_net.name = name

        # class variable BeamMesh.all_beam_nodes is not saved in pickle
        #
        # only the last created beam_network holds all previously created
        beam_net._all_beam_nodes = np.copy(BeamMesh.all_beam_nodes)

        self.beam_network.append(beam_net)

        # sync across beam networks
        for beams in self.beam_network:
            beams._all_beam_nodes = np.copy(BeamMesh.all_beam_nodes)

        # # visualize (debug)
        # import pyvista

        # plotter = pyvista.Plotter()
        # plotter.add_mesh(self.mesh, opacity=0.3)
        # plotter.add_mesh(beam_net)
        # plotter.show()

        return beam_net

    def load_input(self):
        """Use the content in model info to load the input model."""
        self._input = _InputModel(
            part_definitions=self.info.part_definitions,
            input=self.info.input,
            scalar=self.info.scalar,
        )
        return

    def mesh_volume(self, use_wrapper: bool = False):
        """Remesh the input model and fill the volume.

        Parameters
        ----------
        use_wrapper : bool, optional
            Flag for switch to non-manifold mesher, by default False

        Notes
        -----
        Note that when the input surfaces are non-manifold the wrapper tries
        to reconstruct the surface and parts. Inevitably this leads to
        reconstruction errors. Nevertheless, in some instances this approach is
        robuster than meshing from a manifold surface.
        """
        path_to_output_model = os.path.join(self.info.workdir, "simulation_mesh.msh.h5")

        if use_wrapper:
            LOGGER.warning("Meshing from non-manifold model not yet available.")

            fluent_mesh = mesher.mesh_from_non_manifold_input_model(
                model=self._input,
                workdir=self.info.workdir,
                mesh_size=self.info.mesh_size,
                path_to_output=path_to_output_model,
            )
        else:
            fluent_mesh = mesher.mesh_from_manifold_input_model(
                model=self._input,
                workdir=self.info.workdir,
                mesh_size=self.info.mesh_size,
                path_to_output=path_to_output_model,
            )

        # remove empty cell zones
        num_cell_zones1 = len(fluent_mesh.cell_zones)
        fluent_mesh.cell_zones = [cz for cz in fluent_mesh.cell_zones if cz.cells.shape[0] > 0]
        num_cell_zones2 = len(fluent_mesh.cell_zones)
        if num_cell_zones1 > num_cell_zones2:
            LOGGER.warning("Removed {0} cell zones".format(num_cell_zones1 - num_cell_zones2))

        # Use only cell zones that are inside the parts defined in the input.
        fluent_mesh.cell_zones = [
            cz for cz in fluent_mesh.cell_zones if cz.id in self._input.part_ids
        ]

        vtk_grid = fluent_mesh._to_vtk()

        mesh = Mesh(vtk_grid)
        mesh.cell_data["part-id"] = mesh.cell_data["cell-zone-ids"]

        # merge some face zones that Fluent split based on connectivity
        idx_to_remove = []
        for ii, fz in enumerate(fluent_mesh.face_zones):
            if ":" in fz.name:
                basename = fz.name.split(":")[0]
                ref_facezone = next(fz1 for fz1 in fluent_mesh.face_zones if fz1.name == basename)
                LOGGER.debug("Merging {0} with {1}".format(fz.name, ref_facezone.name))
                ref_facezone.faces = np.vstack([ref_facezone.faces, fz.faces])
                idx_to_remove += [ii]

        # remove merged face zone
        fluent_mesh.face_zones = [
            fz for ii, fz in enumerate(fluent_mesh.face_zones) if ii not in idx_to_remove
        ]

        # add face zones to mesh object.
        mesh.boundaries = [
            SurfaceMesh(name=fz.name, triangles=fz.faces, nodes=mesh.nodes, id=fz.id)
            for fz in fluent_mesh.face_zones
            if not "interior" in fz.name
        ]

        self.mesh = mesh

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

    def dump_model(self, filename: pathlib.Path = None) -> None:
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

        with open(filename, "wb") as file:
            pickle.dump(self, file)
        self.info.dump_info()
        return

    def plot_mesh(self, show_edges: bool = True, color_by: str = "part-id"):
        """Plot the volume mesh of the heart model.

        Parameters
        ----------
        show_edges : bool, optional
            Whether to plot the edges, by default True
        color_by : str, optional
            Color by cell/point data, by default "part-id"

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
        plotter.add_mesh(self.mesh, show_edges=show_edges, scalars=color_by)

        plotter.show()
        return

    def plot_part(self, part: Part):
        """Plot a part in mesh.

        Parameters
        ----------
        part : Part
            part to highlight in mesh

        Examples
        --------
        >>> import ansys.heart.preprocessor.models as models
        >>> model = models.HeartModel.load_model("my_model.pickle")
        >>> model.part(model.left_ventricle)
        """
        try:
            import pyvista
        except ImportError:
            LOGGER.warning("pyvista not found. Install with: pip install pyvista")
            return

        mesh = self.mesh

        plotter = pyvista.Plotter()
        plotter.add_mesh(mesh, opacity=0.5, color="white")
        part = mesh.extract_cells(part.element_ids)
        plotter.add_mesh(part, opacity=0.95, color="red")
        plotter.show()
        return

    def plot_fibers(self, n_seed_points: int = 1000):
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
        if not len(self.beam_network) > 0:
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
            for beams in self.beam_network:
                plotter.add_mesh(beams, color="r", line_width=2)
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
        # change the verbosity in the vtk logger and suppress the python logger.
        import logging

        import vtk

        logger = copy.deepcopy(logging.getLogger("pyheart_global"))
        # setting propagate to False is workaround for VTK changing log behavior
        logger.propagate = False

        logger = logging.getLogger()
        logger.disabled = True
        # to suppress vtk errors
        vtk_logger = vtk.vtkLogger
        vtk_logger.SetStderrVerbosity(vtk.vtkLogger.VERBOSITY_OFF)
        with open(filename, "rb") as file:
            model = pickle.load(file)
        logger.disabled = False
        vtk_logger.SetStderrVerbosity(vtk.vtkLogger.VERBOSITY_1)
        return model

    def _set_part_ids(self):
        """Populate part ids."""
        c = 1
        for p in self.parts:
            p.pid = c
            c += 1

    def _add_subparts(self) -> None:
        """Add subparts to parts of type ventricle."""
        for part in self.parts:
            if part.part_type in ["ventricle"]:
                part._add_myocardium_part()
                if "Left ventricle" in part.name:
                    part._add_septum_part()
        return

    def _get_used_element_ids(self) -> np.ndarray:
        """Return array of used element ids."""
        element_ids = np.empty(0, dtype=int)
        for part in self.parts:
            element_ids = np.append(element_ids, part.element_ids)

        return element_ids

    def _add_nodal_areas(self):
        """Compute and add nodal areas to surface nodes."""
        LOGGER.debug("Adding nodal areas")
        for surface in self.mesh.boundaries:
            nodal_areas = vtkmethods.compute_surface_nodal_area_pyvista(surface)
            surface.point_data["nodal_areas"] = nodal_areas

        # compute nodal areas for explicitly named surfaces
        for part in self.parts:
            for surface in part.surfaces:
                if surface.n_points == 0:
                    LOGGER.debug(f"Skipping {surface.name} no points in surface")
                    continue
                if surface.n_cells == 0:
                    LOGGER.debug(f"Skipping {surface.name} no cells in surface.")
                    continue

                nodal_areas = vtkmethods.compute_surface_nodal_area_pyvista(surface)
                surface.point_data["nodal_areas"] = nodal_areas

        # add nodal areas to volume mesh. Note that nodes can be part of
        # multiple surfaces - so we need to perform a summation.
        # interior nodes will have an area of 0.
        self.mesh.point_data["nodal_areas"] = np.zeros(self.mesh.nodes.shape[0])
        for surface in self.mesh.boundaries:
            self.mesh.point_data["nodal_areas"] += surface.point_data["nodal_areas"]
        # self.mesh.write_to_vtk(os.path.join(self.info.workdir, "volume_nodal_areas.vtk"))
        return

    def _add_surface_normals(self):
        """Add surface normal as point data and cell data to all 'named' surfaces in the model.

        Notes
        -----
        Note that we need to flip the normals due to the convention that Fluent Meshing uses.
        That is, normals point into the meshed domain.
        """
        LOGGER.debug("Adding normals to all 'named' surfaces")
        LOGGER.warning("Flipping normals.")
        for part in self.parts:
            for surface in part.surfaces:
                surface_with_normals = surface.compute_normals(
                    cell_normals=True, point_normals=True, inplace=True, flip_normals=True
                )
                surface.cell_data["normals"] = surface_with_normals.cell_data["Normals"]
                surface.point_data["normals"] = surface_with_normals.point_data["Normals"]

        return

    def _update_parts(self):
        """Update the parts using the meshed volume.

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
        self._sync_input_parts_to_model_parts()

        self._extract_septum()
        self._assign_elements_to_parts()
        self._assign_surfaces_to_parts()

        self._validate_parts()
        self._validate_surfaces()

        self._assign_caps_to_parts()
        self._validate_cap_names()

        self._add_nodal_areas()
        self._add_surface_normals()

        self._assign_cavities_to_parts()
        self._extract_apex()

        self.compute_left_ventricle_anatomy_axis()
        self.compute_left_ventricle_aha17()

        if "fiber" not in self.mesh.array_names:
            LOGGER.debug("Adding placeholder for fiber direction.")
            fiber = np.tile([[0.0, 0.0, 1.0]], (self.mesh.n_cells, 1))
            self.mesh.cell_data["fiber"] = fiber

        if "sheet" not in self.mesh.array_names:
            LOGGER.debug("Adding placeholder for sheet direction.")
            sheet = np.tile([[0.0, 1.0, 1.0]], (self.mesh.n_cells, 1))
            self.mesh.cell_data["sheet"] = sheet

        if "uvc_l" not in self.mesh.array_names:
            LOGGER.debug("Add approximate longitudinal coordinates.")
            lv_apex = self.left_ventricle.apex_points[1].xyz
            mv_centroid = [c.centroid for p in self.parts for c in p.caps if "mitral" in c.name][0]
            longitudinal_axis = lv_apex - mv_centroid
            from ansys.heart.preprocessor.mesh.geodisc import rodrigues_rot

            Prot = rodrigues_rot(self.mesh.nodes - lv_apex, longitudinal_axis, [0, 0, -1])
            Prot[:, 2] = Prot[:, 2] - np.min(Prot, axis=0)[2]
            scaling = Prot[:, 2] / np.max(Prot[:, 2])
            self.mesh.point_data["uvc_longitudinal"] = scaling

        return

    def _sync_input_parts_to_model_parts(self):
        """Synchronize the input parts to the model parts.

        Note
        ----
        Checks:
            overwrites the default part ids by those given by user.
        """
        # unassign any part ids.
        for p in self.parts:
            p.pid = None

        for input_part in self._input.parts:
            try:
                idx = self.part_names.index(input_part.name)
                self.parts[idx].pid = input_part.id
            except ValueError:
                LOGGER.debug(f"Failed to find a match for: {input_part.name}")
                continue

        # assign new ids to parts without part id
        for p in self.parts:
            if not p.pid:
                p.pid = max([pid for pid in self.part_ids if pid != None]) + 1
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
        if len(surface_septum) > 1:
            raise ValueError("Expecting only one surface that contains string: 'septum'")
        if len(surface_septum) == 0:
            raise ValueError("No boundary with name: 'septum' found")
        surface_septum = surface_septum[0]

        # extrude septum surface
        faces_septum = connectivity.remove_triangle_layers_from_trimesh(
            surface_septum.triangles, iters=1
        )

        septum_surface_vtk = vtkmethods.create_vtk_surface_triangles(self.mesh.nodes, faces_septum)

        septum_surface_vtk = vtkmethods.smooth_polydata(septum_surface_vtk)

        septum_surface_vtk_extruded = vtkmethods.extrude_polydata(septum_surface_vtk, 20)

        filename_vtk = os.path.join(self.info.workdir, "volume_mesh.vtk")
        self.mesh.write_to_vtk(filename_vtk)
        volume_vtk = vtkmethods.read_vtk_unstructuredgrid_file(filename_vtk)

        element_ids_septum = vtkmethods.cell_ids_inside_enclosed_surface(
            volume_vtk, septum_surface_vtk_extruded
        )

        # assign to septum
        part = next(part for part in self.parts if part.name == "Septum")
        part.element_ids = element_ids_septum
        self.mesh.cell_data["part-id"][element_ids_septum] = part.pid
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

            element_ids = np.where(np.isin(self.mesh.part_ids, part.pid))[0]
            element_ids = element_ids[np.isin(element_ids, used_element_ids, invert=True)]
            part.element_ids = element_ids

        summ = 0
        for part in self.parts:
            LOGGER.debug("Num elements in {0}: {1}".format(part.name, part.element_ids.shape[0]))
            summ = summ + part.element_ids.shape[0]
        LOGGER.debug("Total num elements: {}".format(summ))

        # if summ != self.mesh.tetrahedrons.shape[0]:
        LOGGER.debug(
            "{0} elements assigned to parts - {1} exist in mesh".format(
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
                if "septum" in surface.name:
                    try:
                        septum = [b for b in self.mesh.boundaries if "septum" in b.name]
                        if len(septum) > 1:
                            LOGGER.warning(
                                "Multiple candidate boundaries for septum found, using first one."
                            )
                        boundary_surface = septum[0]
                    except:
                        boundary_surface = None
                if boundary_surface:
                    surface.triangles = boundary_surface.triangles
                    surface.nodes = boundary_surface.nodes
                    surface.id = boundary_surface.id
                else:
                    LOGGER.warning("Could not find matching surface for: {0}".format(surface.name))

        return

    def _assign_caps_to_parts(self, unique_mitral_tricuspid_valve=True) -> None:
        """Use connectivity to obtain cap boundaries and add these to their respective parts.

        Parameters
        ----------
        unique_mitral_tricuspid_valve : bool, optional
            Unique mitral/tricuspid valve defined, by default True
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
                        LOGGER.warning(
                            "Assuming closed edge group: cavity with {0} not be watertight!".format(  # noqa: E501
                                surface.name
                            )
                        )
                        # continue

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

                            # add cap centroid.
                            cap.centroid = np.mean(surf.nodes[cap.node_ids, :], axis=0)

                            # add cap centroid to node list
                            self.mesh.nodes = np.vstack([self.mesh.nodes, cap.centroid])
                            surf.nodes = self.mesh.nodes
                            cap.centroid_id = self.mesh.nodes.shape[0] - 1

                            # get approximate cavity centroid to check normal of cap
                            cavity_centroid = surface.compute_centroid()

                            cap.tessellate1(use_centroid=True)

                            p1 = surf.nodes[cap.triangles[:, 1],] - surf.nodes[cap.triangles[:, 0],]
                            p2 = surf.nodes[cap.triangles[:, 2],] - surf.nodes[cap.triangles[:, 0],]
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
                                cap.tessellate1(use_centroid=True)
                                cap.normal = cap.normal * -1

                            self.mesh._sync_nodes_of_surfaces()

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
                    cap.centroid = cap_ref[0].centroid
                    cap.centroid_id = cap_ref[0].centroid_id
                    cap.tessellate1(use_centroid=True)

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
        for part in self.parts:
            if "atrium" not in part.name and "ventricle" not in part.name:
                continue
            cavity_faces = np.empty((0, 3), dtype=int)

            surfaces = [s for s in part.surfaces if "endocardium" in s.name]
            for surface in surfaces:
                cavity_faces = np.vstack([cavity_faces, surface.triangles])

            for cap in part.caps:
                cavity_faces = np.vstack([cavity_faces, cap.triangles])

            surface = SurfaceMesh(
                name=part.name + " cavity", triangles=cavity_faces, nodes=self.mesh.nodes
            )
            surface.compute_normals(inplace=True)  # recompute normals
            part.cavity = Cavity(surface=surface, name=part.name)
            part.cavity.compute_centroid()

            LOGGER.debug("Volume of cavity: {0} = {1}".format(part.cavity.name, part.cavity.volume))

            part.cavity.surface.write_to_stl(
                os.path.join(self.info.workdir, "-".join(part.cavity.surface.name.lower().split()))
            )

        return

    def _validate_cap_names(self):
        """Validate that caps are attached to right part."""
        for part in self.parts:
            cap_names = [c.name for c in part.caps]
            if part.name == "Left ventricle":
                expected_names = ["mitral", "aortic"]
            elif part.name == "Right ventricle":
                expected_names = ["pulmonary", "tricuspid"]
            elif part.name == "Left atrium":
                expected_names = []
            elif part.name == "Right atrium":
                expected_names = []

            valid_cap_names = True
            for cn in cap_names:
                for en in expected_names:
                    if en in cn:
                        # LOGGER.debug("Breaking...")
                        break
                    LOGGER.error(
                        "Part: {0}. Cap name is {1}, but expecting cap names"
                        "to contain one of {2}".format(part.name, cn, expected_names)
                    )
                    valid_cap_names = False
        return

    def _validate_surfaces(self):
        """Validate that none of the surfaces are empty."""
        is_valid = False
        invalid_surfaces = [s for p in self.parts for s in p.surfaces if s.n_cells == 0]
        if len(invalid_surfaces) == 0:
            is_valid = True
        else:
            for invalid_s in invalid_surfaces:
                LOGGER.error(f"Surface {invalid_s.name} is empty")
                is_valid = False

        self._sync_epicardium_with_part()

        return is_valid

    def _validate_parts(self):
        """Validate that none of the parts are empty."""
        is_valid = False
        invalid_parts = [p for p in self.parts if p.element_ids.shape[0] == 0]
        if len(invalid_parts) == 0:
            is_valid = True
        else:
            for invalid_p in invalid_parts:
                LOGGER.error(f"Part {invalid_p.name} is empty")
                is_valid = False

        return is_valid

    def _sync_epicardium_with_part(self):
        """Clean epicardial surfaces such that these use only nodes of part."""
        for part in self.parts:
            node_ids_part = np.unique(self.mesh.tetrahedrons[part.element_ids, :])
            for surface in part.surfaces:
                if "epicardium" in surface.name:
                    mask = np.all(np.isin(surface.triangles, node_ids_part), axis=1)
                    surface = SurfaceMesh(
                        name=surface.name,
                        triangles=surface.triangles[mask, :],
                        nodes=surface.points,
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

    def get_apex_node_set(
        self,
        part: ["left", "right"] = "left",
        option: ["endocardium", "epicardium", "myocardium"] = "epicardium",
        radius: float = 3,
    ) -> np.ndarray:
        """Get a node set around apex point.

        Parameters
        ----------
        part : left&quot;, &quot;right&quot;], optional
            on which part, by default "left"
        option : endocardium&quot;, &quot;epicardium&quot;, &quot;myocardium&quot;], optional
            on surface or in mesh, by default "epicardium"
        radius : float, optional
            search in radius, by default 3

        Returns
        -------
        np.ndarray
            apex node set
        """
        import scipy.spatial as spatial

        if part == "left":
            part: Part = self.left_ventricle
        elif part == "right":
            part: Part = self.right_ventricle

        point_cloud = self.mesh.points
        point_tree = spatial.cKDTree(point_cloud)
        apex_xyz = part.apex_points[1].xyz  # always from apex at epicardium
        apex_set = point_tree.query_ball_point(apex_xyz, radius)

        if option == "myocardium":
            return np.array(apex_set)
        elif option == "endocardium":
            return np.intersect1d(part.endocardium.node_ids, apex_set)
        elif option == "epicardium":
            return np.intersect1d(part.epicardium.node_ids, apex_set)


class LeftVentricle(HeartModel):
    """Model of just the left ventricle."""

    def __init__(self, info: ModelInfo = None) -> None:
        self.left_ventricle: Part = Part(name="Left ventricle", part_type="ventricle")
        """Left ventricle part."""
        # remove septum - not used in left ventricle only model
        del self.left_ventricle.septum

        self.left_ventricle.has_fiber = True
        self.left_ventricle.is_active = True

        if info:
            super().__init__(info)
        pass


class BiVentricle(HeartModel):
    """Model of the left and right ventricle."""

    def __init__(self, info: ModelInfo = None) -> None:
        self.left_ventricle: Part = Part(name="Left ventricle", part_type="ventricle")
        """Left ventricle part."""
        self.right_ventricle: Part = Part(name="Right ventricle", part_type="ventricle")
        """Right ventricle part."""
        self.septum: Part = Part(name="Septum", part_type="septum")
        """Septum."""

        self.left_ventricle.has_fiber = True
        self.left_ventricle.is_active = True
        self.right_ventricle.has_fiber = True
        self.right_ventricle.is_active = True
        self.septum.has_fiber = True
        self.septum.is_active = True

        if info:
            super().__init__(info)
        pass


class FourChamber(HeartModel):
    """Model of the left/right ventricle and left/right atrium."""

    def __init__(self, info: ModelInfo = None) -> None:
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

        self.left_ventricle.has_fiber = True
        self.left_ventricle.is_active = True
        self.right_ventricle.has_fiber = True
        self.right_ventricle.is_active = True
        self.septum.has_fiber = True
        self.septum.is_active = True

        self.left_atrium.has_fiber = False
        self.left_atrium.is_active = False
        self.right_atrium.has_fiber = False
        self.right_atrium.is_active = False

        if info:
            super().__init__(info)

        pass

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

        AtrioVentricular node is on right artrium endocardium surface and closest to septum.

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

    def compute_av_conduction(self, create_new_nodes=True, beam_length: float = 1.5) -> BeamMesh:
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
            SA_id = self.right_atrium.get_point("SA_node").node_id
        except AttributeError:
            LOGGER.info("SA node is not defined, creating with default option.")
            SA_id = self.compute_SA_node().node_id

        try:
            AV_id = self.right_atrium.get_point("AV_node").node_id
        except AttributeError:
            LOGGER.info("AV node is not defined, creating with default option.")
            AV_id = self.compute_AV_node().node_id

        path_SAN_AVN = right_atrium_endo.geodesic(SA_id, AV_id)

        beam_nodes = path_SAN_AVN.points[1:, :]
        beam_nodes = _refine_line(beam_nodes, beam_length=beam_length)

        # duplicate nodes inside the line, connect only SA node (the first) with 3D
        point_ids = np.linspace(0, len(beam_nodes) - 1, len(beam_nodes), dtype=int)
        point_ids = np.insert(point_ids, 0, SA_id)
        # build connectivity table
        edges = np.vstack((point_ids[:-1], point_ids[1:])).T

        mask = np.ones(edges.shape, dtype=bool)
        mask[0, 0] = False  # SA point at solid

        beam_net = self.add_beam_net(beam_nodes, edges, mask, pid=0, name="SAN_to_AVN")

        return beam_net

    def compute_His_conduction(self, beam_length: float = 1.5):
        """Compute His bundle conduction."""
        start_coord, bifurcation_coord = self._define_hisbundle_start_bifurcation()

        # https://www.researchgate.net/publication/353154291_Morphometric_analysis_of_the_His_bundle_atrioven
        # tricular_fascicle_in_humans_and_other_animal_species_Histological_and_immunohistochemical_study
        # # (1.06 ± 0.6 mm)

        # path start from AV point, to septum start point, then to septum end point
        av_id = None
        for beam in self.beam_network:
            if beam.name == "SAN_to_AVN":
                av_id = beam.edges[-1, -1]
                break

        if av_id is None:
            LOGGER.error(
                "Unable to find the last node of SAN_to_AVN branch, you need to define it."
            )
            exit()

        # create nodes from start to end
        new_nodes = np.array(
            [
                self.right_atrium.get_point("AV_node").xyz,
                start_coord,
                bifurcation_coord,
            ]
        )

        new_nodes = _refine_line(new_nodes, beam_length=beam_length)
        new_nodes = new_nodes[1:, :]

        point_ids = np.concatenate(
            (
                np.array([av_id]),
                +np.linspace(0, len(new_nodes) - 1, len(new_nodes), dtype=int),
            )
        )
        bifurcation_id = point_ids[-1]
        # create beam
        edges = np.vstack((point_ids[:-1], point_ids[1:])).T

        # find end on left endo
        temp_id = pv.PolyData(
            self.mesh.points[self.left_ventricle.endocardium.node_ids, :]
        ).find_closest_point(bifurcation_coord)
        his_end_left_id = self.left_ventricle.endocardium.node_ids[temp_id]
        his_end_left_coord = self.mesh.points[his_end_left_id, :]

        left_his = np.array([bifurcation_coord, his_end_left_coord])
        left_his = _refine_line(left_his, beam_length=beam_length)
        new_nodes = np.vstack((new_nodes, left_his[1:, :]))

        left_his_point_ids = np.concatenate(
            (
                np.array([bifurcation_id]),
                bifurcation_id
                + 1
                + np.linspace(0, len(left_his[1:, :]) - 1, len(left_his[1:, :]), dtype=int),
            )
        )

        edges = np.vstack(
            (edges, np.column_stack((left_his_point_ids[:-1], left_his_point_ids[1:])))
        )
        position_id_His_end_left = np.argwhere(edges == left_his_point_ids[-1])[0]
        # find end on right endo (must on the septum endo)
        temp_id = pv.PolyData(
            self.mesh.points[self.right_ventricle.septum.node_ids, :]
        ).find_closest_point(bifurcation_coord)
        his_end_right_id = self.right_ventricle.septum.node_ids[temp_id]
        his_end_right_coord = self.mesh.points[his_end_right_id, :]

        right_his = np.array([bifurcation_coord, his_end_right_coord])
        right_his = _refine_line(right_his, beam_length=beam_length)
        new_nodes = np.vstack((new_nodes, right_his[1:, :]))
        right_his_point_ids = np.concatenate(
            (
                np.array([bifurcation_id]),
                left_his_point_ids[-1]
                + 1
                + np.linspace(0, len(right_his[1:, :]) - 1, len(right_his[1:, :]), dtype=int),
            )
        )

        edges = np.vstack(
            (edges, np.column_stack((right_his_point_ids[:-1], right_his_point_ids[1:])))
        )
        position_id_His_end_right = np.argwhere(edges == right_his_point_ids[-1])[0]
        # finally
        mask = np.ones(edges.shape, dtype=bool)
        mask[0, 0] = False  # AV point at previous, not offset in creation

        beam_net = self.add_beam_net(new_nodes, edges, mask, pid=0, name="His")
        beam_net.beam_nodes_mask[0, 0] = True  # offset in writer

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
        return Point(
            xyz=his_end_left_coord,
            node_id=beam_net.edges[position_id_His_end_left[0], position_id_His_end_left[1]],
        ), Point(
            xyz=his_end_right_coord,
            node_id=beam_net.edges[position_id_His_end_right[0], position_id_His_end_right[1]],
        )

    def _define_hisbundle_start_bifurcation(self):
        """
        Define start and end points of the bundle of His.

        Start point: create a point in septum, it's closest to AV node.
        End point: create a point in septum, it's close to AV node but randomly chosen.
        """
        AV_node = self.right_atrium.get_point("AV_node")

        septum_point_ids = np.unique(np.ravel(self.mesh.tetrahedrons[self.septum.element_ids]))

        # remove nodes on surface, to make sure His bundle nodes are inside of septum
        septum_point_ids = np.setdiff1d(septum_point_ids, self.left_ventricle.endocardium.node_ids)
        septum_point_ids = np.setdiff1d(septum_point_ids, self.right_ventricle.septum.node_ids)

        septum_pointcloud = pv.PolyData(self.mesh.nodes[septum_point_ids, :])

        # Define start point: closest to artria
        pointcloud_id = septum_pointcloud.find_closest_point(AV_node.xyz)

        His_start_id = septum_point_ids[pointcloud_id]
        His_start = self.mesh.nodes[His_start_id, :]

        # Define end point:  a random close point
        n = 50
        pointcloud_id = septum_pointcloud.find_closest_point(AV_node.xyz, n=n)[n - 1]

        His_bifurcation_id = septum_point_ids[pointcloud_id]
        His_bifurcation = self.mesh.nodes[His_bifurcation_id, :]

        return His_start, His_bifurcation

    def compute_left_right_bundle(self, start_coord, start_id, side: str, beam_length: float = 1.5):
        """Bundle brunch."""
        if side == "Left":
            ventricle = self.left_ventricle
            endo_surface = self.left_ventricle.endocardium
        elif side == "Right":
            ventricle = self.right_ventricle
            face = np.hstack(
                (self.right_ventricle.endocardium.faces, self.right_ventricle.septum.faces)
            )
            endo_surface = pv.PolyData(self.mesh.points, face)

        bundle_branch = endo_surface.geodesic(
            endo_surface.find_closest_point(start_coord),
            endo_surface.find_closest_point(self.mesh.points[ventricle.apex_points[0].node_id]),
        )

        new_nodes = bundle_branch.points
        new_nodes = _refine_line(new_nodes, beam_length=beam_length)
        # exclude first and last (apex) node which belongs to purkinje beam
        new_nodes = new_nodes[1:-1, :]
        point_ids = np.linspace(0, len(new_nodes) - 1, len(new_nodes), dtype=int)
        point_ids = np.insert(point_ids, 0, start_id)
        apex = ventricle.apex_points[0].node_id
        for network in self.beam_network:
            if network.name == side + "-purkinje":
                apex = network.edges[0, 0]
        point_ids = np.append(point_ids, apex)

        edges = np.vstack((point_ids[:-1], point_ids[1:])).T

        mask = np.ones(edges.shape, dtype=bool)
        mask[0, 0] = False  # His end point of previous, not offset at creation
        mask[-1, -1] = False  # Apex point at solid

        beam_net = self.add_beam_net(new_nodes, edges, mask, pid=0, name=side + " bundle branch")
        beam_net.beam_nodes_mask[0, 0] = True
        beam_net.beam_nodes_mask[-1, -1] = True

        return beam_net

    def compute_Bachman_bundle(self):
        """Compute Bachman bundle conduction system."""
        raise NotImplementedError
        return

    def _create_isolation_part(self):
        """Create a new part to isolate between ventricles and atrial."""
        # find interface nodes between ventricles and atrial
        v_ele = np.hstack((self.left_ventricle.element_ids, self.right_ventricle.element_ids))
        a_ele = np.hstack((self.left_atrium.element_ids, self.right_atrium.element_ids))

        ventricles = self.mesh.extract_cells(v_ele)
        atrial = self.mesh.extract_cells(a_ele)

        interface_nids = np.intersect1d(
            ventricles["vtkOriginalPointIds"], atrial["vtkOriginalPointIds"]
        )

        interface_eids = np.where(np.any(np.isin(self.mesh.tetrahedrons, interface_nids), axis=1))[
            0
        ]
        # interface elements on atrial part
        interface_eids = np.intersect1d(interface_eids, a_ele)

        # remove these elements from atrial parts
        self.left_atrium.element_ids = np.setdiff1d(self.left_atrium.element_ids, interface_eids)
        self.right_atrium.element_ids = np.setdiff1d(self.right_atrium.element_ids, interface_eids)

        # find orphan elements of atrial parts and assign to isolation part
        self.mesh["cell_ids"] = np.arange(0, self.mesh.n_cells, dtype=int)
        for atrium in [self.left_atrium, self.right_atrium]:
            clean_obj = self.mesh.extract_cells(atrium.element_ids).connectivity(largest=True)
            connected_cells = clean_obj["cell_ids"]
            orphan_cells = np.setdiff1d(atrium.element_ids, connected_cells)

            # keeep largest connected part for atrial
            atrium.element_ids = connected_cells

            # get orphan cells and set to isolation part
            LOGGER.warning(f"{len(orphan_cells)} orphan cells are found and re-assigned.")
            interface_eids = np.append(interface_eids, orphan_cells)

        # create a new part
        self.add_part("Isolation atrial")
        isolation = self.get_part("Isolation atrial")
        isolation.element_ids = interface_eids
        isolation.has_fiber = True
        isolation.is_active = False


class FullHeart(FourChamber):
    """Model of both ventricles, both atria, aorta and pulmonary artery."""

    def __init__(self, info: ModelInfo = None) -> None:
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

        self.left_ventricle.has_fiber = True
        self.left_ventricle.is_active = True
        self.right_ventricle.has_fiber = True
        self.right_ventricle.is_active = True
        self.septum.has_fiber = True
        self.septum.is_active = True

        self.left_atrium.has_fiber = False
        self.left_atrium.is_active = False
        self.right_atrium.has_fiber = False
        self.right_atrium.is_active = False
        self.aorta.has_fiber = False
        self.aorta.is_active = False
        self.pulmonary_artery.has_fiber = False
        self.pulmonary_artery.is_active = False

        if info:
            super().__init__(info)

        pass


def _refine_line(nodes: np.array, beam_length: float):
    new_nodes = [nodes[0, :]]
    for beam_id in range(len(nodes) - 1):
        point_start = nodes[beam_id, :]
        point_end = nodes[beam_id + 1, :]
        points = _create_line(point_start, point_end, beam_length=beam_length)
        new_nodes = np.vstack((new_nodes, points[1:, :]))
    return new_nodes


if __name__ == "__main__":
    print("Protected")
