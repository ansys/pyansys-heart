# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Module containing classes for the various heart models."""
import copy
import json
import os

# import json
import pathlib
import pickle
import re
from typing import List, Literal, Union

import numpy as np
import pyvista as pv
from scipy.spatial.transform import Rotation as R
import yaml

from ansys.heart.core import LOG as LOGGER
from ansys.heart.preprocessor.input import _InputModel

# from ansys.heart.preprocessor.input import HEART_MODELS
import ansys.heart.preprocessor.mesh.connectivity as connectivity
import ansys.heart.preprocessor.mesh.mesher as mesher
from ansys.heart.preprocessor.mesh.objects import (
    BeamMesh,
    Cap,
    Cavity,
    Mesh,
    Part,
    PartType,
    Point,
    SurfaceMesh,
)
import ansys.heart.preprocessor.mesh.vtkmethods as vtkmethods


class ModelInfo:
    """Contains model information."""

    def __init__(
        self,
        input: Union[pathlib.Path, str, pv.PolyData, pv.UnstructuredGrid] = None,
        scalar: str = "part-id",
        part_definitions: dict = None,
        work_directory: pathlib.Path = ".",
        path_to_simulation_mesh: pathlib.Path = None,
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
        self.path_to_model: str = None
        """Path to model (in .pickle format)."""

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
            Point(name=c.name + "_center", xyz=c.centroid, node_id=c.global_centroid_id)
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

        self.aha_ids = None
        """American Heart Association ID's."""

        self.electrodes: List[Point] = []
        """Electrodes positions for ECG computing."""

        self.beam_network: List[BeamMesh] = []
        """List of beam networks in the mesh."""

        self.electrodes: List[Point] = []
        """Electrodes positions for ECG computing."""

        return

    # def __repr__(self):
    #     """Represent self as string."""
    #     return yaml.dump(self.summary(), sort_keys=False)

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

        if name in [p.name for p in self.parts]:
            LOGGER.error(f"Part {name} has existed.")
            return None

        for part in self.parts:
            try:
                part.element_ids = np.setdiff1d(part.element_ids, eids)
            except:
                LOGGER.error(f"Failed to create part {name}")
                return None

        self.add_part(name)
        new_part: Part = self.get_part(name)

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

    def mesh_volume(
        self,
        use_wrapper: bool = False,
        overwrite_existing_mesh: bool = True,
        path_to_fluent_mesh: str = None,
    ):
        """Remesh the input model and fill the volume.

        Parameters
        ----------
        use_wrapper : bool, optional
            Flag for switch to non-manifold mesher, by default False
        overwrite_existing_mesh : bool, optional
            Flag indicating whether to overwrite the existing .msh.h5 mesh, by default True
        path_to_fluent_mesh : str, optional
            Path to the generated Fluent .msh.h5 mesh, by default None

        Notes
        -----
        Note that when the input surfaces are non-manifold the wrapper tries
        to reconstruct the surface and parts. Inevitably this leads to
        reconstruction errors. Nevertheless, in many instances this approach is
        robuster than meshing from a manifold surface. Moreover, any clear interface
        between parts is potentially lost.
        """
        if not path_to_fluent_mesh:
            path_to_fluent_mesh = os.path.join(self.info.workdir, "simulation_mesh.msh.h5")

        if use_wrapper:
            LOGGER.warning("Meshing from non-manifold model not yet available.")

            fluent_mesh = mesher.mesh_from_non_manifold_input_model(
                model=self._input,
                workdir=self.info.workdir,
                mesh_size=self.info.mesh_size,
                path_to_output=path_to_fluent_mesh,
                overwrite_existing_mesh=overwrite_existing_mesh,
            )
        else:
            fluent_mesh = mesher.mesh_from_manifold_input_model(
                model=self._input,
                workdir=self.info.workdir,
                mesh_size=self.info.mesh_size,
                path_to_output=path_to_fluent_mesh,
                overwrite_existing_mesh=overwrite_existing_mesh,
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

        # remove any unused nodes
        fluent_mesh.clean()

        vtk_grid = fluent_mesh._to_vtk()

        mesh = Mesh(vtk_grid)
        mesh.cell_data["part-id"] = mesh.cell_data["cell-zone-ids"]
        mesh.cell_data["_volume-id"] = mesh.cell_data["cell-zone-ids"]
        for fluent_cell_zone in fluent_mesh.cell_zones:
            mesh._volume_id_to_name[fluent_cell_zone.id] = fluent_cell_zone.name

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
        # TODO: remove the following:
        mesh.boundaries = [
            SurfaceMesh(name=fz.name, triangles=fz.faces, nodes=mesh.nodes, id=fz.id)
            for fz in fluent_mesh.face_zones
            if not "interior" in fz.name
        ]

        for fz in fluent_mesh.face_zones:
            if "interior" not in fz.name:
                surface = SurfaceMesh(
                    name=fz.name, triangles=fz.faces, nodes=fluent_mesh.nodes, id=fz.id
                )

                mesh.add_surface(surface, int(fz.id))
                mesh._surface_id_to_name[int(fz.id)] = fz.name

        self.mesh = mesh.clean()

        return

    def _mesh_fluid_volume(self, remesh_caps: bool = True):
        """Generate a volume mesh of the cavities.

        Parameters
        ----------
        remesh_caps : bool, optional
            Flag indicating whether to remesh the caps of each cavity, by default True
        """
        # get all relevant boundaries for the fluid cavities:
        substrings_include = ["endocardium", "valve-plane", "septum"]
        substrings_include_re = "|".join(substrings_include)

        substrings_exlude = ["pulmonary-valve", "aortic-valve"]
        substrings_exlude_re = "|".join(substrings_exlude)

        boundaries_fluid = [
            b for b in self.mesh.boundaries if re.search(substrings_include_re, b.name)
        ]
        boundaries_exclude = [
            b.name for b in boundaries_fluid if re.search(substrings_exlude_re, b.name)
        ]
        boundaries_fluid = [b for b in boundaries_fluid if b.name not in boundaries_exclude]

        caps = [c._mesh for p in self.parts for c in p.caps]

        if len(boundaries_fluid) == 0:
            LOGGER.debug("Meshing of fluid cavities not possible. No fluid surfaces detected.")
            return

        if len(caps) == 0:
            LOGGER.debug("Meshing of fluid cavities not possible. No caps detected.")
            return

        LOGGER.info("Meshing fluid cavities...")

        # mesh the fluid cavities
        fluid_mesh = mesher.mesh_fluid_cavities(
            boundaries_fluid, caps, self.info.workdir, remesh_caps=remesh_caps
        )

        LOGGER.info(f"Meshed {len(fluid_mesh.cell_zones)} fluid regions...")

        # add part-ids
        cz_ids = np.sort([cz.id for cz in fluid_mesh.cell_zones])

        # ? this offset is arbitrary.
        offset = 10000
        new_ids = np.arange(cz_ids.shape[0]) + offset
        czid_to_pid = {cz_id: new_ids[ii] for ii, cz_id in enumerate(cz_ids)}

        for cz in fluid_mesh.cell_zones:
            cz.id = czid_to_pid[cz.id]

        fluid_mesh._fix_negative_cells()
        fluid_mesh_vtk = fluid_mesh._to_vtk(add_cells=True, add_faces=False)

        fluid_mesh_vtk.cell_data["part-id"] = fluid_mesh_vtk.cell_data["cell-zone-ids"]

        boundaries = [
            SurfaceMesh(name=fz.name, triangles=fz.faces, nodes=fluid_mesh.nodes, id=fz.id)
            for fz in fluid_mesh.face_zones
            if "interior" not in fz.name
        ]

        self.fluid_mesh = Mesh(fluid_mesh_vtk)
        self.fluid_mesh.boundaries = boundaries

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

    def summary(self) -> dict:
        """Get summary information of the model as a ditionary."""
        from ansys.heart.preprocessor.helpers import model_summary

        summary = model_summary(self)
        return summary

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
                LOGGER.info(
                    "\tcap: {:} | # nodes {:d}".format(cap.name, len(cap.global_node_ids_edge))
                )
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

    def dump_model(self, filename: Union[pathlib.Path, str] = None):
        """Save model to .pickle file.

        Parameters
        ----------
        filename : pathlib.Path | str, optional
            Path where the model will be saved, by default None

        Returns
        -------
        str
            Path to where the model is saved.

        Examples
        --------
        >>> model.dump_model("my_heart_model.pickle")

        """
        LOGGER.debug("Writing model to disk")

        if isinstance(filename, pathlib.Path):
            filename = str(filename)

        if not filename:
            filename = os.path.join(self.info.workdir, "heart_model.pickle")

        if os.path.isfile(filename):
            LOGGER.warning(f"Overwriting {filename}")

        with open(filename, "wb") as file:
            pickle.dump(self, file)
        self.info.dump_info()

        self.info.path_to_model = filename

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

        # fiber direction stored in cell data, but cell-to-point filter uses
        # leads to issues, where nan values in any non-volume cell may change
        # the fiber direction in the target point(s).
        mesh = self.mesh.extract_cells_by_type([pv.CellType.TETRA, pv.CellType.HEXAHEDRON])
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
            if part.part_type in [PartType.VENTRICLE]:
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

    # TODO: Should do this on the fly in dynawriter
    def _add_nodal_areas(self):
        """Compute and add nodal areas to surface nodes."""
        raise NotImplementedError("Adding nodal areas is deprecated")
        exit()

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

        self._add_surface_normals()

        self._assign_cavities_to_parts()
        self._update_cap_names()
        self._validate_cap_names()

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

        Notes
        -----
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

    def _extract_septum(self, num_layers_to_remove: int = 1) -> None:
        """Separate the septum elements from the left ventricle.

        Notes
        -----
        Uses the septum surface of the right ventricle
        """
        if not isinstance(self, (BiVentricle, FourChamber, FullHeart)):
            LOGGER.warning("Model type: {0} Not extracting septum elements".format(type(self)))
            return None

        septum_name = [s for s in self.mesh.surface_names if "septum" in s]

        if len(septum_name) > 1:
            raise ValueError("Expecting only one surface that contains string: 'septum'")
        if len(septum_name) == 0:
            raise ValueError("No boundary with name: 'septum' found")
        surface_septum = self.mesh.get_surface_by_name(septum_name[0])

        # extrude septum surface
        new_faces_septum = connectivity.remove_triangle_layers_from_trimesh(
            surface_septum.cast_to_unstructured_grid().cells_dict[pv.CellType.TRIANGLE],
            iters=num_layers_to_remove,
        )

        surface_septum = SurfaceMesh(
            nodes=surface_septum.points, triangles=new_faces_septum
        ).clean()

        septum_surface = surface_septum
        septum_surface.compute_normals()
        septum_surface = septum_surface.smooth()

        # septum_surface_vtk = vtkmethods.smooth_polydata(septum_surface_vtk)

        septum_surface_extruded = vtkmethods.extrude_polydata(septum_surface, 20)
        # septum_surface_vtk_extruded = vtkmethods.extrude_polydata(septum_surface_vtk, 20)

        # only check tetra elements
        volume_vtk = self.mesh.extract_cells_by_type(pv.CellType.TETRA)

        element_ids_septum = vtkmethods.cell_ids_inside_enclosed_surface(
            volume_vtk, septum_surface_extruded
        )
        element_ids_septum = self.mesh._global_tetrahedron_ids[element_ids_septum]

        # assign to septum
        part = next(part for part in self.parts if part.part_type == PartType.SEPTUM)
        part.element_ids = element_ids_septum
        self.mesh.cell_data["part-id"][element_ids_septum] = part.pid
        # manipulate _volume-id
        self.mesh.cell_data["_volume-id"][element_ids_septum] = part.pid
        self.mesh._volume_id_to_name[int(part.pid)] = part.name

        # remove these element ids from the left-ventricle
        part = next(part for part in self.parts if part.name == "Left ventricle")
        mask = np.isin(part.element_ids, element_ids_septum, invert=True)
        part.element_ids = part.element_ids[mask]

        return

    # TODO: since we now use a global mesh we need to
    # TODO: use the global point ids.
    def _extract_apex(self, check_edge: bool = True) -> None:
        """
        Extract the apex for both the endocardium and epicardium of each ventricle.

        Notes
        -----
        Apex defined as the point furthest from the mid-point between caps/valves

        Args:
            check_edge (bool, optional): Checks and corrects if the apex point is on surface edge.
            Defaults to True.
        """
        ventricles = [p for p in self.parts if "ventricle" in p.name]
        surface_substrings = ["endocardium", "epicardium"]
        for ventricle in ventricles:
            # get reference point (center point between two caps)
            cap_centroids = [c.centroid for c in ventricle.caps]
            ref_point = np.mean(np.array(cap_centroids), axis=0)
            for surface_substring in surface_substrings:
                surface = next(s for s in ventricle.surfaces if surface_substring in s.name)
                surface = SurfaceMesh(self.mesh.get_surface(surface.id))
                # surface = next(s for s in self.mesh._surfaces if surface_substring in s.name)
                apical_node_id = surface.node_ids[
                    np.argmax(np.linalg.norm(surface.nodes - ref_point, axis=1))
                ]
                # NOTE: This is the global apical node id when referenced in self.mesh.
                apical_node_id = surface.point_data["_global-point-ids"][apical_node_id]

                if check_edge and np.any(
                    surface.point_data["_global-point-ids"][surface.boundary_edges]
                    == apical_node_id
                ):
                    # remove one layer of elements at the boundary.
                    edgeless_surface = copy.deepcopy(surface)
                    edgeless_surface.triangles = connectivity.remove_triangle_layers_from_trimesh(
                        surface.triangles
                    )

                    apical_node_id = edgeless_surface.find_closest_point(
                        self.mesh.points[apical_node_id, :]
                    )

                    # map to global id
                    apical_node_id = edgeless_surface.point_data["_global-point-ids"][
                        apical_node_id
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
                        xyz=self.mesh.points[apical_node_id, :],
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
            # ! this is valid as long as no additional surfaces are added in self.mesh.
            # ! otherwise (global) element ids may change
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
            "{0}/{1} elements assigned to parts".format(summ, self.mesh.tetrahedrons.shape[0])
        )

        return

    # TODO refactor:
    def _assign_surfaces_to_parts(self) -> None:
        """Assign surfaces generated during remeshing to model parts."""
        for part in self.parts:
            for surface in part.surfaces:
                boundary_name = "-".join(surface.name.lower().split())
                # boundary_surface = self.mesh._get_surface_from_name(boundary_name)
                boundary_surface = SurfaceMesh(self.mesh.get_surface_by_name(boundary_name))
                boundary_surface.id = self.mesh._surface_name_to_id[boundary_name]
                if "septum" in surface.name:
                    try:
                        septum_candidates = [s for s in self.mesh.surface_names if "septum" in s]
                        if len(septum_candidates) > 1:
                            LOGGER.warning(
                                "Multiple candidate surfaces for septum found, using first one."
                            )
                        boundary_surface = self.mesh.get_surface_by_name(septum_candidates[0])
                    except:
                        boundary_surface = None
                if boundary_surface:
                    surface.triangles = boundary_surface.triangles
                    surface.nodes = boundary_surface.nodes
                    surface.id = boundary_surface.id
                else:
                    LOGGER.warning("Could not find matching surface for: {0}".format(surface.name))

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
        idoffset = 1000  # TODO need to improve id checking
        ii = 0

        for part in self.parts:
            if not hasattr(part, "endocardium"):
                continue

            cavity_faces = np.empty((0, 3), dtype=int)

            # select endocardial surfaces
            surfaces = [s for s in part.surfaces if "endocardium" in s.name]

            # NOTE, this is a loop since the right-ventricle endocardium consists
            # of both the "regular" endocardium and the septal endocardium

            # NOTE: pv.merge accepts list of surfaces.
            surface: SurfaceMesh = SurfaceMesh(pv.merge(surfaces))
            surface.name = part.name + " cavity"

            # save this cavity mesh to the centralized mesh object
            surface.id = int(np.sort(self.mesh.surface_ids)[-1] + 1)  # get unique id.

            # Generate patches that close the surface.
            patches = vtkmethods.get_patches_with_centroid(surface)

            LOGGER.debug(f"Generating {len(patches)} caps for {part.name}")

            # TODO: Add patches as a surface to self.mesh instead of modifying nodes directly.
            # TODO: Note that points come from surface, and does not contain all points in the mesh.

            # Create Cap objects with patches
            caps = []
            for patch in patches:
                ii += 1

                cap_name = f"cap_{ii}_{part.name}"

                # create cap: NOTE, mostly for compatibility. Could simplify further
                cap_mesh = SurfaceMesh(patch.clean(), name=cap_name, id=ii + idoffset)

                # Add cap to main mesh.
                self.mesh.add_surface(cap_mesh, id=cap_mesh.id)
                self.mesh._surface_id_to_name[cap_mesh.id] = cap_name

                self.mesh.clean()

                # get the cap mesh from the global mesh:
                # this ensures we can access the _global-point/cell-ids.
                cap_mesh1 = SurfaceMesh(self.mesh.get_surface_by_name(cap_name))
                cap_mesh1.id = cap_mesh.id
                cap_mesh1.name = cap_mesh.name

                cap = Cap(name=cap_name)
                cap._mesh = cap_mesh1
                part.caps.append(cap)

            # TODO: We could do this somewhere else.
            # ! Note that the element ids in cavity.surface don't have any meaning
            # ! and we can only use this for getting some meta-data such as volume
            # ! also it is not updated dynamically.
            # merge patches into cavity surface.
            surface.cell_data["_cap_id"] = 0
            surface_cavity = SurfaceMesh(pv.merge([surface] + [cap._mesh for cap in part.caps]))
            surface_cavity.name = surface.name
            surface_cavity.id = surface.id

            surface_cavity.compute_normals(inplace=True)  # recompute normals

            if not surface_cavity.is_manifold:
                LOGGER.warning("Cavity of {part.name} is not manifold.")

            # add and get from global mesh to get all point/cell data arrays.
            self.mesh.add_surface(surface_cavity, surface_cavity.id)
            self.mesh._surface_id_to_name[surface_cavity.id] = surface_cavity.name
            self.mesh = self.mesh.clean()

            surface_cavity = self.mesh.get_surface(surface_cavity.id)
            surface_cavity.id = surface.id
            surface_cavity.name = surface.name

            # NOTE: Validate if normals are pointing inward.
            part.cavity = Cavity(surface=surface_cavity, name=part.name)
            part.cavity.compute_centroid()

            LOGGER.debug("Volume of cavity: {0} = {1}".format(part.cavity.name, part.cavity.volume))

            part.cavity.surface.save(
                os.path.join(
                    self.info.workdir, "-".join(part.cavity.surface.name.lower().split()) + ".stl"
                )
            )

        return

    # TODO: use `are_connected`` to confirm which other surfaces are connected to the caps.
    def _update_cap_names(self):
        """Try to update the cap names using names of connected boundaries."""
        boundaries_to_check = [
            s for s in self.mesh._surfaces if "valve" in s.name or "inlet" in s.name
        ]
        for part in self.parts:
            for cap in part.caps:
                cap_mesh = self.mesh.get_surface_by_name(cap.name)
                for b in boundaries_to_check:
                    if vtkmethods.are_connected(cap_mesh, b):
                        for split in b.name.split("_"):
                            if "valve" in split or "inlet" in split:
                                break
                        old_cap_name = cap.name

                        cap.name = split.replace("-plane", "").replace("-inlet", "")

                        if "atrium" in part.name and (
                            "tricuspid" in cap.name or "mitral" in cap.name
                        ):
                            cap.name = cap.name + "-atrium"

                        LOGGER.debug(f"Cap {cap.name} connected to {b.name}")
                        # update name to id map:
                        cap_id = self.mesh._surface_name_to_id[old_cap_name]
                        self.mesh._surface_id_to_name[cap_id] = cap.name
                        break

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
                matches = [True for en in expected_names if en in cn]
                if len(matches) == 1:
                    break
                else:
                    LOGGER.error(
                        "Part: {0}. Cap name is {1}, but expecting cap names "
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
            self.mesh._set_global_ids()
            global_node_ids_part = self.mesh.extract_cells(part.element_ids).point_data[
                "_global-point-ids"
            ]

            # ! The only info we use from surface here is the id, and not the mesh info
            # ! we need to go back to the central mesh to obtain an updated copy of
            # ! the corresponding mesh.
            for surface in part.surfaces:
                if "epicardium" in surface.name:
                    # get the surface id.
                    surf_id = self.mesh._surface_name_to_id[surface.name.lower().replace(" ", "-")]
                    global_node_ids_surface = self.mesh.get_surface(surf_id).point_data[
                        "_global-point-ids"
                    ]
                    mask = np.isin(global_node_ids_surface, global_node_ids_part)
                    # do not use any faces that use a node not in the part.
                    mask = np.all(np.isin(surface.triangles, np.argwhere(mask).flatten()), axis=1)

                    LOGGER.debug(f"Removing {np.sum(np.invert(mask))} from {surface.name}")
                    surface.triangles = surface.triangles[mask, :]

                    # add updated mesh to global mesh.
                    # TODO: could just change the _surface-ids in the cell data directly.
                    # TODO: this basically appends this the cells of this mesh at the end of Mesh.
                    self.mesh.remove_surface(surf_id)
                    self.mesh.add_surface(surface, int(surf_id))

        return

    def compute_left_ventricle_anatomy_axis(
        self,
        mv_center: Union[None, np.ndarray] = None,
        av_center: Union[None, np.ndarray] = None,
        first_cut_short_axis=0.2,
    ):
        """Compute the long and short axes of the left ventricle.

        Parameters
        ----------
        mv_center : Union[None, np.ndarray], optional
            mitral valve center, by default None
        av_center : Union[None, np.ndarray], optional
            aortic valve center, by default None
        first_cut_short_axis : float, optional
            relative distance between mv center to apex, by default 0.2
        """
        if mv_center is None:
            try:
                mv_center = next(
                    cap.centroid for cap in self.left_ventricle.caps if cap.name == "mitral-valve"
                )
            except StopIteration:
                LOGGER.error("Cannot define mitral valve center")
                return
        if av_center is None:
            try:
                av_center = next(
                    cap.centroid for cap in self.left_ventricle.caps if cap.name == "aortic-valve"
                )
            except StopIteration:
                LOGGER.error("Cannot define mitral valve center")
                return

        # apex is defined on epicardium
        apex = next(
            ap.xyz for ap in self.left_ventricle.apex_points if ap.name == "apex epicardium"
        )

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

    # TODO: fix this.
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
            if np.dot(n - p_highest, mv_center - p_highest) > 0:
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

    # TODO: fix this.
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

    # TODO: fix this.
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
        set1 = np.unique(free_wall_center["point_ids"][dot >= 0])  # -pi
        set2 = np.unique(free_wall_center["point_ids"][dot < 0])  # pi
        set3 = np.unique(
            np.setdiff1d(septum_center["point_ids"], free_wall_center["point_ids"])
        )  # 0

        # visu
        # mesh["bc"] = np.zeros(mesh.n_points)
        # mesh["bc"][set1] = 1
        # mesh["bc"][set2] = -1
        # mesh["bc"][set3] = 2
        # mesh.set_active_scalars("bc")
        # mesh.plot()

        return set1, set2, set3

    # TODO: fix this.
    def get_apex_node_set(
        self,
        part: Literal["left", "right"] = "left",
        option: Literal["endocardium", "epicardium", "myocardium"] = "epicardium",
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
            return np.intersect1d(
                self.mesh.get_surface(part.endocardium.id).global_node_ids, apex_set
            )
        elif option == "epicardium":
            return np.intersect1d(
                self.mesh.get_surface(part.epicardium.id).global_node_ids, apex_set
            )

    # TODO: fix this.
    def _create_atrioventricular_isolation(self) -> Union[None, Part]:
        """
        Extract a layer of element to isolate between ventricles and atrium.

        Notes
        -----
        These elements are initially belong to atrium.

        Returns
        -------
        Part
            Part of isolation elements.
        """
        if not isinstance(self, FourChamber):
            LOGGER.error("This method is only for FourChamber model.")
            return

        # find interface nodes between ventricles and atrial
        v_ele = np.array([], dtype=int)
        a_ele = np.array([], dtype=int)
        #! Note that this only works since tetrahedrons are located
        #! at start of the mesh object.
        for part in self.parts:
            if part.part_type == PartType.VENTRICLE:
                v_ele = np.append(v_ele, part.element_ids)
            elif part.part_type == PartType.ATRIUM:
                a_ele = np.append(a_ele, part.element_ids)

        ventricles = self.mesh.extract_cells(v_ele)
        atrial = self.mesh.extract_cells(a_ele)

        interface_nids = np.intersect1d(
            ventricles["_global-point-ids"], atrial["_global-point-ids"]
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

            #! Central mesh object not updated. E.g. lose connection between part.element_ids and
            #! model.mesh.volume_ids/.cell_data["_volume-id"]

        # create a new part
        isolation: Part = self.create_part_by_ids(interface_eids, "Isolation atrial")
        isolation.part_type = PartType.ATRIUM
        isolation.fiber = True
        isolation.active = False

        return isolation

    # TODO: fix this.
    def create_atrial_stiff_ring(self, radius: float = 2) -> Union[None, Part]:
        """Create a part for solids close to atrial caps.

        Note
        ----
        Part will be passive and isotropic, material need to be defined

        Parameters
        ----------
        radius : foat, optional
            Influence region, by default 2

        Returns
        -------
        Union[None, Part]
            Part of atrial rings if created
        """
        if not isinstance(self, FourChamber):
            LOGGER.error("This method is only for FourChamber model.")
            return

        # get ring cells from cap node list
        ring_nodes = []
        for cap in self.left_atrium.caps:
            # update cap mesh with up-to-date mesh
            cap._mesh = self.mesh.get_surface(cap._mesh.id)
            if "mitral" not in cap.name:
                ring_nodes.extend(cap.global_node_ids_edge.tolist())
        for cap in self.right_atrium.caps:
            # update cap mesh with up-to-date mesh
            cap._mesh = self.mesh.get_surface(cap._mesh.id)
            if "tricuspid" not in cap.name:
                ring_nodes.extend(cap.global_node_ids_edge.tolist())

        ring_eles = vtkmethods.find_cells_close_to_nodes(self.mesh, ring_nodes, radius=radius)
        #! remove any non-tetrahedron elements
        ring_eles = ring_eles[np.isin(ring_eles, self.mesh._global_tetrahedron_ids)]
        # above search may create orphan elements, pick them to rings
        self.mesh["cell_ids"] = np.arange(0, self.mesh.n_cells, dtype=int)
        unselect_eles = np.setdiff1d(
            np.hstack((self.left_atrium.element_ids, self.right_atrium.element_ids)), ring_eles
        )
        largest = self.mesh.extract_cells(unselect_eles).connectivity(largest=True)
        connected_cells = largest["cell_ids"]
        orphan_cells = np.setdiff1d(unselect_eles, connected_cells)
        if len(orphan_cells) > 0:
            ring_eles = np.hstack((ring_eles, orphan_cells))

        # Create ring part
        ring: Part = self.create_part_by_ids(
            ring_eles, name="base atrial stiff rings"
        )  # TODO name must has 'base', see dynawriter.py L3120
        ring.part_type = PartType.ATRIUM
        ring.fiber = False
        ring.active = False

        return ring


class LeftVentricle(HeartModel):
    """Model of just the left ventricle."""

    def __init__(self, info: ModelInfo = None) -> None:
        self.left_ventricle: Part = Part(name="Left ventricle", part_type=PartType.VENTRICLE)
        """Left ventricle part."""
        # remove septum - not used in left ventricle only model
        del self.left_ventricle.septum

        self.left_ventricle.fiber = True
        self.left_ventricle.active = True

        if info:
            super().__init__(info)
        pass


class BiVentricle(HeartModel):
    """Model of the left and right ventricle."""

    def __init__(self, info: ModelInfo = None) -> None:
        self.left_ventricle: Part = Part(name="Left ventricle", part_type=PartType.VENTRICLE)
        """Left ventricle part."""
        self.right_ventricle: Part = Part(name="Right ventricle", part_type=PartType.VENTRICLE)
        """Right ventricle part."""
        self.septum: Part = Part(name="Septum", part_type=PartType.SEPTUM)
        """Septum."""

        self.left_ventricle.fiber = True
        self.left_ventricle.active = True
        self.right_ventricle.fiber = True
        self.right_ventricle.active = True
        self.septum.fiber = True
        self.septum.active = True

        if info:
            super().__init__(info)
        pass


class FourChamber(HeartModel):
    """Model of the left/right ventricle and left/right atrium."""

    def __init__(self, info: ModelInfo = None) -> None:
        self.left_ventricle: Part = Part(name="Left ventricle", part_type=PartType.VENTRICLE)
        """Left ventricle part."""
        self.right_ventricle: Part = Part(name="Right ventricle", part_type=PartType.VENTRICLE)
        """Right ventricle part."""
        self.septum: Part = Part(name="Septum", part_type=PartType.SEPTUM)
        """Septum."""

        self.left_atrium: Part = Part(name="Left atrium", part_type=PartType.ATRIUM)
        """Left atrium part."""
        self.right_atrium: Part = Part(name="Right atrium", part_type=PartType.ATRIUM)
        """Right atrium part."""

        self.left_ventricle.fiber = True
        self.left_ventricle.active = True
        self.right_ventricle.fiber = True
        self.right_ventricle.active = True
        self.septum.fiber = True
        self.septum.active = True

        self.left_atrium.fiber = False
        self.left_atrium.active = False
        self.right_atrium.fiber = False
        self.right_atrium.active = False

        if info:
            super().__init__(info)

        pass


class FullHeart(FourChamber):
    """Model of both ventricles, both atria, aorta and pulmonary artery."""

    def __init__(self, info: ModelInfo = None) -> None:
        self.left_ventricle: Part = Part(name="Left ventricle", part_type=PartType.VENTRICLE)
        """Left ventricle part."""
        self.right_ventricle: Part = Part(name="Right ventricle", part_type=PartType.VENTRICLE)
        """Right ventricle part."""
        self.septum: Part = Part(name="Septum", part_type=PartType.SEPTUM)
        """Septum."""
        self.left_atrium: Part = Part(name="Left atrium", part_type=PartType.ATRIUM)
        """Left atrium part."""
        self.right_atrium: Part = Part(name="Right atrium", part_type=PartType.ATRIUM)
        """Right atrium part."""

        self.aorta: Part = Part(name="Aorta", part_type=PartType.ARTERY)
        """Aorta part."""
        self.pulmonary_artery: Part = Part(name="Pulmonary artery", part_type=PartType.ARTERY)
        """Pulmonary artery part."""

        self.left_ventricle.fiber = True
        self.left_ventricle.active = True
        self.right_ventricle.fiber = True
        self.right_ventricle.active = True
        self.septum.fiber = True
        self.septum.active = True

        self.left_atrium.fiber = False
        self.left_atrium.active = False
        self.right_atrium.fiber = False
        self.right_atrium.active = False
        self.aorta.fiber = False
        self.aorta.active = False
        self.pulmonary_artery.fiber = False
        self.pulmonary_artery.active = False

        if info:
            super().__init__(info)

        pass


if __name__ == "__main__":
    print("Protected")
