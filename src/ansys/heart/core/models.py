# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
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
import re
from typing import List, Literal, Union

import numpy as np
import pyvista as pv
import yaml

from ansys.heart.core import LOG as LOGGER
from ansys.heart.core.objects import (
    BeamMesh,
    Cap,
    CapType,
    Cavity,
    Mesh,
    Part,
    PartType,
    Point,
    SurfaceMesh,
)
import ansys.heart.core.utils.connectivity as connectivity
import ansys.heart.core.utils.vtk_utils as vtk_utils
from ansys.heart.preprocessor.input import _InputModel
import ansys.heart.preprocessor.mesher as mesher
from ansys.heart.simulator.settings.material.ep_material import EPMaterial
from ansys.heart.simulator.settings.material.material import (
    ISO,
    Mat295,
    MechanicalMaterialModel,
)


def _get_axis_from_field_data(
    mesh: Mesh | pv.UnstructuredGrid, axis_name: Literal["l4cv_axis", "l2cv_axis", "short_axis"]
) -> dict:
    """Get the axis from mesh field data."""
    try:
        return {
            "center": mesh.field_data[axis_name][0],
            "normal": mesh.field_data[axis_name][1],
        }
    except KeyError:
        LOGGER.info(f"Failed to retrieve {axis_name} from mesh field data")
        return None


def _set_field_data_from_axis(
    mesh: Mesh | pv.UnstructuredGrid,
    axis: dict,
    axis_name: Literal["l4cv_axis", "l2cv_axis", "short_axis"],
):
    """Store the axis in mesh field data."""
    if "center" and "normal" not in axis.keys():
        LOGGER.info("Failed to store axis.")
        return None
    data = np.array([value for value in axis.values()])
    if data.shape != (2, 3):
        LOGGER.info("Data has wrong shape, expecting (2,3) shaped data.")
        return None
    mesh.field_data[axis_name] = data
    return mesh


def _set_workdir(workdir: pathlib.Path | str = None) -> str:
    """Set the root working directory.

    Parameters
    ----------
    workdir : pathlib.Path | str, optional
        Path to desired working directory, by default None

    Returns
    -------
    str
        Path to working directory.
    """
    if workdir is None:
        workdir = os.getcwd()

    workdir1 = os.getenv("PYANSYS_HEART_WORKDIR", workdir)
    if workdir1 != workdir:
        LOGGER.info(f"Working directory set to {workdir1}")

    if not os.path.isdir(workdir1):
        LOGGER.info(f"Creating {workdir1}")
        os.makedirs(workdir1, exist_ok=True)
    else:
        LOGGER.warning(f"Working directory {workdir1} already exists.")

    return workdir1


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
    def l4cv_axis(self) -> dict:
        """l4cv axis."""
        return _get_axis_from_field_data(self.mesh, "l4cv_axis")

    @property
    def l2cv_axis(self) -> dict:
        """l2cv axis."""
        return _get_axis_from_field_data(self.mesh, "l2cv_axis")

    @property
    def short_axis(self) -> dict:
        """l2cv axis."""
        return _get_axis_from_field_data(self.mesh, "short_axis")

    @l4cv_axis.setter
    def l4cv_axis(self, axis: dict):
        """Set short axis."""
        _set_field_data_from_axis(self.mesh, axis, "l4cv_axis")

    @l2cv_axis.setter
    def l2cv_axis(self, axis: dict):
        """Set short axis."""
        _set_field_data_from_axis(self.mesh, axis, "l2cv_axis")

    @short_axis.setter
    def short_axis(self, axis: dict):
        """Set short axis."""
        _set_field_data_from_axis(self.mesh, axis, "short_axis")

    @property
    def cap_centroids(self):
        """Return list of cap centroids."""
        return [
            Point(name=c.name + "_center", xyz=c.centroid, node_id=c.global_centroid_id)
            for p in self.parts
            for c in p.caps
        ]

    def __init__(self, working_directory: pathlib.Path | str = None) -> None:
        """Initialize the HeartModel.

        Parameters
        ----------
        working_directory : pathlib.Path | str, optional
            Path to desired working directory, by default None

        Notes
        -----
        Note that if no working directory is specified it will default to the current
        working directory.
        """
        self.workdir = _set_workdir(working_directory)
        """Working directory."""

        self.mesh = Mesh()
        """Computational mesh."""

        self.fluid_mesh = Mesh()
        """Generated fluid mesh."""

        #! TODO: non-functional flag. Remove or replace.
        self._add_blood_pool: bool = False
        """Flag indicating whether to add a blood pool mesh (Experimental)."""

        self._input: _InputModel = None
        """Input model."""

        self._add_subparts()
        """Add any subparts."""

        self._set_part_ids()
        """Set incremental part ids."""

        self.electrodes: List[Point] = []
        """Electrodes positions for ECG computing."""

        self.beam_network: List[BeamMesh] = []
        """List of beam networks in the mesh."""

        self.electrodes: List[Point] = []
        """Electrodes positions for ECG computing."""

        self._part_info = {}
        """Information about all the parts in the model."""

        self._short_axis: dict = None
        """Short axis."""
        self._l2cv_axis: dict = None
        """l2cv axis."""
        self._l4cv_axis: dict = None
        """l4cv axis."""

        return

    def __str__(self):
        """Represent self as string."""
        return yaml.dump(self.summary(), sort_keys=False)

    # TODO: There is overlap with the input module.
    def _get_parts_info(self):
        """Get the id to model map that allows reconstructing the model from a mesh object."""
        for part in self.parts:
            self._part_info.update(part._get_info())
        return self._part_info

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
            LOGGER.error(f"Failed to create {name}. Element list is empty")
            return None

        if name in [p.name for p in self.parts]:
            LOGGER.error(f"Failed to create {name}. Name already exists.")
            return None

        for part in self.parts:
            try:
                part.element_ids = np.setdiff1d(part.element_ids, eids)
            except ValueError:
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
        # TODO: origin of Purkinje and the apex of myiocardium seems logical, but it has more chance
        # TODO: the EP wave will not be triggered.
        # TODO: so I remove it, it means: end of bundle branch connect to apex, origin of Purkinje
        # TODO: is another point on the same location.

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
        edges[mask] += len(self.mesh.points) + len(BeamMesh.all_beam_nodes)

        if len(BeamMesh.all_beam_nodes) == 0:
            BeamMesh.all_beam_nodes = beam_nodes
        else:
            BeamMesh.all_beam_nodes = np.vstack((BeamMesh.all_beam_nodes, beam_nodes))

        # nodes is just for pyvista plot, edges used in writer will be offset
        # TODO: only save necessary nodes, cells, and with a 'global id' array
        beam_net = BeamMesh(
            nodes=np.vstack((self.mesh.points, BeamMesh.all_beam_nodes)),
            edges=edges,
            beam_nodes_mask=mask,
        )
        beam_net.pid = pid
        beam_net.name = name

        #! class variable BeamMesh.all_beam_nodes is not saved in pickle
        #! only the last created beam_network holds all previously created
        #! nodes.
        beam_net._all_beam_nodes = np.copy(BeamMesh.all_beam_nodes)

        self.beam_network.append(beam_net)

        # sync across beam networks
        for beams in self.beam_network:
            beams._all_beam_nodes = np.copy(BeamMesh.all_beam_nodes)

        return beam_net

    def load_input(self, input_vtp: pv.PolyData, part_definitions: dict, scalar: str):
        """Load an input model.

        Parameters
        ----------
        input_vtp : pv.PolyData
            The input surface mesh, represented by a VTK PolyData object.
        part_definitions : dict
            Part definitions of the input model. Each part is enclosed by N number of boundaries.
        scalar : str
            Scalar used to identify boundaries.
        """
        self._input = _InputModel(
            input=input_vtp,
            part_definitions=part_definitions,
            scalar=scalar,
        )
        if self._input is None:
            LOGGER.error("Failed to initialize input model. Please check input arguments.")
            exit()
        return

    # TODO: add working example in docstring, e.g. by using an existing model.
    def mesh_volume(
        self,
        use_wrapper: bool = False,
        overwrite_existing_mesh: bool = True,
        global_mesh_size: float = 1.5,
        path_to_fluent_mesh: str = None,
        mesh_size_per_part: dict = None,
        _global_wrap_size: float = 1.5,
        _wrap_size_per_part: dict = None,
    ) -> Mesh:
        """Remesh the input model and fill the volume.

        Parameters
        ----------
        use_wrapper : bool, optional
            Flag for switch to non-manifold mesher, by default False
        overwrite_existing_mesh : bool, optional
            Flag indicating whether to overwrite the existing .msh.h5 mesh, by default True
        global_mesh_size : float, optional
            Global mesh size used for the generated mesh, by default 1.5
        path_to_fluent_mesh : str, optional
            Path to the generated Fluent .msh.h5 mesh, by default None
        mesh_size_per_part : dict, optional
            Dictionary specifying the target mesh size for each part, by default None.
        _global_wrap_size : float, optional
            Global size used for setting up the size-field for the shrink-wrap algorithm,
            by default None
        _wrap_size_per_part : dict, optional
            Per part size used for setting up the size-field for the shrink-wrap algorithm,
            by default None

        Examples
        --------
        >>> from ansys.heart.core.models import HeartModel
        >>> model = HeartModel()
        >>> model.load_input(geom, part_definitions, scalar)
        >>> # mesh the volume with a global size of 1.5 and size of 1 for the left ventricle.
        >>> model.mesh_volume(
        ...     use_wrapper=True,
        ...     global_mesh_size=1.5,
        ...     path_to_fluent_mesh="simulation-mesh.msh.h5",
        ...     mesh_size_per_part={"Left ventricle": 1},
        ... )

        Notes
        -----
        When the input surfaces are non-manifold the wrapper tries
        to reconstruct the surface and parts. Inevitably this leads to
        reconstruction errors. Nevertheless, in many instances this approach is
        robuster than meshing from a manifold surface. Moreover, any clear interface
        between parts is potentially lost.
        When mesh_size_per_part is incomplete, remaining part sizes default to the
        global mesh size. This is an experimental setting. Any wrap sizes given
        as input argument are ignored when the wrapper is not used.
        """
        if not path_to_fluent_mesh:
            path_to_fluent_mesh = os.path.join(self.workdir, "simulation_mesh.msh.h5")

        if use_wrapper:
            self.mesh = mesher.mesh_from_non_manifold_input_model(
                model=self._input,
                workdir=self.workdir,
                global_mesh_size=global_mesh_size,
                path_to_output=path_to_fluent_mesh,
                overwrite_existing_mesh=overwrite_existing_mesh,
                mesh_size_per_part=mesh_size_per_part,
                _global_wrap_size=_global_wrap_size,
                _wrap_size_per_part=_wrap_size_per_part,
            )
        else:
            LOGGER.warning("Meshing from manifold model is experimental.")
            self.mesh = mesher.mesh_from_manifold_input_model(
                model=self._input,
                workdir=self.workdir,
                mesh_size=global_mesh_size,
                path_to_output=path_to_fluent_mesh,
                overwrite_existing_mesh=overwrite_existing_mesh,
            )

        filename = os.path.join(self.workdir, "volume-mesh-post-meshing.vtu")
        self.mesh.save(filename)

        return self.mesh

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
            b for b in self.mesh._surfaces if re.search(substrings_include_re, b.name)
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
        fluid_mesh = mesher._mesh_fluid_cavities(
            boundaries_fluid, caps, self.workdir, remesh_caps=remesh_caps
        )

        LOGGER.info(f"Meshed {len(fluid_mesh.cell_zones)} fluid regions...")

        # add part-ids
        cz_ids = np.sort([cz.id for cz in fluid_mesh.cell_zones])

        # TODO: this offset is arbitrary.
        offset = 10000
        new_ids = np.arange(cz_ids.shape[0]) + offset
        czid_to_pid = {cz_id: new_ids[ii] for ii, cz_id in enumerate(cz_ids)}

        for cz in fluid_mesh.cell_zones:
            cz.id = czid_to_pid[cz.id]

        fluid_mesh._fix_negative_cells()
        fluid_mesh_vtk = fluid_mesh._to_vtk(add_cells=True, add_faces=False)

        fluid_mesh_vtk.cell_data["_volume-id"] = fluid_mesh_vtk.cell_data["cell-zone-ids"]

        boundaries = [
            SurfaceMesh(name=fz.name, triangles=fz.faces, nodes=fluid_mesh.nodes, id=fz.id)
            for fz in fluid_mesh.face_zones
            if "interior" not in fz.name
        ]

        self.fluid_mesh = Mesh(fluid_mesh_vtk)
        for boundary in boundaries:
            self.fluid_mesh.add_surface(boundary, boundary.id, boundary.name)

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
        """Get summary information of the model as a dictionary."""
        from ansys.heart.core.utils.misc import model_summary

        summary = model_summary(self)
        return summary

    def plot_mesh(self, show_edges: bool = True, color_by: str = "_volume-id"):
        """Plot the volume mesh of the heart model.

        Parameters
        ----------
        show_edges : bool, optional
            Whether to plot the edges, by default True
        color_by : str, optional
            Color by cell/point data, by default "_volume-id"

        Examples
        --------
        >>> import ansys.heart.preprocessor.models as models
        >>> model = models.HeartModel.load_model("heart_model.pickle")
        >>> model.plot_mesh(show_edges=True)
        """
        plotter = pv.Plotter()
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
        mesh = self.mesh

        plotter = pv.Plotter()
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
        plotter = pv.Plotter()

        # fiber direction is stored in cell data, but the cell-to-point filter
        # leads to issues, where nan values in any non-volume cell may change
        # the fiber direction in the target point(s).
        mesh = self.mesh.extract_cells_by_type([pv.CellType.TETRA, pv.CellType.HEXAHEDRON])
        mesh = mesh.ctp()
        streamlines = mesh.streamlines(vectors="fiber", source_radius=75, n_points=n_seed_points)
        if streamlines.n_cells == 0:
            LOGGER.error(
                "Failed to generate streanlines with radius {source_radius} and {n_seed_points}"
            )
            return None
        tubes = streamlines.tube()
        plotter.add_mesh(mesh, opacity=0.5, color="white")
        plotter.add_mesh(tubes, color="white")
        plotter.show()
        return plotter

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
            import matplotlib as matplotlib
        except ImportError:
            LOGGER.warning("matplotlib not found. Install matplotlib with: pip install matplotlib")
            return

        surfaces_to_plot = [s for p in self.parts for s in p.surfaces]
        valves = [b for b in self.mesh._surfaces if "valve" in b.name or "border" in b.name]
        surfaces_to_plot = surfaces_to_plot + valves

        color_map = matplotlib.pyplot.get_cmap("tab20", len(surfaces_to_plot))

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
            plotter = pv.Plotter()
            plotter.add_mesh(self.mesh, color="w", opacity=0.3)
            for beams in self.beam_network:
                plotter.add_mesh(beams, color="r", line_width=2)
            plotter.show()
        except Exception:
            LOGGER.warning("Failed to plot mesh.")
        return

    def save_model(self, filename: str):
        """Save the model and necessary info to reconstruct.

        Parameters
        ----------
        filename : str
            Path to the model

        Notes
        -----
        The mesh of the heart model will be saved as .vtu file, and
        an additional partinfo.json file will be written to reconstruct
        the heart model from the VTU file.

        Examples
        --------
        >>> model.save_model("my-heart-model.vtu")

        """
        extension = pathlib.Path(filename).suffix
        if extension != "":
            mesh_path = filename.replace(extension, ".vtu")
            map_path = filename.replace(extension, ".partinfo.json")
        else:
            mesh_path = filename + ".vtu"
            map_path = filename + ".partinfo.json"

        self.mesh.save(mesh_path)

        with open(map_path, "w") as f:
            json.dump(self._get_parts_info(), f, indent=4)

        return

    # TODO: could consider having this as a static method.
    # TODO: Right now this only reconstructs the surfaces and parts that
    # TODO: are defined in the HeartModel classes:
    # TODO: LeftVentricle, BiVentricle, FourChamber and FullHeart.
    # TODO: Should consider to also reconstruct the parts that are not explicitly
    # TODO: defined in the class.
    def load_model_from_mesh(self, filename_mesh: str, filename_part_info: str):
        """Load model from an existing VTU file and part info dictionary.

        Parameters
        ----------
        filename_mesh : str
            Path to the VTU file containing the mesh.
        filename_part_info : str
            Path to the JSON file that contains the part info to reconstruct the model.

        Examples
        --------
        >>> from ansys.heart.preprocessor.models import FullHeart
        >>> model: FullHeart = FullHeart()
        >>> model.load_model_from_mesh("mesh.vtu", "mesh.partinfo.json")

        """
        # try to load the mesh.
        self.mesh.load_mesh(filename_mesh)

        # open part info
        with open(filename_part_info, "r") as f:
            self._part_info = json.load(f)
            part_info = self._part_info

        # try to reconstruct parts from part info
        # TODO: @mhoeijm alternatively we can use parts defined in part_info,
        # TODO: but we lose autocomplete. e.g.
        # TODO: use keys in part info: for part_name in part_info.keys():
        # TODO: init part by:
        # TODO: part = Part(part_1.name, PartType(part_info[part_1.name]["part-type"]))
        for part_1 in self.parts:
            try:
                list(part_info.keys()).index(part_1.name)
            except ValueError:
                LOGGER.warning(f"{part_1.name} not in part info")
                continue

            #! try to add surfaces to part by using the pre-defined surfaces
            #! Should part-info define the entire heart model and part attributes?
            for surface in part_1.surfaces:
                surface1 = self.mesh.get_surface_by_name(surface.name)
                if not surface1:
                    continue
                super(SurfaceMesh, surface).__init__(surface1)
                surface.id = surface1.id
                surface.name = surface1.name

            part_1.pid = part_info[part_1.name]["part-id"]

            try:
                part_1.element_ids = np.argwhere(
                    np.isin(self.mesh.cell_data["_volume-id"], part_1.pid)
                ).flatten()
            except Exception:
                LOGGER.warning(f"Failed to set element ids for {part_1.name}")
                pass

            # try to initialize cavity object.
            if part_info[part_1.name]["cavity"] != {}:
                cavity_name = list(part_info[part_1.name]["cavity"].keys())[0]
                cavity_id = list(part_info[part_1.name]["cavity"].values())[0]
                part_1.cavity = Cavity(surface=self.mesh.get_surface(cavity_id), name=cavity_name)

            if part_info[part_1.name]["caps"] != {}:
                for cap_name, cap_id in part_info[part_1.name]["caps"].items():
                    #! note that we sasume cap name equals cap type here.
                    cap = Cap(cap_name, cap_type=CapType(cap_name))
                    cap._mesh = self.mesh.get_surface(cap_id)
                    part_1.caps.append(cap)

            # TODO: add non-standard part by setattr(self, part_name_n, part)

        # NOTE: #? Wrap in try-block?
        # NOTE: #? add validation method to make sure all essential components are present?
        try:
            self._extract_apex()
        except Exception:
            LOGGER.warning("Failed to extract apex. Consider setting apex manually.")

        if any(v is None for v in [self.short_axis, self.l4cv_axis, self.l2cv_axis]):
            LOGGER.warning("Heart not defined in the VTU file.")
            try:
                LOGGER.warning("Computing heart axis...")
                self._define_anatomy_axis()
            except Exception:
                LOGGER.error(
                    "Failed to extract heart axis. Consider computing and setting them manually."
                )
        else:
            LOGGER.info("Heart axis defined in the VTU file is reused...")

        return

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
        return

    def _get_used_element_ids(self) -> np.ndarray:
        """Return array of used element ids."""
        element_ids = np.empty(0, dtype=int)
        for part in self.parts:
            element_ids = np.append(element_ids, part.element_ids)

        return element_ids

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

        self._assign_cavities_to_parts()
        self._update_cap_types()
        self._validate_cap_names()

        self._extract_apex()
        self._define_anatomy_axis()

        if "fiber" not in self.mesh.array_names:
            LOGGER.debug("Adding placeholder for fiber direction.")
            fiber = np.tile([[0.0, 0.0, 1.0]], (self.mesh.n_cells, 1))
            self.mesh.cell_data["fiber"] = fiber

        if "sheet" not in self.mesh.array_names:
            LOGGER.debug("Adding placeholder for sheet direction.")
            sheet = np.tile([[0.0, 1.0, 1.0]], (self.mesh.n_cells, 1))
            self.mesh.cell_data["sheet"] = sheet

        self._get_parts_info()

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
                p.pid = max([pid for pid in self.part_ids if pid is not None]) + 1

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

        septum_name = [
            s
            for s in self.mesh.surface_names
            if "right" in s.lower() and "ventricle" in s.lower() and "septum" in s.lower()
        ]

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

        septum_surface_extruded = vtk_utils.extrude_polydata(septum_surface, 20)

        # only check tetra elements
        volume_vtk = self.mesh.extract_cells_by_type(pv.CellType.TETRA)

        element_ids_septum = vtk_utils.cell_ids_inside_enclosed_surface(
            volume_vtk, septum_surface_extruded
        )
        element_ids_septum = self.mesh._global_tetrahedron_ids[element_ids_septum]

        # assign to septum
        part = next(part for part in self.parts if part.part_type == PartType.SEPTUM)
        part.element_ids = element_ids_septum
        # manipulate _volume-id
        self.mesh.cell_data["_volume-id"][element_ids_septum] = part.pid
        self.mesh._volume_id_to_name[int(part.pid)] = part.name

        # remove these element ids from the left-ventricle
        part = next(part for part in self.parts if part.name == "Left ventricle")
        mask = np.isin(part.element_ids, element_ids_septum, invert=True)
        part.element_ids = part.element_ids[mask]

        return

    def _extract_apex(self, check_edge: bool = True) -> None:
        """Extract the apex of the ventricles.

        Notes
        -----
        Apex is the defined as the point furthest from the mid-point between cap/valves.

        Parameters
        ----------
        check_edge : bool, optional
            Checks and corrects if the apical point is on the edge of a surface, by default True
        """
        ventricles = [p for p in self.parts if "ventricle" in p.name]
        surface_substrings = ["endocardium", "epicardium"]
        for ventricle in ventricles:
            # get reference point (center point between two caps).
            cap_centroids = [c.centroid for c in ventricle.caps]
            ref_point = np.mean(np.array(cap_centroids), axis=0)
            for surface_substring in surface_substrings:
                surface_id = next((s.id for s in ventricle.surfaces if surface_substring in s.name))
                surface = self.mesh.get_surface(surface_id)

                apical_node_id = surface.node_ids_triangles[
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
                        f"Initial apical point is on edge of {surface.name}, the next closest point is used"  # noqa: E501
                    )

                # assign apex point
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
        # get element ids of each part.
        used_element_ids = self._get_used_element_ids()
        for part in self.parts:
            if len(part.element_ids) > 0:
                LOGGER.warning(
                    "Part {0} seems to already have elements assigned: skipping".format(part.name)
                )
                continue
            # ! this is valid as long as no additional surfaces are added in self.mesh.
            # ! otherwise (global) element ids may change
            element_ids = np.where(np.isin(self.mesh.cell_data["_volume-id"], part.pid))[0]
            element_ids = element_ids[np.isin(element_ids, used_element_ids, invert=True)]
            part.element_ids = element_ids

        summ = 0
        for part in self.parts:
            LOGGER.info("Num elements in {0}: {1}".format(part.name, part.element_ids.shape[0]))
            summ = summ + part.element_ids.shape[0]
        LOGGER.info("Total num elements: {}".format(summ))

        LOGGER.info(
            "{0}/{1} elements assigned to parts".format(summ, self.mesh.tetrahedrons.shape[0])
        )

        return

    # TODO: refactor:
    def _assign_surfaces_to_parts(self) -> None:
        """Assign surfaces generated during remeshing to model parts."""
        for part in self.parts:
            for surface in part.surfaces:
                boundary_name = "-".join(surface.name.lower().split())
                boundary_surface = self.mesh.get_surface_by_name(boundary_name)

                if "septum" in surface.name.lower() and "right ventricle" in surface.name.lower():
                    try:
                        septum_candidates = [s for s in self.mesh.surface_names if "septum" in s]
                        if len(septum_candidates) > 1:
                            LOGGER.warning(
                                "Multiple candidate surfaces for septum found, using first one."
                            )
                        boundary_surface = self.mesh.get_surface_by_name(septum_candidates[0])
                    except Exception:
                        boundary_surface = None

                if boundary_surface:
                    #! change boundary name in self.mesh to align with heart model: note that
                    #! we may want to do this in another place.
                    self.mesh._surface_id_to_name[boundary_surface.id] = surface.name
                    super(SurfaceMesh, surface).__init__(boundary_surface)
                    surface.id = boundary_surface.id

                else:
                    LOGGER.info("Could not find matching surface for: {0}".format(surface.name))

        return

    def _assign_cavities_to_parts(self) -> None:
        """Create cavities based on endocardium surfaces and cap definitions."""
        # construct cavities with endocardium and caps
        starting_id = int(np.sort(self.mesh.surface_ids)[-1] + 1)
        starting_cap_id = int((starting_id // 1000 + 2) * 1000)
        iii = 0

        for ii, part in enumerate(self.parts):
            if not hasattr(part, "endocardium"):
                continue

            # select endocardial surfaces
            # NOTE, this is a loop since the right-ventricle endocardium consists
            # of both the "regular" endocardium and the septal endocardium.
            surfaces = [
                self.mesh.get_surface(s.id)
                for s in part.surfaces
                if "endocardium" in s.name and s.n_cells > 0
            ]
            if len(surfaces) == 0:
                LOGGER.warning(f"Skipping part {part.name}: only empty surfaces present.")
                continue

            surface: SurfaceMesh = SurfaceMesh(pv.merge(surfaces))
            surface.name = part.name + " cavity"

            # save this cavity mesh to the centralized mesh object
            surface.id = starting_id + ii  # get unique id.

            # Generate patches that close the surface.
            patches = vtk_utils.get_patches_with_centroid(surface)

            LOGGER.debug(f"Generating {len(patches)} caps for {part.name}")

            # TODO: Note that points come from surface, and does not contain all points in the mesh.
            # Create Cap objects with patches.
            for patch in patches:
                iii += 1

                cap_name = f"cap_{iii}_{part.name}"

                # create cap: NOTE, mostly for compatibility. Could simplify further
                cap_mesh = SurfaceMesh(patch.clean(), name=cap_name, id=starting_cap_id + iii)

                # Add cap to main mesh.
                self.mesh.add_surface(cap_mesh, id=cap_mesh.id, name=cap_name)

                self.mesh.clean()

                # get the cap mesh from the global mesh:
                # this ensures we can access the _global-point/cell-ids.
                cap_mesh1 = self.mesh.get_surface_by_name(cap_name)

                cap = Cap(name=cap_name)
                cap._mesh = cap_mesh1
                part.caps.append(cap)

            # TODO: We could do this somewhere else.
            # ! Note that the element/cell ids in cavity.surface don't have any meaning
            # ! and we can only use this for getting some meta-data such as volume
            # ! also it is not updated dynamically.
            # merge patches into cavity surface.
            surface.cell_data["_cap_id"] = 0
            surface_cavity = SurfaceMesh(pv.merge([surface] + [cap._mesh for cap in part.caps]))
            surface_cavity.name = surface.name
            surface_cavity.id = surface.id

            #! Force normals of cavity surface to point inward.
            surface_cavity.force_normals_inwards()

            # add and get from global mesh to get all point/cell data arrays.
            self.mesh.add_surface(surface_cavity, id=surface_cavity.id, name=surface_cavity.name)
            self.mesh = self.mesh.clean()

            surface_cavity = self.mesh.get_surface(surface_cavity.id)

            part.cavity = Cavity(surface=surface_cavity, name=surface_cavity.name)
            part.cavity.compute_centroid()

            LOGGER.debug("Volume of cavity: {0} = {1}".format(part.cavity.name, part.cavity.volume))

            part.cavity.surface.save(
                os.path.join(
                    self.workdir, "-".join(part.cavity.surface.name.lower().split()) + ".stl"
                )
            )

        return

    def _update_cap_types(self):
        """Try to update the cap types using names of connected boundaries."""
        boundaries_to_check = [
            s for s in self.mesh._surfaces if "valve" in s.name or "inlet" in s.name
        ]
        for part in self.parts:
            for cap in part.caps:
                cap_mesh = self.mesh.get_surface_by_name(cap.name)
                for b in boundaries_to_check:
                    if vtk_utils.are_connected(cap_mesh, b):
                        for split in b.name.split("_"):
                            if "valve" in split or "inlet" in split:
                                break

                        cap_name = split.replace("-plane", "").replace("-inlet", "")
                        cap.type = CapType(cap_name)

                        if "atrium" in part.name and (
                            cap.type in [CapType.TRICUSPID_VALVE, CapType.MITRAL_VALVE]
                        ):
                            cap_name = cap_name + "-atrium"
                            cap.type = CapType(cap.type.value + "-atrium")

                        cap.name = cap_name

                        LOGGER.debug(f"Cap {cap.type.value} connected to {b.name}")
                        # update name to id map:
                        self.mesh._surface_id_to_name[cap_mesh.id] = cap.type.value
                        break

        return

    def _validate_cap_names(self):
        """Validate that caps are attached to right part."""
        for part in self.parts:
            cap_types = [c.type for c in part.caps]
            if part.name == "Left ventricle":
                expected_cap_types = [
                    CapType.MITRAL_VALVE,
                    CapType.AORTIC_VALVE,
                    CapType.COMBINED_MITRAL_AORTIC_VALVE,
                ]
            elif part.name == "Right ventricle":
                expected_cap_types = [CapType.PULMONARY_VALVE, CapType.TRICUSPID_VALVE]
            elif part.name == "Left atrium":
                expected_cap_types = [
                    CapType.LEFT_ATRIUM_APPENDAGE,
                    CapType.LEFT_INFERIOR_PULMONARY_VEIN,
                    CapType.LEFT_SUPERIOR_PULMONARY_VEIN,
                    CapType.RIGHT_INFERIOR_PULMONARY_VEIN,
                    CapType.RIGHT_SUPERIOR_PULMONARY_VEIN,
                    CapType.MITRAL_VALVE_ATRIUM,
                ]
            elif part.name == "Right atrium":
                expected_cap_types = [
                    CapType.TRICUSPID_VALVE_ATRIUM,
                    CapType.SUPERIOR_VENA_CAVA,
                    CapType.INFERIOR_VENA_CAVA,
                ]

            unexpected_captypes = [ctype for ctype in cap_types if ctype not in expected_cap_types]
            if len(unexpected_captypes) > 0:
                LOGGER.error(
                    "Part: {0}. Cap types {1} not in expected cap types:{2}".format(
                        part.name, unexpected_captypes, expected_cap_types
                    )
                )

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
                    surf_id = self.mesh._surface_name_to_id[surface.name]
                    global_node_ids_surface = self.mesh.get_surface(surf_id).point_data[
                        "_global-point-ids"
                    ]
                    mask = np.isin(global_node_ids_surface, global_node_ids_part)
                    # do not use any faces that use a node not in the part.
                    mask = np.all(np.isin(surface.triangles, np.argwhere(mask).flatten()), axis=1)

                    LOGGER.debug(f"Removing {np.sum(np.invert(mask))} faces from {surface.name}")
                    surface.triangles = surface.triangles[mask, :]

                    # add updated mesh to global mesh.
                    # TODO: could just change the _surface-ids in the cell data directly.
                    # TODO: this basically appends the cells of surface at the end of Mesh.
                    self.mesh.remove_surface(surf_id)
                    self.mesh.add_surface(surface, int(surf_id), name=surface.name)

        return

    def _define_anatomy_axis(self):
        """Define long and short axes from left ventricle landmarks."""
        from ansys.heart.core.utils.landmark_utils import compute_anatomy_axis

        mv_center = next(
            cap.centroid for cap in self.left_ventricle.caps if cap.type == CapType.MITRAL_VALVE
        )

        av_center = next(
            cap.centroid for cap in self.left_ventricle.caps if cap.type == CapType.AORTIC_VALVE
        )

        apex = next(
            ap.xyz for ap in self.left_ventricle.apex_points if ap.name == "apex epicardium"
        )

        l4cv, l2cv, short = compute_anatomy_axis(
            mv_center, av_center, apex, first_cut_short_axis=0.2
        )

        self.l4cv_axis = l4cv
        self.l2cv_axis = l2cv
        self.short_axis = short

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
                self.mesh.get_surface(part.endocardium.id).global_node_ids_triangles, apex_set
            )
        elif option == "epicardium":
            return np.intersect1d(
                self.mesh.get_surface(part.epicardium.id).global_node_ids_triangles, apex_set
            )

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
        # TODO: move this method to FourChamber class.
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
            LOGGER.warning(f"{len(orphan_cells)} orphan cells are re-assigned.")
            interface_eids = np.append(interface_eids, orphan_cells)

            #! Central mesh object not updated. E.g. lose connection between part.element_ids and
            #! model.mesh.volume_ids/.cell_data["_volume-id"]

        if interface_eids.shape[0] == 0:
            LOGGER.warning(
                """Atria and ventricles do not seem to be
                connected, not generating a separate part for isolation."""
            )
            return None

        # create a new part
        isolation: Part = self.create_part_by_ids(interface_eids, "Atrioventricular isolation")
        isolation.part_type = PartType.ATRIUM
        isolation.fiber = True
        isolation.active = False
        isolation.ep_material = EPMaterial.Insulator()

        return isolation

    def create_stiff_ventricle_base(
        self,
        threshold_left_ventricle: float = 0.9,
        threshold_right_ventricle: float = 0.95,
        stiff_material: MechanicalMaterialModel = Mat295(
            rho=0.001, iso=ISO(itype=1, beta=2, kappa=10, mu1=0.1, alpha1=2)
        ),
    ) -> None | Part:
        """Use universal coordinates to generate a stiff base region.

        Parameters
        ----------
        threshold_left_ventricle : float, optional
            uvc_l larger than threshold will be set as stiff material, by default 0.9
        threshold_right_ventricle : float, optional
            a uvc_l value larger than this threshold in the right ventricle will be set to a stiff
            material, by default 0.95
        stiff_material : MechanicalMaterialModel, optional
            material to assign, by default MAT295(rho=0.001,
            iso=ISO(itype=1, beta=2, kappa=10, mu1=0.1, alpha1=2)

        Returns
        -------
        Part
            Part associated with the stiff base region.
        """
        try:
            v = self.mesh.point_data_to_cell_data()["apico-basal"]
        except KeyError:
            LOGGER.error("Array named 'apico-basal' cannot be found, cannot create base part.")
            LOGGER.error("Please call simulator.compute_uhc() first.")
            return

        eids = np.intersect1d(
            np.where(v > threshold_left_ventricle)[0], self.left_ventricle.element_ids
        )
        if not isinstance(self, LeftVentricle):
            # uvc-L of RV is generally smaller, *1.05 to be comparable with LV
            eid_r = np.intersect1d(
                np.where(v > threshold_right_ventricle)[0],
                self.right_ventricle.element_ids,
            )
            eids = np.hstack((eids, eid_r))

        part: Part = self.create_part_by_ids(eids, "base")
        part.part_type = PartType.VENTRICLE
        part.fiber = False
        part.active = False
        part.meca_material = stiff_material
        # assign default EP material as for ventricles
        part.ep_material = EPMaterial.Active()

        return part

    def create_atrial_stiff_ring(self, radius: float = 2) -> None | Part:
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
        # TODO: @mhoeijm move this to FourChamber class
        if not isinstance(self, FourChamber):
            LOGGER.error("This method is only for FourChamber model.")
            return

        # get ring cells from cap node list
        ring_nodes = []
        for cap in self.left_atrium.caps:
            # update cap mesh with up-to-date mesh
            cap._mesh = self.mesh.get_surface(cap._mesh.id)
            if cap.type is not CapType.MITRAL_VALVE_ATRIUM:
                ring_nodes.extend(cap.global_node_ids_edge.tolist())
        for cap in self.right_atrium.caps:
            # update cap mesh with up-to-date mesh
            cap._mesh = self.mesh.get_surface(cap._mesh.id)
            if cap.type is not CapType.TRICUSPID_VALVE_ATRIUM:
                ring_nodes.extend(cap.global_node_ids_edge.tolist())

        ring_eles = vtk_utils.find_cells_close_to_nodes(self.mesh, ring_nodes, radius=radius)
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
        ring: Part = self.create_part_by_ids(ring_eles, name="atrial stiff rings")
        ring.part_type = PartType.ATRIUM
        ring.fiber = False
        ring.active = False
        # assign default EP material
        ring.ep_material = EPMaterial.Active()

        return ring


class LeftVentricle(HeartModel):
    """Model of just the left ventricle."""

    def __init__(self, working_directory: pathlib.Path | str = None) -> None:
        self.left_ventricle: Part = Part(name="Left ventricle", part_type=PartType.VENTRICLE)
        """Left ventricle part."""
        # remove septum - not used in left ventricle only model
        del self.left_ventricle.septum

        self.left_ventricle.fiber = True
        self.left_ventricle.active = True

        super().__init__(working_directory=working_directory)
        pass


class BiVentricle(HeartModel):
    """Model of the left and right ventricle."""

    def __init__(self, working_directory: pathlib.Path | str = None) -> None:
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

        super().__init__(working_directory=working_directory)
        pass


class FourChamber(HeartModel):
    """Model of the left/right ventricle and left/right atrium."""

    def __init__(self, working_directory: pathlib.Path | str = None) -> None:
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

        super().__init__(working_directory=working_directory)

        pass


class FullHeart(FourChamber):
    """Model of both ventricles, both atria, aorta and pulmonary artery."""

    def __init__(self, working_directory: pathlib.Path | str = None) -> None:
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

        self.aorta.ep_material = EPMaterial.Insulator()
        self.pulmonary_artery.ep_material = EPMaterial.Insulator()

        super().__init__(working_directory=working_directory)

        pass


if __name__ == "__main__":
    print("Protected")
