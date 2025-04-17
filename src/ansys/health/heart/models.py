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
import pathlib
import re
from typing import List, Literal, Union

import numpy as np
import pyvista as pv
import yaml

from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.exceptions import InvalidHeartModelError
from ansys.health.heart.objects import (
    Cap,
    CapType,
    Cavity,
    Mesh,
    Part,
    PartType,
    Point,
    SurfaceMesh,
)
from ansys.health.heart.pre.conduction_beams import (
    ConductionPath,
)
from ansys.health.heart.pre.input import _InputModel
import ansys.health.heart.pre.mesher as mesher
from ansys.health.heart.settings.material.ep_material import EPMaterial
from ansys.health.heart.settings.material.material import (
    ISO,
    Mat295,
    MechanicalMaterialModel,
)
import ansys.health.heart.utils.connectivity as connectivity
import ansys.health.heart.utils.vtk_utils as vtk_utils


def _get_axis_from_field_data(
    mesh: Mesh | pv.UnstructuredGrid, axis_name: Literal["l4cv_axis", "l2cv_axis", "short_axis"]
) -> dict:
    """Get axis from the mesh field data."""
    try:
        return {
            "center": mesh.field_data[axis_name][0],
            "normal": mesh.field_data[axis_name][1],
        }
    except KeyError:
        LOGGER.info(f"Failed to retrieve {axis_name} from the mesh field data")
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
        LOGGER.info("Data has wrong shape. Expecting (2,3) shaped data.")
        return None
    mesh.field_data[axis_name] = data
    return mesh


def _set_workdir(workdir: pathlib.Path | str = None) -> str:
    """Set the root working directory.

    Parameters
    ----------
    workdir : pathlib.Path | str, default: None
        Path to the desired working directory.

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
        """List of parts."""
        parts = []
        for key, value in self.__dict__.items():
            attribute = getattr(self, key)
            if isinstance(attribute, Part):
                parts.append(attribute)
        return parts

    @property
    def part_names(self) -> List[str]:
        """List of part names."""
        part_names = []
        for part in self.parts:
            part_names.append(part.name)
        return part_names

    @property
    def part_ids(self) -> List[int]:
        """List of used part IDs."""
        return [part.pid for part in self.parts]

    @property
    def surfaces(self) -> List[SurfaceMesh]:
        """List of all defined surfaces."""
        return [s for p in self.parts for s in p.surfaces]

    @property
    def surface_names(self) -> List[str]:
        """List of all defined surface names."""
        return [s.name for s in self.surfaces]

    @property
    def surface_ids(self) -> List[str]:
        """List of all defined surface IDs."""
        return [s.id for s in self.surfaces]

    @property
    def cavities(self) -> List[Cavity]:
        """List of all cavities in the model."""
        return [part.cavity for part in self.parts if part.cavity]

    @property
    def part_name_to_part_id(self) -> dict:
        """Dictionary that maps the part name to the part ID."""
        return {p.name: p.pid for p in self.parts}

    @property
    def part_id_to_part_name(self) -> dict:
        """Dictionary that maps the part ID to the part name."""
        return {p.pid: p.name for p in self.parts}

    @property
    def surface_name_to_surface_id(self) -> dict:
        """Dictionary that maps the surface name to the surface ID."""
        return {s.name: s.id for p in self.parts for s in p.surfaces}

    @property
    def surface_id_to_surface_name(self) -> dict:
        """Dictionary that maps the surface ID to the surface name."""
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
        """List of cap centroids."""
        return [
            Point(name=c.name + "_center", xyz=c.centroid, node_id=c.global_centroid_id)
            for p in self.parts
            for c in p.caps
        ]

    def __init__(self, working_directory: pathlib.Path | str = None) -> None:
        """Initialize the heart model.

        Parameters
        ----------
        working_directory : pathlib.Path | default: None
            Path to desired working directory.

        Notes
        -----
        If no working directory is specified, the current working directory is used.
        """
        self.workdir = _set_workdir(working_directory)
        """Working directory."""

        self.mesh = Mesh()
        """Computational mesh."""

        self.fluid_mesh = Mesh()
        """Generated fluid mesh."""

        #! TODO: non-functional flag. Remove or replace.
        self._add_blood_pool: bool = False
        """Flag indicating whether a blood pool mesh is added. (Experimental)"""

        self._input: _InputModel = None
        """Input model."""

        self._add_subparts()
        """Add any subparts."""

        self._set_part_ids()
        """Set incremental part IDs."""

        self.electrodes: List[Point] = []
        """Electrodes positions for ECG computing."""

        self._conduction_paths: list[ConductionPath] = []
        """Conduction paths list."""

        self._conduction_mesh: Mesh = Mesh()
        """Mesh of merged conduction paths."""

        self.electrodes: List[Point] = []
        """Electrodes positions for ECG computing."""

        self._part_info = {}
        """Information about all parts in the model."""

        self._short_axis: dict = None
        """Short axis."""
        self._l2cv_axis: dict = None
        """l2cv axis."""
        self._l4cv_axis: dict = None
        """l4cv axis."""

        return

    @property
    def conduction_paths(self):
        """Return the list of conduction path."""
        return self._conduction_paths

    @property
    def conduction_mesh(self):
        """Return the conduction mesh."""
        return self._conduction_mesh

    def assign_conduction_paths(self, beams: ConductionPath | list[ConductionPath]):
        """Assign conduction paths to the model.

        Parameters
        ----------
        beams : ConductionBeams | list[ConductionBeams]
            list of conduction beams.

        Notes
        -----
        If conduction paths are already defined, they will be removed.
        """
        if len(self._conduction_paths) > 0:
            LOGGER.warning("Removing previously defined conduction beams.")
            self._conduction_paths: list[ConductionPath] = []
            self._conduction_mesh: Mesh = Mesh()

        if isinstance(beams, ConductionPath):
            beams = [beams]

        for beam in beams:
            self._conduction_paths.append(beam)

            # merge beam into conduction_system
            merge_ids, target_ids = self._find_merge_points(beam)
            self._conduction_mesh = self._safe_line_merge(
                self._conduction_mesh, beam.mesh, merge_ids, target_ids
            )

        # deduce the IDs of the conduction mesh in final mesh
        self._conduction_mesh.point_data["_shifted_id"] = self._shifted_id()

    def _find_merge_points(self, beam: ConductionPath):
        registered_name = [c.name for c in self._conduction_paths]

        merge_ids = []
        target_ids = []

        if beam._up is not None and beam._up in registered_name:
            target = next(i for i in self._conduction_paths if i.name == beam._up)

            LOGGER.info(
                f"merge first node of {beam.name.value} into closet point of {target.name.value}"
            )
            merge_ids.append(0)

            target_mesh = self._conduction_mesh.extract_cells(
                self._conduction_mesh["_line-id"] == target.id
            )
            sub_id = target_mesh.find_closest_point(beam.mesh.points[0])
            id2 = target_mesh["vtkOriginalPointIds"][sub_id]
            target_ids.append(id2)

        if beam._down is not None and beam._down in registered_name:
            target = next(i for i in self._conduction_paths if i.name == beam._down)
            LOGGER.info(
                f"merge last node of {beam.name.value} into closet point of {target.name.value}"
            )
            merge_ids.append(-1)
            target_mesh = self._conduction_mesh.extract_cells(
                self._conduction_mesh["_line-id"] == target.id
            )
            sub_id = target_mesh.find_closest_point(beam.mesh.points[-1])
            id2 = target_mesh["vtkOriginalPointIds"][sub_id]
            target_ids.append(id2)

        return merge_ids, target_ids

    def _shifted_id(self) -> np.ndarray:
        """
        Deduce node IDs after merging to solid mesh.

        TODO: move it to Mesh.
        """
        from scipy import spatial

        kdtree = spatial.cKDTree(self.mesh.points)

        is_connected = self._conduction_mesh["_is-connected"].astype(bool)
        querry_points = self._conduction_mesh.points[is_connected]
        dst, solid_id = kdtree.query(querry_points)
        LOGGER.info(f"Maximal distance from solid-beam connected node:{np.max(dst)}")

        shifted_ids = np.linspace(
            0, self._conduction_mesh.n_points - 1, num=self._conduction_mesh.n_points, dtype=int
        )
        # for connected nodes, replace by solid mesh ID
        shifted_ids[is_connected] = solid_id
        # for beam-only nodes, shift their IDs
        for i in range(self._conduction_mesh.n_points):
            if not is_connected[i]:
                shifted_ids[i] += self.mesh.n_points - np.sum(is_connected[:i])

        return shifted_ids

    @staticmethod
    def _safe_line_merge(base: Mesh, add_mesh: Mesh, mereged_id: list, target_id: list) -> Mesh:
        """
        Merge lines by explicitly define IDs to be merged.

        TODO: move it to Mesh.
        """

        def get_lines(m: pv.UnstructuredGrid | pv.PolyData):
            if m.GetNumberOfCells() == 0:
                return np.empty(shape=(0, 2), dtype=np.int_)
            if isinstance(m, pv.UnstructuredGrid):
                return m.cells.reshape(-1, 3)[:, 1:]
            elif isinstance(m, pv.PolyData):
                return m.lines.reshape(-1, 3)[:, 1:]

        base_points = base.points
        base_lines = get_lines(base)

        if base.GetNumberOfCells() == 0:
            point_data = np.empty(shape=(0,))
            cell_data = np.empty(shape=(0,))
        else:
            point_data = base.point_data["_is-connected"]
            cell_data = base.cell_data["_line-id"]

        if mereged_id == []:
            # no merge
            new_points = add_mesh.points
            new_point_data = add_mesh.point_data["_is-connected"]
            new_lines = get_lines(add_mesh) + len(base_points)

        elif mereged_id == [0]:
            new_points = add_mesh.points[1:]
            # first node is merged, lead to an offset of all lines
            new_lines = get_lines(add_mesh) + len(base_points) - 1
            # replace first node
            new_lines[0, 0] = target_id[0]

            # point data
            new_point_data = add_mesh.point_data["_is-connected"][1:]

        elif mereged_id == [0, -1]:
            # first node is merged, lead to an offset of all lines
            # last node is just dropped
            new_points = add_mesh.points[1:-1]
            new_lines = get_lines(add_mesh) + len(base_points) - 1
            # replace first node
            new_lines[0, 0] = target_id[0]
            # replace last node
            new_lines[-1, 1] = target_id[1]

            # point data
            new_point_data = add_mesh.point_data["_is-connected"][1:-1]
        else:
            NotImplementedError("Do not handle this merge lines.")

        merged = pv.PolyData()
        merged.points = np.vstack((base_points, new_points))
        merged_lines = np.vstack((base_lines, new_lines))
        merged.lines = np.hstack(
            (2 * np.ones(len(merged_lines), dtype=int)[:, np.newaxis], merged_lines)
        )
        merged.cell_data["_line-id"] = np.hstack((cell_data, add_mesh.cell_data["_line-id"]))
        merged.point_data["_is-connected"] = np.hstack((point_data, new_point_data))

        return Mesh(merged)

    def __str__(self):
        """Represent self as string."""
        return yaml.dump(self.summary(), sort_keys=False)

    # TODO: There is overlap with the input module.
    def _get_parts_info(self):
        """Get the ID to the model map that allows reconstructing the model from a mesh object."""
        for part in self.parts:
            self._part_info.update(part._get_info())
        return self._part_info

    def create_part_by_ids(self, eids: List[int], name: str) -> Union[None, Part]:
        """Create a part by element IDs.

        Parameters
        ----------
        eids : List[int]
            List of element IDs.
        name : str
            Part name.

        Returns
        -------
        Union[None, Part]
           Part if successful
        """
        if len(eids) == 0:
            LOGGER.error(f"Failed to create {name}. Element list is empty.")
            return None

        if name in [p.name for p in self.parts]:
            LOGGER.error(f"Failed to create {name}. Name already exists.")
            return None

        for part in self.parts:
            try:
                part.element_ids = np.setdiff1d(part.element_ids, eids)
            except ValueError:
                LOGGER.error(f"Failed to create part {name}.")
                return None

        self.add_part(name)
        new_part: Part = self.get_part(name)

        new_part.element_ids = eids

        return new_part

    def load_input(self, input_vtp: pv.PolyData, part_definitions: dict, scalar: str):
        """Load an input model.

        Parameters
        ----------
        input_vtp : pv.PolyData
            Input surface mesh represented by a VTK ``PolyData`` object.
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
            LOGGER.error("Failed to initialize input model. Check the input arguments.")
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
        use_wrapper : bool, default: False
            Whether to use the non-manifold mesher.
        overwrite_existing_mesh : bool, default: True
            Whether to overwrite the existing MSH.H5 mesh.
        global_mesh_size : float, default: 1.5
            Global mesh size for the generated mesh.
        path_to_fluent_mesh : str, default: None
            Path to the generated Fluent MSH.H5 mesh.
        mesh_size_per_part : dict, default: None
            Dictionary specifying the target mesh size for each part.
        _global_wrap_size : float, default: 1.5
            Global size for setting up the size-field for the shrink-wrap algorithm.
        _wrap_size_per_part : dict, default: None
            Per part size for setting up the size-field for the shrink-wrap algorithm.

        Examples
        --------
        >>> from ansys.health.heart.models import HeartModel
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
        When the input surfaces are non-manifold, the wrapper tries
        to reconstruct the surface and parts. Inevitably this leads to
        reconstruction errors. Nevertheless, in many instances this approach is
        more robust than meshing from a manifold surface. Moreover, any clear interface
        between parts is potentially lost.

        When the ``mesh_size_per_part`` is incomplete, remaining part sizes default to the
        global mesh size. This is an experimental setting. Any wrap sizes given
        as input arguments are ignored when the wrapper is not used.
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
        remesh_caps : bool, default: True
            Whether to remesh the caps of each cavity.
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
            LOGGER.debug(
                "Meshing of fluid cavities is not possible. No fluid surfaces are detected."
            )
            return

        if len(caps) == 0:
            LOGGER.debug("Meshing of fluid cavities is not possible. No caps are detected.")
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
        """Get a specific part based on a part name."""
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
        from ansys.health.heart.utils.misc import model_summary

        summary = model_summary(self)
        return summary

    def plot_mesh(self, show_edges: bool = True, color_by: str = "_volume-id"):
        """Plot the volume mesh of the heart model.

        Parameters
        ----------
        show_edges : bool, default: True
            Whether to plot the edges.
        color_by : str, default: ``'_volume-id'``
            Color by cell/point data.

        Examples
        --------
        >>> import ansys.health.heart.models as models
        >>> model = models.HeartModel.load_model("heart_model.pickle")
        >>> model.plot_mesh(show_edges=True)
        """
        plotter = pv.Plotter()
        plotter.add_mesh(self.mesh, show_edges=show_edges, scalars=color_by)

        plotter.show()
        return

    def plot_part(self, part: Part):
        """Plot a part in the mesh.

        Parameters
        ----------
        part : Part
            Part to highlight in the mesh.

        Examples
        --------
        >>> import ansys.health.heart.models as models
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
        plot_raw_mesh : bool, default: False
            Whether to plot the streamlines on the raw mesh.
        n_seed_points : int, default: 1000
            Number of seed points. While the default is ``1000``, using ``5000``
            is recommended.

        Examples
        --------
        >>> import ansys.health.heart.models as models
        >>> model = models.HeartModel.load_model("my_model.pickle")
        >>> model.plot_fibers(n_seed_points=5000)
        """
        plotter = pv.Plotter()

        # fiber direction is stored in cell data, but the cell-to-point filter
        # leads to issues, where nan values in any non-volume cell may change
        # the fiber direction in the target points.
        mesh = self.mesh.extract_cells_by_type([pv.CellType.TETRA, pv.CellType.HEXAHEDRON])
        mesh = mesh.ctp()
        streamlines = mesh.streamlines(vectors="fiber", source_radius=75, n_points=n_seed_points)
        if streamlines.n_cells == 0:
            LOGGER.error(
                "Failed to generate streamlines with radius {source_radius} and {n_seed_points}."
            )
            return None
        tubes = streamlines.tube()
        plotter.add_mesh(mesh, opacity=0.5, color="white")
        plotter.add_mesh(tubes, color="white")
        plotter.show()
        return plotter

    def plot_surfaces(self, show_edges: bool = True):
        """Plot all surfaces in the model.

        Examples
        --------
        Import modules and load model.
        >>> import ansys.health.heart.models as models
        >>> model = models.HeartModel.load_model("my_model.pickle")
        Plot the model
        >>> model.plot(show_edges=True)
        """
        try:
            import matplotlib as matplotlib
        except ImportError:
            LOGGER.warning("Matplotlib is not found. Install with 'pip install matplotlib'.")
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
        if self._conduction_mesh is None or self._conduction_mesh.number_of_cells == 0:
            LOGGER.info("No conduction system was found.")
            return

        try:
            plotter = pv.Plotter()
            plotter.add_mesh(self.mesh, color="w", opacity=0.1)
            self._conduction_mesh.set_active_scalars("_line-id")
            beams = self._conduction_mesh
            plotter.add_mesh(beams, line_width=2)
            plotter.show()
        except Exception:
            LOGGER.warning("Failed to plot the mesh.")
        return

    def save_model(self, filename: str):
        """Save the model and necessary information to reconstruct.

        Parameters
        ----------
        filename : str
            Path to the model.

        Notes
        -----
        The mesh of the heart model is saved as a VTU file. An
        additional ``partinfo.json`` file is written to reconstruct
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
        """Load a model from an existing VTU file and part information dictionary.

        Parameters
        ----------
        filename_mesh : str
            Path to the VTU file containing the mesh.
        filename_part_info : str
            Path to the JSON file that contains the part information for reconstructing the model.

        Examples
        --------
        >>> from ansys.health.heart.models import FullHeart
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
                LOGGER.warning(f"{part_1.name} is not in the part information.")
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
                LOGGER.warning(f"Failed to set element IDs for {part_1.name}.")
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
            LOGGER.warning("Heart is not defined in the VTU file.")
            try:
                LOGGER.warning("Computing heart axis...")
                self._define_anatomy_axis()
            except Exception:
                LOGGER.error(
                    "Failed to extract heart axis. Consider computing and setting manually."
                )
        else:
            LOGGER.info("Heart axis defined in the VTU file is reused...")

        return

    def _set_part_ids(self):
        """Populate part IDs."""
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
        """Get an array of used element IDs."""
        element_ids = np.empty(0, dtype=int)
        for part in self.parts:
            element_ids = np.append(element_ids, part.element_ids)

        return element_ids

    def _update_parts(self):
        """Update the parts using the meshed volume.

        Notes
        -----
        1. Extracts septum.
        2. Updates parts to include element IDs of the respective part.
        3. Assigns surfaces to each part.
        4. Extracts the closing caps.
        5. Creates cavities.
        6. Extracts apical points.
        7. Computes left-ventricle axis.
        8. Computes left-ventricle 17 segments.
        9. Adds nodal areas.
        10. Adds surface normals to boundaries.
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
            Overwrites the default part IDs by those given by the user.
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

        This method uses the septum surface of the right ventricle.
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
            raise InvalidHeartModelError(
                "Expecting only one surface that contains string: 'septum'"
            )
        if len(septum_name) == 0:
            raise InvalidHeartModelError("No boundary found with name: 'septum'")
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

        # remove these element ID from the left-ventricle
        part = next(part for part in self.parts if part.name == "Left ventricle")
        mask = np.isin(part.element_ids, element_ids_septum, invert=True)
        part.element_ids = part.element_ids[mask]

        return

    def _extract_apex(self, check_edge: bool = True) -> None:
        """Extract the apex of the ventricles.

        Notes
        -----
        Apex is defined as the point furthest from the mid-point between the cap/valves.

        Parameters
        ----------
        check_edge : bool, default: True
            Whether to check and correct if the apical point is on the edge of a surface.
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
                        f"Initial apical point is on edge of {surface.name}. The next closest point is used."  # noqa: E501
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
        """Get the element IDs of each part and assign these to the ``Part`` objects."""
        # get element IDs of each part.
        used_element_ids = self._get_used_element_ids()
        for part in self.parts:
            if len(part.element_ids) > 0:
                LOGGER.warning(
                    "Part {0} seems to already have elements assigned. Skipping.".format(part.name)
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
                                "Multiple candidate surfaces for septum exist. Using the first one."
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
        idoffset = 1000  # TODO: need to improve id checking
        ii = 0

        for part in self.parts:
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
                LOGGER.warning(f"Skipping part {part.name}. Only empty surfaces are present.")
                continue

            surface: SurfaceMesh = SurfaceMesh(pv.merge(surfaces))
            surface.name = part.name + " cavity"

            # save this cavity mesh to the centralized mesh object
            surface.id = int(np.sort(self.mesh.surface_ids)[-1] + 1)  # get unique ID.

            # Generate patches that close the surface.
            patches = vtk_utils.get_patches_with_centroid(surface)

            LOGGER.debug(f"Generating {len(patches)} caps for {part.name}")

            # TODO: Note that points come from surface, and does not contain all points in the mesh.
            # Create Cap objects with patches.
            for patch in patches:
                ii += 1

                cap_name = f"cap_{ii}_{part.name}"

                # create cap: NOTE, mostly for compatibility. Could simplify further
                cap_mesh = SurfaceMesh(patch.clean(), name=cap_name, id=ii + idoffset)

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
        """Try to update the cap types using the names of the connected boundaries."""
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
        """Validate that caps are attached to the right part."""
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
                    "Part: {0}. Cap types {1} are not in expected cap types:{2}".format(
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
                LOGGER.error(f"Surface {invalid_s.name} is empty.")
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
                LOGGER.error(f"Part {invalid_p.name} is empty.")
                is_valid = False

        return is_valid

    def _sync_epicardium_with_part(self):
        """Clean epicardial surfaces such that these use only nodes of the part."""
        for part in self.parts:
            self.mesh._set_global_ids()
            global_node_ids_part = self.mesh.extract_cells(part.element_ids).point_data[
                "_global-point-ids"
            ]

            # ! The only information we use from surface here is the ID, not the mesh information.
            # ! We need to go back to the central mesh to obtain an updated copy of
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

                    LOGGER.debug(f"Removing {np.sum(np.invert(mask))} faces from {surface.name}.")
                    surface.triangles = surface.triangles[mask, :]

                    # add updated mesh to global mesh.
                    # TODO: could just change the _surface-ids in the cell data directly.
                    # TODO: this basically appends the cells of surface at the end of Mesh.
                    self.mesh.remove_surface(surf_id)
                    self.mesh.add_surface(surface, int(surf_id), name=surface.name)

        return

    def _define_anatomy_axis(self):
        """Define long and short axes from left ventricle landmarks."""
        from ansys.health.heart.utils.landmark_utils import compute_anatomy_axis

        try:
            left_ventricle: Part = self.left_ventricle

        except AttributeError:
            LOGGER.info("Left ventricle part does not exist to build anatomical axis.")
            self.l2cv_axis = self.l4cv_axis = self.short_axis = {}
            return

        mv_center = next(
            cap.centroid for cap in left_ventricle.caps if cap.type == CapType.MITRAL_VALVE
        )
        av_center = next(
            cap.centroid for cap in left_ventricle.caps if cap.type == CapType.AORTIC_VALVE
        )
        apex = next(ap.xyz for ap in left_ventricle.apex_points if ap.name == "apex epicardium")

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
        """Get a nodeset around the apex point.

        Parameters
        ----------
        part : Literal["left", "right"], default: "left"
            On which part.
        option : Literal["endocardium", "epicardium", "myocardium"], default: "epicardium"
            On surface or in mesh.
        radius : float, default: 3
            Search in radius.

        Returns
        -------
        np.ndarray
            Apex nodeset
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
        Extract a layer of element to isolate between the ventricles and atrium.

        Notes
        -----
        These elements initially belong to the atrium.

        Returns
        -------
        Part
            Part of isolation elements.
        """
        # TODO: move this method to FourChamber class.
        if not isinstance(self, FourChamber):
            LOGGER.error("This method is only for the four-chamber heart model.")
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
            clean_obj = self.mesh.extract_cells(atrium.element_ids).connectivity(
                extraction_mode="largest"
            )
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
                connected. Not generating a separate part for isolation."""
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
        threshold_left_ventricle : float, default: 0.9
            If the ``uvc_l`` value is larger than this threshold in the left ventricle,
            it is set as stiff material.
        threshold_right_ventricle : float, default: 0.95
            If the ``uvc_l`` value is larger than this threshold in the right ventricle,
            it is set to a stiff
            material.
        stiff_material : MechanicalMaterialModel, default: MAT295(rho=0.001,
            iso=ISO(itype=1, beta=2, kappa=10, mu1=0.1, alpha1=2)
            Material to assign.

        Returns
        -------
        Part
            Part associated with the stiff base region.
        """
        try:
            v = self.mesh.point_data_to_cell_data()["apico-basal"]
        except KeyError:
            LOGGER.error("Array named 'apico-basal' is not found. Cannot create base part.")
            LOGGER.error("Call simulator.compute_uhc() first.")
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
        """Create a part for solids close to the atrial caps.

        Note
        ----
        Part created is passive and isotropic. The material must be defined.

        Parameters
        ----------
        radius : foat, default: 2
            Influence region.

        Returns
        -------
        Union[None, Part]
            Part of atrial rings if created.
        """
        # TODO: @mhoeijm move this to FourChamber class
        if not isinstance(self, FourChamber):
            LOGGER.error("This method is only for the four-chamber heart model.")
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
        largest = self.mesh.extract_cells(unselect_eles).connectivity(extraction_mode="largest")
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
    """Model of only the left ventricle."""

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
    """Model of the left and right ventricles."""

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
    """Model of both ventricles, both atria, the aorta, and the pulmonary artery."""

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
