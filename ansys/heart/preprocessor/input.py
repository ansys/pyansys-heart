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

"""Module to manage input.

Notes
-----
This module manages the different types of input that can be handled and include:
1. User specified boundary mesh. This will require remeshing.

Methods are provided to validate the volume and boundary mesh objects (pyvista objects),
and to get the necessary parts or boundaries for each respective model.
"""

import copy
import os
from pathlib import Path
from typing import List, Union

from ansys.heart.custom_logging import LOGGER
from ansys.heart.preprocessor.mesh.vtkmethods import add_solid_name_to_stl
import numpy as np
import pyvista as pv
import yaml

# BOUNDARIES_PER_HEART_PART stores the reference id's of the various parts and
# lists the required boundaries that enclose the part
# E.g. Left ventricle is enclosed by the endocardium, epicardium, and aortic/mitral valve interfaces
# hence to create a left-ventricle we need to provide these boundaries. The mesher will
# then generate the volume mesh for that part.
BOUNDARIES_PER_HEART_PART = {
    "Left ventricle myocardium": {
        "id": 1,
        "enclosed_by_boundaries": {
            "left-ventricle-endocardium": 1,
            "left-ventricle-epicardium": 2,
            "interface@left-ventricle_aortic-valve": 3,
            "interface@left-ventricle_mitral-valve": 4,
        },
    },
    "Right ventricle myocardium": {
        "id": 2,
        "enclosed_by_boundaries": {
            "right-ventricle-endocardium": 6,
            "right-ventricle-epicardium": 7,
            "interface@right-ventricle_tricuspid-valve": 8,
            "interface@right-ventricle_pulmonary-valve": 9,
        },
    },
    "Septum": {
        "id": 3,
        "enclosed_by_boundaries": {
            "left-ventricle-septum": 10,
            "right-ventricle-septum": 11,
            "interface@left-ventricle_septum": 12,
        },
    },
    "Left atrium myocardium": {
        "id": 4,
        "enclosed_by_boundaries": {
            "left-atrium-endocardium": 13,
            "left-atrium-epicardium": 14,
            "interface@left-atrium_left-atrium-appendage-inlet": 15,
            "interface@left-atrium_left-superior-pulmonary-vein-inlet": 16,
            "interface@left-atrium_left-inferior-pulmonary-vein-inlet": 17,
            "interface@right-atrium_right-inferior-pulmonary-vein-inlet": 18,
            "interface@right-atrium_right-superior-pulmonary-vein-inlet": 19,
        },
    },
    "Right atrium myocardium": {
        "id": 4,
        "enclosed_by_boundaries": {
            "right-atrium-endocardium": 20,
            "right-atrium-epicardium": 21,
            "interface@right-atrium_superior-vena-cava-inlet": 22,
            "interface@right-atrium_inferior-vena-cava-inlet": 23,
        },
    },
    "Aorta wall": {"id": 6, "enclosed_by_boundaries": {"aorta-wall": 24}},
    "Pulmonary artery wall": {"id": 7, "enclosed_by_boundaries": {"pulmonary-artery-wall": 25}},
}

# the different types of "base" models supported
HEART_MODELS = {
    "LeftVentricle": ["Left ventricle myocardium"],
    "BiVentricle": ["Left ventricle myocardium", "Right ventricle myocardium", "Septum"],
    "FourChamber": [
        "Left ventricle myocardium",
        "Right ventricle myocardium",
        "Septum",
        "Left atrium myocardium",
        "Right atrium myocardium",
    ],
    "FullHeart": [
        "Left ventricle myocardium",
        "Right ventricle myocardium",
        "Septum",
        "Left atrium myocardium",
        "Right atrium myocardium",
        "Aorta wall",
        "Pulmonary artery wall",
    ],
}


class _InputBoundary(pv.PolyData):
    def __init__(
        self,
        var_inp=None,
        faces=None,
        n_faces=None,
        lines=None,
        n_lines=None,
        strips=None,
        n_strips=None,
        deep=False,
        force_ext=None,
        force_float=True,
        id=None,
        name: str = "",
    ) -> None:
        super().__init__(
            var_inp, faces, n_faces, lines, n_lines, strips, n_strips, deep, force_ext, force_float
        )
        self.id = id
        """ID of boundary."""
        self.name = name
        """Name of boundary."""

    @property
    def triangles(self):
        """Returns all triangles."""
        if not self.is_all_triangles:
            return
        else:
            return self.faces.reshape(self.n_cells, 4)[:, 1:]

    def __repr__(self):
        return f"Name:{self.name}\nid:{self.id}\n{super().__repr__()}"


class _InputPart:
    def __init__(self, name="", id=None, boundaries: List[_InputBoundary] = []) -> None:
        if not isinstance(boundaries, list):
            raise TypeError("Boundaries should be a list.")
        self.name = name
        """Name of part."""
        self.id = id
        """id of part."""
        self.boundaries: List[_InputBoundary] = boundaries
        """list of boundaries that enclose the part."""
        pass

    @property
    def boundary_ids(self):
        return [b.id for b in self.boundaries]

    @property
    def boundary_names(self):
        return [b.name for b in self.boundaries]

    @property
    def combined_boundaries(self):
        """Combined boundaries."""
        combined = pv.PolyData()
        for b in self.boundaries:
            combined += b
        return combined

    @property
    def is_manifold(self):
        """Flag indicating whether the part is manifold (watertight)."""
        return self.combined_boundaries.is_manifold

    def __repr__(self) -> str:
        return f"Name: {self.name}\nid:{self.id}\n" + "Number of boundaries:{len(self.boundaries)}"

    def write_boundaries(self, writedir: Union[str, Path] = ".", extension: str = ".stl"):
        """Write all boundaries of the part."""
        filenames = []
        for b in self.boundaries:
            filename = b.name.lower().replace(" ", "_")
            filename = os.path.join(writedir, filename + extension)
            b.save(filename)
            filenames += filename


class _InputModel:
    """Class to manage the different types of input.

    Notes
    -----
    Supported inputs include:
    1. [NotImplemented] Unstructured grid file or object with part-ids
    2. [NotImplemented] Multiblock VTK file or object with a single UnstructuredGrid block or
    with multiple PolyData objects
    3. PolyData file or object with boundary-ids
    """

    def __init__(
        self,
        input: Union[Path, str, pv.UnstructuredGrid, pv.PolyData] = None,
        scalar: str = None,
        part_definitions: dict = None,
    ) -> None:
        self.input_polydata: pv.PolyData = None
        """Input boundary."""
        self.part_definitions = part_definitions
        """Part definitions."""

        self._parts: List[_InputPart] = []

        # try to read the input.
        if isinstance(input, (Path, str)):
            LOGGER.info(f"Reading {input}...")
            if not os.path.isfile(input):
                raise FileNotFoundError(f"File {input} not found.")

        boundary_is_set = False
        try:
            self.input_polydata = pv.PolyData(input)
            boundary_is_set = True
        except:
            NotImplementedError(f"Failed to load file {input}. Other file types not supported yet.")
            return

        if self.input_polydata and scalar:
            LOGGER.debug(f"Renaming {scalar} to boundary-id")
            self.input_polydata.rename_array(scalar, "boundary-id")

        if part_definitions == None:
            return
            # raise NotImplementedError("Default part definitions not yet implemented.")

        self._add_parts(part_definitions)
        self._validate_input()
        self._validate_uniqueness()

        return

    @property
    def parts(self):
        """List of defined parts."""
        return self._parts

    @property
    def part_ids(self):
        """List of part ids."""
        return [p.id for p in self._parts]

    @property
    def part_names(self):
        """List of part names."""
        return [p.name for p in self._parts]

    @property
    def boundaries(self):
        """List of boundaries."""
        return [b for p in self._parts for b in p.boundaries]

    @property
    def boundary_names(self):
        """List of defined boundary names."""
        return [b.name for b in self.boundaries]

    @property
    def boundary_ids(self):
        """List of boundary ids."""
        return [b.id for b in self.boundaries]

    @property
    def as_single_polydata(self):
        """Combine all given boundaries into single polydata."""
        all_boundaries = pv.PolyData()
        for b in self.boundaries:
            all_boundaries += b
        return all_boundaries

    def __repr__(self):
        """Represent self."""
        return (
            "Input boundary mesh:\n" + str(self.input_polydata) + yaml.dump(self.part_definitions)
        )

    def _add_parts(self, part_definitions: dict):
        """Update the list of parts based on the part definitions."""
        is_visited = np.full(self.input_polydata.n_cells, False)

        for part_name in part_definitions.keys():
            boundaries = []
            for boundary_name, boundary_id in part_definitions[part_name][
                "enclosed_by_boundaries"
            ].items():
                mask = np.isin(self.input_polydata["boundary-id"], boundary_id)
                is_visited[mask] = True
                polydata = self.input_polydata.extract_cells(mask).extract_surface()

                b = _InputBoundary(polydata, name=boundary_name, id=boundary_id)
                boundaries.append(b)

            self._parts.append(
                _InputPart(
                    name=part_name, id=part_definitions[part_name]["id"], boundaries=boundaries
                )
            )

        # truncate input.
        if not np.all(is_visited):
            LOGGER.warning(
                "Not all faces are assigned to a boundary. Removing unused faces from input."
            )
            self.input_polydata = self.input_polydata.extract_cells(is_visited)

        return

    def _validate_if_parts_manifold(self):
        """Check if all parts are manifold (watertight)."""
        for p in self.parts:
            if not p.is_manifold:
                return False
        return True

    def _validate_uniqueness(self):
        """Validate whether there are any boundaries with duplicate ids or names."""
        is_valid = True
        # check id to name map
        mapper = {}
        for b in self.boundaries:
            if b.id not in list(mapper.keys()):
                mapper[b.id] = b.name
            else:
                if b.name != mapper[b.id]:
                    LOGGER.error(
                        "Boundary with id {0} has name {1} but expecting name {2}".format(
                            b.id, b.name, mapper[b.id]
                        )
                    )
                    is_valid = False
        # repeat for name to id map
        mapper = {}
        for b in self.boundaries:
            if b.name not in list(mapper.keys()):
                mapper[b.name] = b.id
            else:
                if b.id != mapper[b.name]:
                    LOGGER.error(
                        "Boundary with name {0} has id {1} but expecting id {2}".format(
                            b.name, b.id, mapper[b.name]
                        )
                    )
                    is_valid = False
        if not is_valid:
            LOGGER.warning("Please specify unique boundary name/id combination.")
        return is_valid

    def _validate_input(self):
        """Validate whether the provided scalars or list of scalars yield non-empty meshes."""
        if len(self.parts) == 0:
            LOGGER.warning("No parts defined, nothing to validate.")
            return None
        for part in self.parts:
            for b in part.boundaries:
                if not np.any(np.isin(self.input_polydata.cell_data["boundary-id"], b.id)):
                    raise ValueError(f"Boundary {b.name} with scalar {b.id} is empty.")
        return

    def plot(self, show_edges: bool = True):
        """Plot all boundaries."""
        try:
            import matplotlib as mpl
        except ImportError:
            LOGGER.error("Failed to import matplotlib. Install with pip install matplotlib.")
            return
        import matplotlib as mpl

        cmap = mpl.colormaps["tab20b"]

        plotter = pv.Plotter()
        for ii, b in enumerate(self.boundaries):
            plotter.add_mesh(
                b, show_edges=show_edges, color=cmap(ii / len(self.boundaries)), label=b.name
            )

        plotter.add_legend()
        plotter.show()

    def write_part_boundaries(
        self,
        writedir: Union[str, Path] = ".",
        extension: str = ".stl",
        avoid_duplicates: bool = True,
    ):
        """Write boundaries of all parts."""
        saved = []
        for b in self.boundaries:
            if avoid_duplicates:
                if b.id in saved:
                    continue
            filename = os.path.join(writedir, b.name.lower().replace(" ", "_") + extension)
            b.save(filename)
            add_solid_name_to_stl(filename, b.name, file_type="binary")
            saved.append(b.id)


def _invert_dict(d: dict):
    """Invert the input dictionary."""
    return {v: k for k, v in d.items()}


def _get_required_parts(model_type: str) -> dict:
    """Get a dict of required parts for the given model."""
    try:
        part_names = HEART_MODELS[model_type]
    except KeyError:
        raise KeyError(
            "{0} invalid heart model type. Valid model types include: {1}".format(
                model_type, list(HEART_MODELS.keys())
            )
        )
    parts = {}
    for p in part_names:
        parts[p] = _get_part_name_to_part_id_map()[p]
    return parts


def _get_part_name_to_part_id_map() -> dict:
    """Get map that maps the part names to the part ids."""
    mapper = {}
    for k, value in BOUNDARIES_PER_HEART_PART.items():
        mapper[k] = value["id"]
    return mapper


def _get_part_id_to_part_name_map() -> dict:
    """Get map that maps the part ids to the part names."""
    return _invert_dict(_get_part_name_to_part_id_map())


def _get_boundary_name_to_boundary_id_map() -> dict:
    """Get the map that maps the boundary name to the boundary id."""
    mapper = {}
    for part_name, part_subdict in BOUNDARIES_PER_HEART_PART.items():
        mapper.update(part_subdict["enclosed_by_boundaries"])
    return mapper


def _get_boundary_id_to_boundary_name_map() -> dict:
    """Get the map that maps the boundary name to the boundary id."""
    return _invert_dict(_get_boundary_name_to_boundary_id_map())


def _get_required_boundaries(model_type: str) -> List[str]:
    """Return a list of boundaries required for the given model."""
    parts = _get_required_parts(model_type)
    required_boundaries = []
    for p in parts:
        required_boundaries += BOUNDARIES_PER_HEART_PART[p]["enclosed_by_boundaries"]
    return required_boundaries


class InputManager:
    """Class to manage the different types of input.

    Notes
    -----
    Supported inputs include:
    1. Unstructured grid file or object with part-ids
    2. Multiblock VTK file or object with a single UnstructuredGrid block or
    with multiple PolyData objects
    3. PolyData file or object with boundary-ids
    """

    def __init__(
        self,
        input: Union[Path, str, pv.UnstructuredGrid, pv.PolyData] = None,
        scalar: str = None,
        name_to_id_map: dict = None,
    ) -> None:
        """Read provided input volume or boundary mesh.

        Parameters
        ----------
        input :  Union[Path, str, pv.UnstructuredGrid, pv.PolyData], optional
            An input volume mesh or boundary mesh, either as a pyvista
            UnstructuredGrid object, pyvista PolyData object or path to vtk-like file,
            by default None
        scalar : str, optional
            Scalar array to use for either the part ids or boundary ids, by default None
        name_to_id_map : dict, optional
            Map indicating which part/boundary name corresponds to which part/boundary id,
            by default None

        Examples
        --------
        Reading a UnstructuredGrid from a file and give the part-name to part-id map

        >>> mesh_file = "unstructured_grid.vtu" # unstructured grid where 'tags'
        ...                                       cell data represents the part-ids
        >>> input = InputManager(mesh_file, scalar="tags",
        ...             name_to_id_map={"Left ventricle myocardium" : 3,
        ...                             "Right ventricle myocardium" : 1})

        Reading a boundary mesh (PolyData) from a file and explicitly give the boundary
        name to boundary-id map

        >>> mesh_file = "boundary_mesh.vtk" # PolyData where 'cell-tags' represents the boundary-ids
        >>> input = InputManager(mesh_file, scalar="cell-tags",
            ...     name_to_id_map = {
        ...             "left-ventricle-endocardium": 3,
        ...             "left-ventricle-epicardium": 6,
        ...             "interface@left-ventricle_aortic-valve": 1,
        ...             "interface@left-ventricle_mitral-valve": 2})
        """
        # Try to populate these attributes during initialization.
        self.input_volume: pv.UnstructuredGrid = None
        """Input volume mesh."""
        self.input_boundary: pv.PolyData = None
        """Input boundary."""
        self._part_id_mapping = (
            _get_part_name_to_part_id_map() if not name_to_id_map else name_to_id_map
        )
        """Maps part-ids to part-names."""
        self._boundary_id_mapping = (
            _get_boundary_name_to_boundary_id_map() if not name_to_id_map else name_to_id_map
        )
        """Maps boundary-names to boundary-ids."""

        # try to read the input.
        if isinstance(input, (Path, str)):
            LOGGER.info(f"Reading {input}...")
            if not os.path.isfile(input):
                raise FileNotFoundError(f"File {input} not found.")

        volume_is_set = False
        boundary_is_set = False
        try:
            self.input_boundary = pv.PolyData(input)
            boundary_is_set = True
        except:
            try:
                self.input_volume = pv.UnstructuredGrid(input)
                volume_is_set = True
            except:
                pass

        if not volume_is_set and not boundary_is_set:
            try:
                multi_block = pv.MultiBlock(input)
                if len(multi_block) == 1 and isinstance(multi_block[0], pv.UnstructuredGrid):
                    self.input_volume = multi_block[0]
                    volume_is_set = True
                elif len(multi_block) > 0:
                    raise NotImplementedError(
                        "Support for Multi-Block PolyData not yet implemented."
                    )
                    for ii, block in enumerate(multi_block):
                        if not isinstance(block, pv.PolyData):
                            raise ValueError("Expecting PolyData in MultiBlock with size > 1")
                        if ii == 0:
                            boundary = block
                        else:
                            boundary = boundary.merge(block)
                    boundary_is_set = True
            except:
                raise ImportError(f"Failed to load {input} as volume or boundary.")

        # change array names if scalar is given.
        if self.input_volume and scalar:
            LOGGER.debug(f"Renaming {scalar} to part-id")
            self.input_volume.rename_array(scalar, "part-id")
        if self.input_boundary and scalar:
            LOGGER.debug(f"Renaming {scalar} to boundary-id")
            self.input_boundary.rename_array(scalar, "boundary-id")

        # validate
        self.validate()

        # reorder
        if volume_is_set and name_to_id_map:
            self._reorder_part_ids(name_to_id_map)

        if boundary_is_set and name_to_id_map:
            self._reorder_boundary_ids(name_to_id_map)

        pass

    def __repr__(self):
        """Represent self."""
        return (
            "Input volume mesh:\n"
            + str(self.input_volume)
            + "Input boundary mesh:\n"
            + str(self.input_boundary)
        )

    def _reorder_part_ids(self, part_name_to_part_id: dict):
        """Reorder the input part ids such that they correspond with BOUNDARIES_PER_HEART_PART."""
        old_ids = copy.deepcopy(self.input_volume.cell_data["part-id"])
        new_ids = self.input_volume.cell_data["part-id"]

        if not np.all(np.isin(old_ids, list(part_name_to_part_id.values()))):
            raise ValueError(
                "Unable to map all part ids to part name: please extend dictionary"
                + "with all defined part-ids",
            )

        max_defined_id = np.max(list(_get_part_name_to_part_id_map().values()))
        self._part_id_mapping = {}

        for key, old_id in part_name_to_part_id.items():
            if key not in list(BOUNDARIES_PER_HEART_PART.keys()):
                target_id = max_defined_id + 1
                max_defined_id += 1
            else:
                target_id = BOUNDARIES_PER_HEART_PART[key]["id"]

            mask = old_ids == old_id
            new_ids[mask] = target_id
            self._part_id_mapping[key] = target_id

        self.input_volume.cell_data["part-id"] = new_ids

        return

    def _reorder_boundary_ids(self, boundary_name_to_boundary_id: dict):
        """Reorder the input part ids such that they correspond with BOUNDARIES_PER_HEART_PART."""
        old_ids = copy.deepcopy(self.input_boundary.cell_data["boundary-id"])
        new_ids = self.input_boundary.cell_data["boundary-id"]
        # reference map
        ref_map = _get_boundary_name_to_boundary_id_map()
        max_defined_id = max(ref_map.values())
        self._boundary_id_mapping = {}

        if not np.all(np.isin(old_ids, list(boundary_name_to_boundary_id.values()))):
            raise ValueError(
                "Unable to map all boundary ids to boundary name: please extend dictionary"
            )

        for key, old_id in boundary_name_to_boundary_id.items():
            if key not in ref_map.keys():
                target_id = max_defined_id + 1
                max_defined_id += 1
            else:
                target_id = ref_map[key]

            mask = old_ids == old_id
            new_ids[mask] = target_id
            self._boundary_id_mapping[key] = target_id

        self.input_boundary.cell_data["boundary-id"] = new_ids
        return

    def _validate_volume_mesh(self):
        """Perform some validation steps on the volume mesh."""
        if not "part-id" in self.input_volume.cell_data.keys():
            raise KeyError("Missing 'part-id' array in cell data.")
        return

    def _validate_boundary_mesh(self):
        """Perform some validation steps on the boundary mesh."""
        if not "boundary-id" in self.input_boundary.cell_data.keys():
            raise KeyError("Missing 'boundary-d' in cell-data.")

        if not self.input_boundary.is_manifold:
            raise ImportWarning("Input boundary has gaps and is not watertight.")
        return

    def export_boundaries(self, format: str, folder: Union[Path, str] = ".") -> None:
        """Export the boundaries as separate stls."""
        from ansys.heart.preprocessor.mesh.misc import add_solid_name_to_stl

        boundary_ids = np.unique(self.input_boundary.cell_data["boundary-id"])
        id_to_name = _get_boundary_id_to_boundary_name_map()

        for id in boundary_ids:
            boundary = self.input_boundary.threshold([id - 1e-3, id + 1e-3], scalars="boundary-id")
            boundary_name = id_to_name[id]
            file_path = os.path.join(folder, boundary_name + ".stl")
            boundary.extract_surface().save(file_path)
            add_solid_name_to_stl(file_path, boundary_name, file_type="binary")

    def validate(self):
        """Validate the given input."""
        if self.input_volume:
            self._validate_volume_mesh()
        if self.input_boundary:
            self._validate_boundary_mesh()

        return

    def is_valid_input(self):
        """Validate if the model has the proper boundaries or parts defined."""
        is_valid = False
        try:
            self._validate_volume_mesh()
            is_valid = True
        except:
            pass
        try:
            self._validate_boundary_mesh()
            is_valid = True
        except:
            pass
        return is_valid

    def get_required_parts_and_boundaries(self, model_type: str) -> dict:
        """Return a dictionary of the required parts and boundaries for a specific model."""
        parts = _get_required_parts(model_type)
        boundaries = _get_required_boundaries(model_type)
        LOGGER.info({"Parts": parts, "Boundaries": boundaries})

        return {"Parts": parts, "Boundaries": boundaries}

    def get_input(self):
        """Return the validated input volume or boundary."""
        if isinstance(self.input_volume, pv.UnstructuredGrid):
            return self.input_volume
        elif isinstance(self.input_boundary, pv.PolyData):
            return self.input_boundary
