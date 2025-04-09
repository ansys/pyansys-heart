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

"""Module to manage input.

Notes
-----
This module manages the different types of input that can be handled and include:
1. User specified boundary mesh. This will require remeshing.

Methods are provided to validate the volume and boundary mesh objects (pyvista objects),
and to get the necessary parts or boundaries for each respective model.
"""

import os
from pathlib import Path

import numpy as np
import pyvista as pv
import yaml

from ansys.heart.core import LOG as LOGGER
from ansys.heart.core.exceptions import InvalidInputModelTypeError
from ansys.heart.core.utils.vtk_utils import add_solid_name_to_stl

# BOUNDARIES_PER_HEART_PART stores the reference id's of the various parts and
# lists the required boundaries that enclose the part
# E.g. Left ventricle is enclosed by the endocardium, epicardium, and aortic/mitral valve interfaces
# hence to create a left-ventricle we need to provide these boundaries. The mesher will
# then generate the volume mesh for that part.
_BOUNDARIES_PER_HEART_PART = {
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
_HEART_MODELS = {
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

    def __repr__(self):
        return f"Name:{self.name}\nid:{self.id}\n{super().__repr__()}"


class _InputPart:
    def __init__(self, name="", id=None, boundaries: list[_InputBoundary] = []) -> None:
        if not isinstance(boundaries, list):
            raise TypeError("Boundaries should be a list.")
        self.name = name
        """Name of part."""
        self.id = id
        """id of part."""
        self.boundaries: list[_InputBoundary] = boundaries
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

    def write_boundaries(self, writedir: str | Path = ".", extension: str = ".stl") -> None:
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
        input: Path | str | pv.UnstructuredGrid | pv.PolyData = None,
        scalar: str = None,
        part_definitions: dict = None,
    ) -> None:
        self.input_polydata: pv.PolyData = None
        """Input boundary."""
        self.part_definitions = part_definitions
        """Part definitions."""

        self._parts: list[_InputPart] = []

        # try to read the input.
        if isinstance(input, (Path, str)):
            LOGGER.info(f"Reading {input}...")
            if not os.path.isfile(input):
                raise FileNotFoundError(f"File {input} not found.")

        try:
            self.input_polydata = pv.PolyData(input)
        except Exception as e:
            raise NotImplementedError(f"Failed to load file {input}. {e}")

        if part_definitions is None:
            LOGGER.error("Please specify part definitions.")
            return None
        if scalar is None:
            LOGGER.error(
                "Please specify a scalar that is used to identify the enclosing boundaries."
            )
            return None

        if scalar != "boundary-id":
            if scalar in self.input_polydata.cell_data.keys():
                LOGGER.info(f"Renaming {scalar} to boundary-id")
                self.input_polydata.rename_array(scalar, "boundary-id")
            else:
                LOGGER.error(f"Failed to rename {scalar} to boundary-id")
                return None

        self.part_definitions = self._add_parts(part_definitions)
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

    def _add_parts(self, part_definitions: dict) -> dict:
        """Update the list of parts based on the part definitions.

        Parameters
        ----------
        part_definitions : dict
            Dictionary defining the part names and list of boundaries it is enclosed by.

        Returns
        -------
        dict
            Updated part definitions.

        Notes
        -----
        The part definitions only change if multiple ids are given for a single surface.
        The first item in the list defines the new boundary id, and the others are merged.

        """
        is_visited = np.full(self.input_polydata.n_cells, False)

        for part_name in part_definitions.keys():
            boundaries = []
            for boundary_name, boundary_id in part_definitions[part_name][
                "enclosed_by_boundaries"
            ].items():
                mask = np.isin(self.input_polydata["boundary-id"], boundary_id)
                is_visited[mask] = True
                polydata = self.input_polydata.extract_cells(mask).extract_surface()

                # take on first value in list
                if isinstance(boundary_id, list):
                    boundary_id = boundary_id[0]

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

        return part_definitions

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
        except ImportError as error:
            LOGGER.error(
                f"Failed to import matplotlib. Install with pip install matplotlib. {error}"
            )
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
        writedir: str | Path = ".",
        extension: str = ".stl",
        avoid_duplicates: bool = True,
        add_name_to_header: bool = True,
    ) -> tuple[list, list]:
        """Write boundaries of all parts."""
        saved = []
        boundary_names = []
        filenames = []
        ii = 0
        for b in self.boundaries:
            if avoid_duplicates:
                if b.id in saved:
                    continue
            filename = os.path.join(writedir, b.name.lower().replace(" ", "_") + extension)
            b.save(filename)
            if add_name_to_header:
                boundary_name = b.name
            else:
                boundary_name = ""

            add_solid_name_to_stl(filename, boundary_name, file_type="binary")
            saved.append(b.id)

            boundary_names += [boundary_name]
            filenames += [filename]
            ii += 1

        return filenames, boundary_names


def _invert_dict(d: dict) -> dict:
    """Invert the input dictionary."""
    return {v: k for k, v in d.items()}


def _get_required_parts(model_type: str) -> dict:
    """Get a dict of required parts for the given model."""
    try:
        part_names = _HEART_MODELS[model_type]
    except KeyError:
        raise InvalidInputModelTypeError(
            "{0} invalid heart model type. Valid model types include: {1}".format(
                model_type, list(_HEART_MODELS.keys())
            )
        )
    parts = {}
    for p in part_names:
        parts[p] = _get_part_name_to_part_id_map()[p]
    return parts


def _get_part_name_to_part_id_map() -> dict:
    """Get map that maps the part names to the part ids."""
    mapper = {}
    for k, value in _BOUNDARIES_PER_HEART_PART.items():
        mapper[k] = value["id"]
    return mapper


def _get_part_id_to_part_name_map() -> dict:
    """Get map that maps the part ids to the part names."""
    return _invert_dict(_get_part_name_to_part_id_map())


def _get_boundary_name_to_boundary_id_map() -> dict:
    """Get the map that maps the boundary name to the boundary id."""
    mapper = {}
    for part_name, part_subdict in _BOUNDARIES_PER_HEART_PART.items():
        mapper.update(part_subdict["enclosed_by_boundaries"])
    return mapper


def _get_boundary_id_to_boundary_name_map() -> dict:
    """Get the map that maps the boundary name to the boundary id."""
    return _invert_dict(_get_boundary_name_to_boundary_id_map())


def _get_required_boundaries(model_type: str) -> list[str]:
    """Return a list of boundaries required for the given model."""
    parts = _get_required_parts(model_type)
    required_boundaries = []
    for p in parts:
        required_boundaries += _BOUNDARIES_PER_HEART_PART[p]["enclosed_by_boundaries"]
    return required_boundaries
