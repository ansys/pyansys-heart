"""Module to manage input.

Notes
-----
This module manages the different types of input that can be handled and include:
1. User specified (volume) mesh.
2. User specified surface mesh. This will require remeshing.

Methods are provided to validate the volume and surface mesh objects (pyvista objects),
and to get the necessary parts or surfaces for each respective model.
"""

import copy
import os
from pathlib import Path
from typing import List, Union

from ansys.heart.custom_logging import LOGGER
import numpy as np
import pyvista as pv

# SURFACES_PER_HEART_PART stores the reference id's of the various parts and
# lists the required surfaces that enclose the part
# E.g. Left ventricle is enclosed by the endocardium, epicardium, and aortic/mitral valve interfaces
# hence to create a left-ventricle we need to provide these surfaces. The mesher will
# then generate the volume mesh for that part.
SURFACES_PER_HEART_PART = {
    "Left ventricle myocardium": {
        "id": 1,
        "enclosed_by_surfaces": {
            "left-ventricle-endocardium": 1,
            "left-ventricle-epicardium": 2,
            "interface@left-ventricle_aortic-valve": 3,
            "interface@left-ventricle_mitral-valve": 4,
        },
    },
    "Right ventricle myocardium": {
        "id": 2,
        "enclosed_by_surfaces": {
            "right-ventricle-endocardium": 5,
            "right-ventricle-epicardium": 6,
            "interface@right-ventricle_tricuspid-valve": 7,
            "interface@right-ventricle_pulmonary-valve": 8,
        },
    },
    "Left atrium myocardium": {
        "id": 3,
        "enclosed_by_surfaces": {
            "left-atrium-endocardium": 9,
            "left-atrium-epicardium": 10,
            "interface@left-atrium_left-atrium-appendage-inlet": 11,
            "interface@left-atrium_left-superior-pulmonary-vein-inlet": 12,
            "interface@left-atrium_left-inferior-pulmonary-vein-inlet": 13,
            "interface@right-atrium_right-inferior-pulmonary-vein-inlet": 14,
            "interface@right-atrium_right-superior-pulmonary-vein-inlet": 15,
        },
    },
    "Right atrium myocardium": {
        "id": 4,
        "enclosed_by_surfaces": {
            "right-atrium-endocardium": 16,
            "right-atrium-epicardium": 17,
            "interface@right-atrium_superior-vena-cava-inlet": 18,
            "interface@right-atrium_inferior-vena-cava-inlet": 19,
        },
    },
    "Aorta wall": {"id": 5, "enclosed_by_surfaces": {"aorta-wall": 20}},
    "Pulmonary artery wall": {"id": 6, "enclosed_by_surfaces": {"pulmonary-artery-wall": 21}},
}

# the different types of "basee" models supported
HEART_MODELS = {
    "LeftVentricle": ["Left ventricle myocardium"],
    "BiVentricle": ["Left ventricle myocardium", "Right ventricle myocardium"],
    "FourChamber": [
        "Left ventricle myocardium",
        "Right ventricle myocardium",
        "Left atrium myocardium",
        "Right atrium myocardium",
    ],
    "FullHeart": [
        "Left ventricle myocardium",
        "Right ventricle myocardium",
        "Left atrium myocardium",
        "Right atrium myocardium",
        "Aorta wall",
        "Pulmonary artery wall",
    ],
}


def _invert_dict(d: dict):
    """Invert the input dictionary."""
    return {v: k for k, v in d.items()}


def _get_required_parts(model_type: str) -> List[str]:
    """Get a list of required parts for the given model."""
    try:
        part_names = HEART_MODELS[model_type]
    except KeyError:
        raise KeyError(
            "{0} invalid heart model type. Valid model types include: {1}".format(
                model_type, list(HEART_MODELS.keys())
            )
        )
    return part_names


def _get_part_name_to_part_id_map() -> dict:
    """Get map that maps the part names to the part ids."""
    mapper = {}
    for k, value in SURFACES_PER_HEART_PART.items():
        mapper[k] = value["id"]
    return mapper


def _get_part_id_to_part_name_map() -> dict:
    """Get map that maps the part ids to the part names."""
    return _invert_dict(_get_part_name_to_part_id_map())


def _get_surface_name_to_surface_id_map() -> dict:
    """Get the map that maps the surface name to the surface id."""
    mapper = {}
    for part_name, part_subdict in SURFACES_PER_HEART_PART.items():
        mapper.update(part_subdict["enclosed_by_surfaces"])
    return mapper


def _get_surface_id_to_surface_name_map() -> dict:
    """Get the map that maps the surface name to the surface id."""
    return _invert_dict(_get_surface_name_to_surface_id_map())


def _get_required_surfaces(model_type: str) -> List[str]:
    """Return a list of surfaces required for the given model."""
    parts = _get_required_parts(model_type)
    required_surfaces = []
    for p in parts:
        required_surfaces += SURFACES_PER_HEART_PART[p]["enclosed_by_surfaces"]
    return required_surfaces


class InputManager:
    """Class to manage the different types of input.

    Notes
    -----
    Supported inputs include:
    1. Unstructured grid file or object with part-ids
    2. Multiblock VTK file or object with a single UnstructuredGrid block or
    with multiple PolyData objects
    3. PolyData file or object with surface-ids
    """

    def __init__(
        self,
        input: Union[Path, str, pv.UnstructuredGrid, pv.PolyData] = None,
        scalar: str = None,
        name_to_id_map: dict = None,
    ) -> None:
        """Read provided input volume or surface mesh.

        Parameters
        ----------
        input :  Union[Path, str, pv.UnstructuredGrid, pv.PolyData], optional
            An input volume mesh or surface mesh, either as a pyvista
            UnstructuredGrid object, pyvista PolyData object or path to vtk-like file,
            by default None
        scalar : str, optional
            Scalar array to use for either the part ids or surface ids, by default None
        name_to_id_map : dict, optional
            Map indicating which part/surface name corresponds to which part/surface id,
            by default None

        Examples
        --------
        Reading a UnstructuredGrid from a file and give the part-name to part-id map

        >>> mesh_file = "unstructured_grid.vtu" # unstructured grid where 'tags'
        ...                                       cell data represents the part-ids
        >>> input = InputManager(mesh_file, scalar="tags",
        ...             name_to_id_map={"Left ventricle" : 3, "Right ventricle" : 1})

        Reading a surface mesh (PolyData) from a file and explicitly give the surface
        name to surface-id map

        >>> mesh_file = "surface_mesh.vtk" # PolyData where 'cell-tags' represents the surface-ids
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
        self.input_surface: pv.PolyData = None
        """Input surface."""
        self._part_id_mapping = (
            _get_part_name_to_part_id_map() if not name_to_id_map else name_to_id_map
        )
        """Maps part-ids to part-names."""
        self._surface_id_mapping = (
            _get_surface_name_to_surface_id_map() if not name_to_id_map else name_to_id_map
        )
        """Maps surface-names to surface-ids."""

        # try to read the input.
        if isinstance(input, (Path, str)):
            LOGGER.info(f"Reading {input}...")
            if not os.path.isfile(input):
                raise FileNotFoundError(f"File {input} not found.")

        volume_set = False
        surface_set = False
        try:
            self.input_surface = pv.PolyData(input)
            surface_set = True
        except:
            try:
                self.input_volume = pv.UnstructuredGrid(input)
                volume_set = True
            except:
                pass

        if not volume_set and not surface_set:
            try:
                multi_block = pv.MultiBlock(input)
                if len(multi_block) == 1 and isinstance(multi_block[0], pv.UnstructuredGrid):
                    self.input_volume = multi_block[0]
                    volume_set = True
                elif len(multi_block) > 0:
                    raise NotImplementedError(
                        "Support for Multi-Block PolyData not yet implemented."
                    )
                    for ii, block in enumerate(multi_block):
                        if not isinstance(block, pv.PolyData):
                            raise ValueError("Expecting PolyData in MultiBlock with size > 1")
                        if ii == 0:
                            surface = block
                        else:
                            surface = surface.merge(block)
                    surface_set = True
            except:
                raise ImportError(f"Failed to load {input} as volume or surface.")

        # change array names if scalar is given.
        if self.input_volume and scalar:
            LOGGER.debug(f"Renaming {scalar} to part-id")
            self.input_volume.rename_array(scalar, "part-id")
        if self.input_surface and scalar:
            LOGGER.debug(f"Renaming {scalar} to surface-id")
            self.input_surface.rename_array(scalar, "surface-id")

        # validate
        self.validate()

        # reorder
        if volume_set and name_to_id_map:
            self._reorder_part_ids(name_to_id_map)

        if surface_set and name_to_id_map:
            self._reorder_surface_ids(name_to_id_map)

        pass

    def __repr__(self):
        """Represent self."""
        return (
            "Input volume mesh:\n"
            + str(self.input_volume)
            + "Input surface mesh:\n"
            + str(self.input_surface)
        )

    def _reorder_part_ids(self, part_name_to_part_id: dict):
        """Reorder the input part ids such that they correspond with SURFACES_PER_HEART_PART."""
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
            if key not in list(SURFACES_PER_HEART_PART.keys()):
                target_id = max_defined_id + 1
                max_defined_id += 1
            else:
                target_id = SURFACES_PER_HEART_PART[key]["id"]

            mask = old_ids == old_id
            new_ids[mask] = target_id
            self._part_id_mapping[key] = target_id

        self.input_volume.cell_data["part-id"] = new_ids

        return

    def _reorder_surface_ids(self, surface_name_to_surface_id: dict):
        """Reorder the input part ids such that they correspond with SURFACES_PER_HEART_PART."""
        old_ids = copy.deepcopy(self.input_surface.cell_data["surface-id"])
        new_ids = self.input_surface.cell_data["surface-id"]
        # reference map
        ref_map = _get_surface_name_to_surface_id_map()
        max_defined_id = max(ref_map.values())
        self._surface_id_mapping = {}

        if not np.all(np.isin(old_ids, list(surface_name_to_surface_id.values()))):
            raise ValueError(
                "Unable to map all surface ids to surface name: please extend dictionary"
            )

        for key, old_id in surface_name_to_surface_id.items():
            if key not in ref_map.keys():
                target_id = max_defined_id + 1
                max_defined_id += 1
            else:
                target_id = ref_map[key]

            mask = old_ids == old_id
            new_ids[mask] = target_id
            self._surface_id_mapping[key] = target_id

        self.input_surface.cell_data["surface-id"] = new_ids
        return

    def _validate_volume_mesh(self):
        """Perform some validation steps on the volume mesh."""
        if not "part-id" in self.input_volume.cell_data.keys():
            raise KeyError("Missing 'part-id' array in cell data.")
        # change part-id such that it always corresponds to the default part-ids.
        return

    def _validate_surface_mesh(self):
        """Perform some validation steps on the surface mesh."""
        if not "surface-id" in self.input_surface.cell_data.keys():
            raise KeyError("Missing 'surface-id' in cell-data.")

        # change surface-id such that it always corresponds to the default surface-ids.
        return

    def export_surfaces(self, format: str, folder: Union[Path, str] = ".") -> None:
        """Export the surfaces as separate stls."""
        from ansys.heart.preprocessor.mesh.misc import add_solid_name_to_stl

        surface_ids = np.unique(self.input_surface.cell_data["surface-id"])
        id_to_name = _get_surface_id_to_surface_name_map()

        for id in surface_ids:
            surface = self.input_surface.threshold([id - 1e-3, id + 1e-3], scalars="surface-id")
            surface_name = id_to_name[id]
            file_path = os.path.join(folder, surface_name + ".stl")
            surface.extract_surface().save(file_path)
            add_solid_name_to_stl(file_path, surface_name, file_type="binary")

    def validate(self):
        """Validate the given input."""
        if self.input_volume:
            self._validate_volume_mesh()
        if self.input_surface:
            self._validate_surface_mesh()

        return

    def is_valid_input(self):
        """Validate if the model has the right surfaces or parts defined for further processing."""
        is_valid = False
        try:
            self._validate_volume_mesh()
            is_valid = True
        except:
            pass
        try:
            self._validate_surface_mesh()
            is_valid = True
        except:
            pass
        return is_valid

    def get_required_parts_and_surfaces(self, model_type: str) -> dict:
        """Return a dictionary of the required parts and surfaces for a specific model."""
        parts = _get_required_parts(model_type)
        surfaces = _get_required_surfaces(model_type)
        LOGGER.info({"Parts": parts, "Surfaces": surfaces})

        return {"Parts": parts, "Surfaces": surfaces}

    def get_input(self):
        """Return the validated input volume or surface."""
        if isinstance(self.input_volume, pv.UnstructuredGrid):
            return self.input_volume
        elif isinstance(self.input_surface, pv.PolyData):
            return self.input_surface


PART_NAME_ID_MAPPING_DATABASES = {
    "Strocchi2020": {
        "Left ventricle myocardium": 1,
        "Right ventricle myocardium": 2,
        "Left atrium myocardium": 3,
        "Right atrium myocardium": 4,
        "Aorta wall": 5,
        "Pulmonary artery wall": 6,
        "Left atrial appendage border": 7,
        "Left superior pulmonary vein border": 8,
        "Left inferior pulmonary vein border": 9,
        "Right inferior pulmonary vein border": 10,
        "Right superior pulmonary vein border": 11,
        "Superior vena cava border": 12,
        "Inferior vena cava border": 13,
        "Mitral valve plane": 14,
        "Tricuspid valve plane": 15,
        "Aortic valve plane": 16,
        "Pulmonary valve plane": 17,
        "Left atrium appendage inlet": 18,
        "Left superior pulmonary vein inlet": 19,
        "Left inferior pulmonary vein inlet": 20,
        "Right inferior pulmonary vein inlet": 21,
        "Right superior pulmonary vein inlet": 22,
        "Superior vena cava inlet": 23,
        "Inferior vena cava inlet": 24,
    },
    "Rodero2021": {
        "Left ventricle myocardium": 1,
        "Right ventricle myocardium": 2,
        "Left atrium myocardium": 3,
        "Right atrium myocardium": 4,
        "Aorta wall": 5,
        "Pulmonary artery wall": 6,
        "Left atrial appendage border": 18,
        "Left superior pulmonary vein border": 21,
        "Left inferior pulmonary vein border": 20,
        "Right inferior pulmonary vein border": 19,
        "Right superior pulmonary vein border": 22,
        "Superior vena cava border": 23,
        "Inferior vena cava border": 24,
        "Mitral valve plane": 7,
        "Tricuspid valve plane": 8,
        "Aortic valve plane": 9,
        "Pulmonary valve plane": 10,
        "Left atrium appendage inlet": 11,
        "Left superior pulmonary vein inlet": 12,
        "Left inferior pulmonary vein inlet": 13,
        "Right inferior pulmonary vein inlet": 14,
        "Right superior pulmonary vein inlet": 15,
        "Superior vena cava inlet": 16,
        "Inferior vena cava inlet": 17,
    },
    "Strocchi2020_simplified": {
        "Left ventricle myocardium": 1,
        "Right ventricle myocardium": 2,
        "Mitral valve plane": 11,
        "Aortic valve plane": 12,
        "Pulmonary valve plane": 21,
        "Tricuspid valve plane": 22,
    },
}
SURFACE_NAME_ID_MAPPING_DATABASES = {
    "LabeledSurface": {
        "Left ventricle endocardium": 1,
        "Left ventricle epicardium": 2,
        "Left ventricle myocardium mitral valve": 3,
        "Left ventricle myocardium aortic valve": 4,
        "Right ventricle endocardium": 5,
        "Right ventricle epicardium": 6,
        "Right ventricle myocardium pulmonary valve": 7,
        "Right ventricle myocardium tricuspid valve": 8,
    },
}
