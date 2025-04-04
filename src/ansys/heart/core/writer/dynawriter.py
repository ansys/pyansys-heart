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

"""Module contain. classes for writing LS-DYNA keywords based.

Notes
-----
Uses a HeartModel (from ansys.heart.core.models).

"""

import copy
from enum import Enum
import json

# import missing keywords
import os
import time
from typing import Callable, List, Literal, NamedTuple, Union

import numpy as np
import pandas as pd
import pyvista as pv
import scipy.spatial as spatial

from ansys.dyna.core.keywords import keywords
from ansys.heart.core import LOG as LOGGER
from ansys.heart.core.models import (
    BiVentricle,
    FourChamber,
    FullHeart,
    HeartModel,
    LeftVentricle,
)
from ansys.heart.core.objects import Cap, CapType, Part, PartType, SurfaceMesh, _ConductionType
from ansys.heart.core.utils.vtk_utils import compute_surface_nodal_area_pyvista
from ansys.heart.core.writer import custom_keywords as custom_keywords
from ansys.heart.core.writer.define_function_templates import (  # noqa F401
    _define_function_0d_system,
    _ed_load_template,
)
from ansys.heart.core.writer.heart_decks import (
    BaseDecks,
    ElectroMechanicsDecks,
    ElectrophysiologyDecks,
    FiberGenerationDecks,
    MechanicsDecks,
    PurkinjeGenerationDecks,
)
from ansys.heart.core.writer.keyword_utils import (
    add_beams_to_kw,
    add_nodes_to_kw,
    create_define_curve_kw,
    create_define_sd_orientation_kw,
    create_discrete_elements_kw,
    create_element_shell_keyword,
    create_element_solid_ortho_keyword,
    create_elemetn_solid_keyword,
    create_node_keyword,
    create_node_set_keyword,
    create_segment_set_keyword,
    fast_element_writer,
    get_list_of_used_ids,
)
from ansys.heart.core.writer.material_keywords import MaterialHGOMyocardium, MaterialNeoHook
from ansys.heart.simulator.settings.material.ep_material import CellModel, EPMaterial
from ansys.heart.simulator.settings.material.material import (
    Mat295,
    MechanicalMaterialModel,
    NeoHookean,
)
from ansys.heart.simulator.settings.settings import SimulationSettings, Stimulation


class _BoundaryConditionType(Enum):
    """Boundary condition type."""

    FIX = "fix"
    ROBIN = "Robin"


class CVInteraction(NamedTuple):
    """Template to define control volume interaction."""

    id: int
    cvid1: int
    cvid2: int
    lcid: int
    name: str
    parameters: dict


class ControlVolume(NamedTuple):
    """Template to define control volume."""

    part: Part
    id: int
    Interactions: list[CVInteraction]


class BaseDynaWriter:
    """Base class that contains essential features for all LS-DYNA heart models."""

    def __init__(self, model: HeartModel, settings: SimulationSettings = None) -> None:
        """Initialize writer by loading a HearModel and the desired settings.

        Parameters
        ----------
        model : HeartModel
            HeartModel object which contains the necessary information for the writer,
            such as nodes, elements, and parts
        settings : SimulationSettings, optional
            Simulation settings used to create the LS-DYNA model.
            Loads defaults if None, by default None

        Example
        -------
        TODO: add example
        """
        self.model = model
        """Model information necessary for creating the LS-DYNA .k files."""

        self.kw_database = BaseDecks()

        # These are general attributes useful for keeping track of ids:
        self.max_node_id: int = 0
        """Max node id."""
        self._used_part_ids: List[int] = []

        self.section_ids = []
        """List of used section ids."""
        self.mat_ids = []
        """List of used mat ids."""
        # self.volume_mesh = {
        #     "nodes": np.empty(0),
        #     "tetra": np.empty(0),
        #     "cell_data": {},
        #     "point_data": {},
        # }
        self.volume_mesh = model.mesh
        """Volume mesh information."""

        # keeps track of some element id offsets
        self.id_offset = {
            "part": 0,
            "section": 0,
            "material": 0,
            "vector": 0,
            "element": {"solid": 0, "discrete": 0, "shell": 0},
        }
        """Id offset for several relevant keywords."""

        #! Do we really need the below?
        for part in self.model.parts:
            if not part.pid:
                part.pid = np.max([p.pid for p in self.model.parts if p.pid]) + 1

        self.id_offset["part"] = np.max(self.model.part_ids)

        # ! Removed the below since the part ids in self.model.parts are already defined.
        # for part in self.model.parts:
        #     id += 1
        #     # cannot use get_unique_part_id() because it checks in Deck()
        #     # part.pid = self.get_unique_part_id()
        #     # part.pid = id
        # """Assign part id for heart parts."""

        """List of .k files to include in main. This is derived from the Decks classes."""
        self.include_files = []

        if not settings:
            self.settings = SimulationSettings()
            """Simulation settings."""
            LOGGER.warning("No settings provided - loading default values.")
            self.settings.load_defaults()

        else:
            self.settings = settings
            """Simulation settings."""

        self.settings.to_consistent_unit_system()
        self._check_settings()

        return

    def _check_settings(self):
        """Check if required settings are available."""
        import ansys.heart.simulator.settings.settings as sett

        subsettings_classes = [
            getattr(self.settings, attr).__class__
            for attr in self.settings.__dict__
            if isinstance(getattr(self.settings, attr), sett.Settings)
        ]

        if isinstance(self, MechanicsDynaWriter):
            if sett.Mechanics not in subsettings_classes:
                raise ValueError("Expecting mechanics settings.")

        elif isinstance(self, FiberGenerationDynaWriter):
            if sett.Fibers not in subsettings_classes:
                raise ValueError("Expecting fiber settings.")

        elif isinstance(self, PurkinjeGenerationDynaWriter):
            if sett.Purkinje not in subsettings_classes:
                raise ValueError("Expecting Purkinje settings.")

        elif isinstance(self, ElectrophysiologyDynaWriter):
            if sett.Electrophysiology not in subsettings_classes:
                raise ValueError("Expecting electrophysiology settings.")
        elif isinstance(self, LaplaceWriter):
            pass
        else:
            raise NotImplementedError(
                f"Checking settings for {self.__class__.__name__} not yet implemented."
            )

        return

    def _update_node_db(self, ids: np.ndarray = None):
        """Update node database.

        Parameters
        ----------
        ids : np.ndarray, optional
            0-based ids of nodes to write, by default None
        """
        LOGGER.debug("Updating node keywords...")
        node_kw = keywords.Node()
        if ids is not None:
            nodes = np.vstack([ids + 1, self.model.mesh.points[ids, :].T]).T
            node_kw = add_nodes_to_kw(nodes, node_kw)
        else:
            node_kw = add_nodes_to_kw(self.model.mesh.points, node_kw)

        self.kw_database.nodes.append(node_kw)

        return

    def _update_parts_db(self):
        """Loop over parts defined in the model and creates keywords."""
        LOGGER.debug("Updating part keywords...")

        # add parts with a dataframe
        section_id = self.get_unique_section_id()

        # get list of cavities from model
        for part in self.model.parts:
            # material ID = part ID
            part.mid = part.pid

            part_df = pd.DataFrame(
                {
                    "heading": [part.name],
                    "pid": [part.pid],
                    "secid": [section_id],
                    "mid": [part.mid],
                }
            )
            part_kw = keywords.Part()
            part_kw.parts = part_df

            self.kw_database.parts.append(part_kw)

        # set up section solid for cavity myocardium
        section_kw = keywords.SectionSolid(secid=section_id, elform=13)

        self.kw_database.parts.append(section_kw)

        return

    def _update_segmentsets_db(self, add_caps: bool = False, add_cavities: bool = True):
        """Update the segment set database."""
        # NOTE 0: add all surfaces as segment sets
        # NOTE 1: need to more robustly check segids that are already used?

        # add closed cavity segment sets
        if add_cavities:
            cavities = [p.cavity for p in self.model.parts if p.cavity]
            for cavity in cavities:
                #! Get up to date surface mesh of cavity.
                surface = self.model.mesh.get_surface(cavity.surface.id)
                segset_id = self.get_unique_segmentset_id()

                #! recompute normals: point normals may have changed
                #! do we need some check to ensure normals are pointing inwards?
                #! Could use surface.force_normals_inwards()
                surface.force_normals_inwards()

                cavity.surface._seg_set_id = segset_id
                kw = create_segment_set_keyword(
                    segments=surface.triangles_global + 1,
                    segid=cavity.surface._seg_set_id,  # TODO: replace
                    title=surface.name,
                )
                # append this kw to the segment set database
                self.kw_database.segment_sets.append(kw)

        # write surfaces as segment sets
        for part in self.model.parts:
            for surface in part.surfaces:
                surface_global = self.model.mesh.get_surface(surface.id)
                if not surface_global:
                    LOGGER.debug(f"Failed to create segment set for {surface.name}")
                    continue
                if surface_global.n_cells == 0:
                    LOGGER.debug(f"Failed to create segment set for {surface.name}. Empty mesh.")
                    continue

                segset_id = self.get_unique_segmentset_id()
                surface._seg_set_id = segset_id

                kw = create_segment_set_keyword(
                    segments=surface_global.triangles_global + 1,
                    segid=segset_id,
                    title=surface.name,
                )
                # append this kw to the segment set database
                self.kw_database.segment_sets.append(kw)

        if add_caps:
            # create corresponding segment sets
            caps = [cap for part in self.model.parts for cap in part.caps]
            for cap in caps:
                cap_mesh = self.model.mesh.get_surface(cap._mesh.id)
                segid = self.get_unique_segmentset_id()
                cap._mesh._seg_set_id = segid
                cap._seg_set_id = segid
                segset_kw = create_segment_set_keyword(
                    segments=cap_mesh.triangles_global + 1,
                    segid=cap._seg_set_id,
                    title=cap.name,
                )
                self.kw_database.segment_sets.append(segset_kw)
        return

    def _filter_bc_nodes(self, surface: SurfaceMesh):
        """Remove one or more nodes from tetrahedrons having all nodes in the boundary.

        Notes
        -----
        The removed node must be connected with at least 1 node outside the boundary, see #656.

        Parameters
        ----------
        surface : SurfaceMesh
            Boundary surface to be analysed.

        Returns
        -------
        node_ids : np.ndarray
            Array of boundary nodes after problematic node removal.
        """
        # getting elements in active parts
        element_ids = np.array([], dtype=int)
        node_ids = surface.global_node_ids_triangles

        for part in self.model.parts:
            element_ids = np.append(element_ids, part.element_ids)

        element_ids = np.unique(element_ids)
        active_tets = self.model.mesh.tetrahedrons[element_ids]

        # make sure not all nodes of the same elements are in the boundary
        node_mask = np.zeros(self.model.mesh.number_of_points, dtype=int)
        # tag boundary nodes with value 1
        node_mask[node_ids] = 1

        tet_mask = np.array(
            [
                node_mask[active_tets[:, 0]],
                node_mask[active_tets[:, 1]],
                node_mask[active_tets[:, 2]],
                node_mask[active_tets[:, 3]],
            ]
        )

        # getting tets with 4 nodes in boundary
        issue_tets = np.where(np.sum(tet_mask, axis=0) == 4)[0]

        # getting corresponding nodes
        issue_nodes = active_tets[issue_tets, :]

        # counting node appearances
        u_active_tets, tet_count_active = np.unique(active_tets, return_counts=True)
        u_issue_nodes, tet_count_issue = np.unique(issue_nodes, return_counts=True)

        # finding issue nodes that belong to at least one non-issue tet
        removable_mask = np.array(
            [
                tet_count_active[np.where(u_active_tets == ii)[0][0]]
                != tet_count_issue[np.where(u_issue_nodes == ii)[0][0]]
                for ii in issue_nodes.flatten()
            ]
        ).reshape(-1, 4)

        # removing the first issue node belonging to at least one non-issue tet (for each tet)
        column_idxs = np.argmax(removable_mask, axis=1)
        nodes_toremove = np.unique(
            [issue_nodes[ii, column_idxs[ii]] for ii in range(len(issue_tets))]
        )

        # checking that there are no nodes that only belong to non-issue tets
        if not np.all(np.any(removable_mask, axis=1)):
            # removing all such nodes and all their neighbors
            unsolvable_nodes = np.unique(issue_nodes[np.where(~np.any(removable_mask, axis=1))[0]])
            #! NOTE: surface.point_neighbors uses local indexing, so should get local index
            #! from global indices.
            local_point_ids = np.where(
                np.isin(surface.point_data["_global-point-ids"], unsolvable_nodes)
            )[0]
            local_unsolvable_nodes = np.unique(
                [
                    neighbor
                    for ii, node in enumerate(unsolvable_nodes)
                    for neighbor in surface.point_neighbors(local_point_ids[ii])
                ]
            )
            global_unsolvable_nodes = surface.point_data["_global-point-ids"][
                local_unsolvable_nodes
            ]
            nodes_toremove = np.append(nodes_toremove, global_unsolvable_nodes)

        node_ids = np.setdiff1d(node_ids, nodes_toremove)

        for cell in issue_tets:
            LOGGER.warning(
                f"All nodes of cell {cell + 1} are in nodeset of {surface.name},"
                + " removing at least one node."
            )

        return node_ids

    def _update_nodesets_db(
        self, remove_duplicates: bool = True, remove_one_node_from_cell: bool = False
    ):
        """Update the node set database.

        Parameters
        ----------
        remove_duplicates : bool, optional
            Remove nodes if they are used in other nodeset, by default True
        remove_one_node_from_cell : bool, optional
            Remove a node if a cell has all nodes in nodeset, by default False

        Notes
        -----
            In FiberGenerationWriter, we do not allow all nodes of same element in one nodeset.
        """
        # formats endo, epi- and septum nodeset keywords, do for all surfaces
        # for each surface in each part add the respective node-set
        # Use same ID as surface
        # TODO: check if database already contains nodesets (there will be duplicates otherwise)
        used_node_ids = np.empty(0, dtype=int)

        # add node-set for each cap
        for part in self.model.parts:
            for cap in part.caps:
                # update cap mesh:
                cap._mesh = self.model.mesh.get_surface(cap._mesh.id)
                if remove_duplicates:
                    node_ids = np.setdiff1d(cap.global_node_ids_edge, used_node_ids)
                else:
                    node_ids = cap.global_node_ids_edge

                if len(node_ids) == 0:
                    LOGGER.debug(
                        "Nodes already used. Skipping node set for {0}".format(
                            part.name + " " + cap.name
                        )
                    )
                    continue

                cap._node_set_id = self.get_unique_nodeset_id()

                kw = create_node_set_keyword(
                    node_ids + 1, node_set_id=cap._node_set_id, title=cap.name
                )
                self.kw_database.node_sets.append(kw)

                # node_set_id = node_set_id + 1

                used_node_ids = np.append(used_node_ids, node_ids)

        # add node-set for each surface
        for part in self.model.parts:
            for surface in part.surfaces:
                #! get up-to-date version of the surface.
                surface1 = self.model.mesh.get_surface(surface.id)
                if surface1.n_cells == 0:
                    LOGGER.debug(f"Failed to create node set for {surface.name}. Empty mesh.")
                    continue

                if remove_one_node_from_cell:
                    node_ids = self._filter_bc_nodes(surface1)
                else:
                    node_ids = surface1.global_node_ids_triangles
                if remove_duplicates:
                    node_ids = np.setdiff1d(node_ids, used_node_ids)

                surface._node_set_id = self.get_unique_nodeset_id()
                kw = create_node_set_keyword(
                    node_ids + 1, node_set_id=surface._node_set_id, title=surface.name
                )

                used_node_ids = np.append(used_node_ids, node_ids)

                self.kw_database.node_sets.append(kw)

    def _get_unique_id(self, keyword: str, return_used_ids: bool = False) -> int:
        """Get unique id of given keyword.

        Parameters
        ----------
        keyword : str
            Keyword string: valid inputs include:
            ["SECTION", "PART", "MAT", "SET_SEGMENT", "SET_NODE", "CURVE", ...]

        Returns
        -------
        int
            Next unique id
        """
        used_ids = [0]
        for key in self.kw_database.__dict__.keys():
            db = self.kw_database.__dict__[key]
            used_ids = np.append(used_ids, get_list_of_used_ids(db, keyword))
        used_ids = np.array(used_ids, dtype=int)
        _, counts = np.unique(used_ids, return_counts=True)
        if np.any(counts > 1):
            raise ValueError("{0} Duplicate ids found for: {1}".format(counts, keyword))

        if return_used_ids:
            return np.max(used_ids) + 1, used_ids
        else:
            return np.max(used_ids) + 1

    def get_unique_part_id(self) -> int:
        """Suggest a unique non-used part id."""
        return self._get_unique_id("PART")

    def get_unique_mat_id(self) -> int:
        """Suggest a unique non-used material id."""
        return self._get_unique_id("MAT")

    def get_unique_section_id(self) -> int:
        """Suggest a unique non-used section id."""
        return self._get_unique_id("SECTION")

    def get_unique_segmentset_id(self) -> int:
        """Suggest a unique non-used segment set id."""
        return self._get_unique_id("SET_SEGMENT")

    def get_unique_nodeset_id(self) -> int:
        """Suggest a unique non-used node set id."""
        return self._get_unique_id("SET_NODE")

    def get_unique_partset_id(self) -> int:
        """Suggest a unique non-used node set id."""
        return self._get_unique_id("SET_PART")

    def get_unique_curve_id(self) -> int:
        """Suggest a unique curve-id."""
        return self._get_unique_id("DEFINE_CURVE")

    def _get_list_of_includes(self):
        """Get a list of files to include in main.k. Omit any empty decks."""
        for deckname, deck in vars(self.kw_database).items():
            if deckname == "main":
                continue
            # skip if no keywords are present in the deck
            if len(deck.keywords) == 0:
                LOGGER.debug("No keywords in deck: {0}".format(deckname))
                continue
            self.include_files.append(deckname)
        return

    def _add_includes(self):
        """Add *INCLUDE keywords."""
        for include_file in self.include_files:
            filename_to_include = include_file + ".k"
            self.kw_database.main.append(keywords.Include(filename=filename_to_include))

        return

    def export(self, export_directory: str):
        """Write the model to files."""
        tstart = time.time()
        LOGGER.info("Writing all LS-DYNA .k files...")

        # is this reachable??
        if not export_directory:
            export_directory = os.path.join(
                self.model.workdir,
                self.__class__.__name__.lower().replace("dynawriter", ""),
            )

        if not os.path.isdir(export_directory):
            os.makedirs(export_directory)

        # export .k files
        self.export_databases(export_directory)

        # export settings
        self.settings.save(os.path.join(export_directory, "simulation_settings.yml"))

        tend = time.time()
        LOGGER.debug("Time spent writing files: {:.2f} s".format(tend - tstart))

        return

    def export_databases(self, export_directory: str):
        """Export each of non-empty databases to a specified directory."""
        if not export_directory:
            export_directory = self.model.info.working_directory

        for deckname, deck in vars(self.kw_database).items():
            # skip empty databases:
            if deck.keywords == []:
                continue
            LOGGER.info("Writing: {}".format(deckname))

            filepath = os.path.join(export_directory, deckname + ".k")

            if deckname == "solid_elements":
                if os.path.isfile(filepath):
                    os.remove(filepath)
                for element_kw in deck.keywords:
                    fast_element_writer(element_kw, filepath)
                with open(filepath, "a") as f:
                    f.write("*END\n")

            else:
                deck.export_file(filepath)

        return

    def _keep_ventricles(self):
        """Remove any non-ventricular parts."""
        LOGGER.debug("Just keeping ventricular-parts for fiber/purkinje generation")
        parts_to_keep = [
            p.name for p in self.model.parts if p.part_type in [PartType.VENTRICLE, PartType.SEPTUM]
        ]
        self._keep_parts(parts_to_keep)
        return

    def _keep_parts(self, parts_to_keep: List[str]):
        """Remove parts by a list of part names."""
        parts_to_remove = [part for part in self.model.part_names if part not in parts_to_keep]
        for part_to_remove in parts_to_remove:
            LOGGER.warning(f"Removing: {part_to_remove}")
            self.model.remove_part(part_to_remove)
        return

    def _update_solid_elements_db(self, add_fibers: bool = True):
        """
        Create Solid (ortho) elements for all parts.

        Parameters
        ----------
        add_fibers: bool, True
            if add fiber in general.
        """
        LOGGER.debug("Updating solid element keywords...")

        if add_fibers:
            cell_data_fields = self.model.mesh.cell_data.keys()
            if "fiber" not in cell_data_fields or "sheet" not in cell_data_fields:
                raise KeyError("Mechanics writer requires fiber and sheet fields")

        # create elements for each part
        for part in self.model.parts:
            if add_fibers and part.fiber:
                part_add_fibers = True
            else:
                part_add_fibers = False

            LOGGER.debug(
                "\tAdding elements for {0} | adding fibers: {1}".format(part.name, part_add_fibers)
            )
            #! This only works since tetrahedrons are at start of model.mesh, and surface
            #! cells are added behind these tetrahedrons.
            tetrahedrons = self.model.mesh.tetrahedrons[part.element_ids, :] + 1
            num_elements = tetrahedrons.shape[0]

            # element_ids = np.arange(1, num_elements + 1, 1) + solid_element_count
            part_ids = np.ones(num_elements, dtype=int) * part.pid

            # format the element keywords
            if not part_add_fibers:
                kw_elements = keywords.ElementSolid()
                elements = pd.DataFrame(
                    {
                        "eid": part.element_ids + 1,
                        "pid": part_ids,
                        "n1": tetrahedrons[:, 0],
                        "n2": tetrahedrons[:, 1],
                        "n3": tetrahedrons[:, 2],
                        "n4": tetrahedrons[:, 3],
                        "n5": tetrahedrons[:, 3],
                        "n6": tetrahedrons[:, 3],
                        "n7": tetrahedrons[:, 3],
                        "n8": tetrahedrons[:, 3],
                    }
                )
                kw_elements.elements = elements

            elif part_add_fibers:
                fiber = self.volume_mesh.cell_data["fiber"][part.element_ids]
                sheet = self.volume_mesh.cell_data["sheet"][part.element_ids]

                # normalize fiber and sheet directions:
                # norm = np.linalg.norm(fiber, axis=1)
                # fiber = fiber / norm[:, None]
                # norm = np.linalg.norm(sheet, axis=1)
                # sheet = sheet / norm[:, None]

                kw_elements = create_element_solid_ortho_keyword(
                    elements=tetrahedrons,
                    a_vec=fiber,
                    d_vec=sheet,
                    e_id=part.element_ids + 1,
                    part_id=part_ids,
                    element_type="tetra",
                )

            # add elements to database
            self.kw_database.solid_elements.append(kw_elements)
            # solid_element_count = solid_element_count + num_elements

        return


class MechanicsDynaWriter(BaseDynaWriter):
    """Class for preparing the input for a mechanics LS-DYNA simulation."""

    def __init__(
        self,
        model: HeartModel,
        settings: SimulationSettings = None,
    ) -> None:
        super().__init__(model=model, settings=settings)

        self.kw_database = MechanicsDecks()
        """Collection of keyword decks relevant for mechanics."""

        self.system_model_name = self.settings.mechanics.system.name
        """Name of system model to use."""

        self.set_flow_area: bool = True
        """If flow area is set for control volume."""
        return

    @property
    def system_model_name(self):
        """System model name.

        Notes
        -----
        Valid options include:
        ["ConstantPreloadWindkesselAfterload",
        "ClosedLoop].

        """
        return self._system_model

    @system_model_name.setter
    def system_model_name(self, value: str):
        if value not in [
            "ConstantPreloadWindkesselAfterload",
            "ClosedLoop",
        ]:
            raise ValueError("System model not valid")
        self._system_model = value

    def update(self, dynain_name: str = None, robin_bcs: list[Callable] = None):
        """Update the keyword database.

        Parameters
        ----------
        dynain_name : str, optional
            dynain file from stress free configuration computation, by default None
        robin_bcs : list[Callable], optional
            A list of lambda functions to apply Robin-type BCs, by default None

        Notes
        -----
        Do not need to write mesh files if dynain file is given.
        """
        self._update_main_db()

        self._add_damping()

        self._update_parts_db()
        self._update_material_db(add_active=True)
        self._update_segmentsets_db(add_caps=True)
        self._update_nodesets_db()

        if dynain_name is None:
            # write mesh
            self._update_node_db()
            self._update_solid_elements_db(add_fibers=True)
            # write cap shells with mesh
            self._update_cap_elements_db(add_mesh=True)
        else:
            self.kw_database.main.append(keywords.Include(filename=dynain_name))
            # cap mesh has been defined in dynain file
            self._update_cap_elements_db(add_mesh=False)

        # for boundary conditions
        if robin_bcs is None:
            # default BC
            self._add_cap_bc(bc_type=_BoundaryConditionType.ROBIN)
        else:
            # loop for every Robin BC function
            for robin_bc in robin_bcs:
                self.kw_database.boundary_conditions.extend(robin_bc())

        self._add_pericardium_bc()

        # for control volume
        system_settings = copy.deepcopy(self.settings.mechanics.system)
        system_settings._remove_units()

        if isinstance(self.model, LeftVentricle):
            lcid = self.get_unique_curve_id()
            system_map = [
                ControlVolume(
                    part=self.model.left_ventricle,
                    id=1,
                    Interactions=[
                        CVInteraction(
                            id=1,
                            cvid1=1,
                            cvid2=0,
                            lcid=lcid,
                            name="constant_preload_windkessel_afterload_left",
                            parameters=system_settings.left_ventricle,
                        )
                    ],
                )
            ]
        # Four chamber with active atrial
        elif isinstance(self, ElectroMechanicsDynaWriter) and isinstance(self.model, FourChamber):
            lcid = self.get_unique_curve_id()
            system_map = [
                ControlVolume(
                    part=self.model.left_ventricle,
                    id=1,
                    Interactions=[
                        CVInteraction(
                            id=1,
                            cvid1=1,
                            cvid2=0,
                            lcid=lcid,
                            name="afterload_windkessel_left",
                            parameters=system_settings.left_ventricle,
                        ),
                    ],
                ),
                ControlVolume(
                    part=self.model.right_ventricle,
                    id=2,
                    Interactions=[
                        CVInteraction(
                            id=2,
                            cvid1=2,
                            cvid2=0,
                            lcid=lcid + 1,
                            name="afterload_windkessel_right",
                            parameters=system_settings.right_ventricle,
                        ),
                    ],
                ),
                ControlVolume(
                    part=self.model.left_atrium,
                    id=3,
                    Interactions=[
                        CVInteraction(
                            id=3,
                            cvid1=3,
                            cvid2=0,
                            lcid=lcid + 2,
                            name="constant_flow_left_atrium",
                            parameters={"flow": -83.0},  # ~5 L/min
                        ),
                        CVInteraction(
                            id=4,
                            cvid1=3,
                            cvid2=1,
                            lcid=lcid + 3,
                            name="valve_mitral",
                            parameters={"Rv": 1e-6},
                        ),
                    ],
                ),
                ControlVolume(
                    part=self.model.right_atrium,
                    id=4,
                    Interactions=[
                        CVInteraction(
                            id=5,
                            cvid1=4,
                            cvid2=0,
                            lcid=lcid + 4,
                            name="constant_flow_right_atrium",
                            parameters={"flow": -83.0},  # ~5 L/min
                        ),
                        CVInteraction(
                            id=6,
                            cvid1=4,
                            cvid2=2,
                            lcid=lcid + 5,
                            name="valve_tricuspid",
                            parameters={"Rv": 1e-6},
                        ),
                    ],
                ),
            ]
        else:  # BiVentricle model or higher
            lcid = self.get_unique_curve_id()
            system_map = [
                ControlVolume(
                    part=self.model.left_ventricle,
                    id=1,
                    Interactions=[
                        CVInteraction(
                            id=1,
                            cvid1=1,
                            cvid2=0,
                            lcid=lcid,
                            name="constant_preload_windkessel_afterload_left",
                            parameters=system_settings.left_ventricle,
                        )
                    ],
                ),
                ControlVolume(
                    part=self.model.right_ventricle,
                    id=2,
                    Interactions=[
                        CVInteraction(
                            id=2,
                            cvid1=2,
                            cvid2=0,
                            lcid=lcid + 1,
                            name="constant_preload_windkessel_afterload_right",
                            parameters=system_settings.right_ventricle,
                        )
                    ],
                ),
            ]

        self._update_controlvolume_db(system_map)

        self._get_list_of_includes()
        self._add_includes()

        return

    def export(self, export_directory: str):
        """Write the model to files."""
        super().export(export_directory)

        # TODO: Close loop is only available from a customized LSDYNA executable
        # add system json in case of closed loop. For open-loop this is already
        # added in the control volume database
        if (
            self.system_model_name == "ClosedLoop"
            and self.__class__.__name__ == "MechanicsDynaWriter"
        ):
            # exports system model
            path_system_model_settings = os.path.join(
                export_directory, "system_model_settings.json"
            )
            with open(path_system_model_settings, "w") as outfile:
                json.dump(self.system_model_json, indent=4, fp=outfile)

        return

    def _update_main_db(self):
        """Update the main .k file."""
        LOGGER.debug("Updating main keywords...")

        self.kw_database.main.append("$$- Unit system: g-mm-ms-N-MPa-mJ -$$")
        self.kw_database.main.title = self.model.__class__.__name__

        if isinstance(self, ZeroPressureMechanicsDynaWriter):
            settings = self.settings.stress_free
            self._add_solution_controls()
            self._add_export_controls(settings.analysis.dt_d3plot.m)

        elif isinstance(self, (MechanicsDynaWriter, ElectroMechanicsDynaWriter)):
            settings = self.settings.mechanics
            self._add_solution_controls(
                end_time=settings.analysis.end_time.m,
                dtmin=settings.analysis.dtmin.m,
                dtmax=settings.analysis.dtmax.m,
            )
            self._add_export_controls(
                dt_output_d3plot=settings.analysis.dt_d3plot.m,
                dt_output_icvout=settings.analysis.dt_icvout.m,
            )

        return

    def _add_solution_controls(
        self,
        end_time: float = 5000,
        dtmin: float = 1.0,
        dtmax: float = 10.0,
        simulation_type: str = "quasi-static",
    ):
        """Add solution controls, output controls and solver settings."""
        # add termination keywords
        self.kw_database.main.append(keywords.ControlTermination(endtim=end_time))

        # add implicit controls
        if simulation_type == "quasi-static":
            imass = 1
            gamma = 0.6
            beta = 0.38
        elif simulation_type == "static":
            imass = 0
            gamma = 0.5
            beta = 0.25
        else:
            raise ValueError(
                "Simulation type not recognized: Please choose either quasi-static or static"
            )

        # prefill_time = self.parameters["Material"]["Myocardium"]["Active"]["Prefill"]
        self.kw_database.main.append(
            keywords.ControlImplicitDynamics(
                imass=imass,
                gamma=gamma,
                beta=beta,
                # active dynamic process only after prefilling
                # tdybir=prefill_time,
            )
        )

        self.kw_database.main.append("$$ Disable auto step due 0D model $$")
        self.kw_database.main.append(
            keywords.ControlImplicitAuto(iauto=0, dtmin=dtmin, dtmax=dtmax)
        )

        # add general implicit controls
        self.kw_database.main.append(
            keywords.ControlImplicitGeneral(imflag=1, dt0=dtmax)
        )  # imflag=1 means implicit

        # add implicit solution controls

        self.kw_database.main.append(
            keywords.ControlImplicitSolution(
                # maxref=35,
                dctol=0.02,
                ectol=1e6,
                rctol=1e3,
                abstol=-1e-20,
                dnorm=1,
                # diverg=2,
                lstol=-0.9,
                lsmtd=5,
                # d3itctl=1,
                nlprint=3,
                nlnorm=4,
            )
        )

        # add implicit solver controls
        self.kw_database.main.append(custom_keywords.ControlImplicitSolver(autospc=2))

        self.kw_database.main.append(keywords.ControlAccuracy(osu=1, inn=4, iacc=1))
        return

    def _add_export_controls(self, dt_output_d3plot: float = 0.05, dt_output_icvout: float = 0.001):
        """Add solution controls to the main simulation.

        Parameters
        ----------
        dt_output_d3plot : float, optional
            Writes full D3PLOT results at this time-step spacing, by default 0.05
        dt_output_icvout : float, optional
            Writes control volume results at this time-step spacing, by default 0.001
        """
        # add output control
        self.kw_database.main.append(keywords.ControlOutput(npopt=1, neecho=1, ikedit=0, iflush=0))

        # add export controls
        self.kw_database.main.append(keywords.DatabaseIcvout(dt=dt_output_icvout, binary=2))
        self.kw_database.main.append(keywords.DatabaseAbstat(dt=dt_output_icvout, binary=2))

        self.kw_database.main.append(keywords.DatabaseGlstat(dt=0.1, binary=2))

        self.kw_database.main.append(keywords.DatabaseMatsum(dt=0.1, binary=2))

        # # frequency of full results
        # lcid = self.get_unique_curve_id()
        # time = [
        #     0,
        #     self.parameters["Material"]["Myocardium"]["Active"]["Prefill"] * 0.99,
        #     self.parameters["Material"]["Myocardium"]["Active"]["Prefill"],
        #     self.parameters["Time"]["End Time"],
        # ]
        # step = [10 * dt_output_d3plot, 10 * dt_output_d3plot, dt_output_d3plot, dt_output_d3plot]
        # kw_curve = create_define_curve_kw(
        #     x=time,
        #     y=step,
        #     curve_name="d3plot out control",
        #     curve_id=lcid,
        #     lcint=0,
        # )
        # self.kw_database.main.append(kw_curve)

        self.kw_database.main.append(
            keywords.DatabaseBinaryD3Plot(
                dt=dt_output_d3plot,
                # lcdt=lcid, ioopt=1
            )
        )

        self.kw_database.main.append(
            keywords.DatabaseExtentBinary(neiph=27, strflg=1, maxint=0, resplt=1)
        )

        return

    def _add_damping(self):
        """Add damping to the main file."""
        lcid_damp = self.get_unique_curve_id()
        # mass damping
        kw_damp = keywords.DampingGlobal(lcid=lcid_damp)

        kw_damp_curve = create_define_curve_kw(
            x=[0, 10e25],  # to create a constant curve
            y=self.settings.mechanics.analysis.global_damping.m * np.array([1, 1]),
            curve_name="global damping [ms^-1]",
            curve_id=lcid_damp,
            lcint=0,
        )
        self.kw_database.main.append(kw_damp)
        self.kw_database.main.append(kw_damp_curve)

        # stiff damping
        for part in self.model.parts:
            self.kw_database.main.append(f"$$ {part.name} stiffness damping [ms]")
            kw = keywords.DampingPartStiffness(
                pid=part.pid, coef=self.settings.mechanics.analysis.stiffness_damping.m
            )
            self.kw_database.main.append(kw)
        return

    def _update_material_db(self, add_active: bool = True, em_couple: bool = False):
        #
        for part in self.model.parts:
            if isinstance(part.meca_material, MechanicalMaterialModel.DummyMaterial):
                # assign material for part if it's empty
                LOGGER.info(f"Material of {part.name} will be assigned automatically.")
                if part.fiber:
                    part.meca_material = self.settings.get_mechanical_material(
                        required_type="anisotropic", ep_coupled=em_couple
                    )
                    # disable active module
                    if not part.active:
                        part.meca_material.active = None

                else:
                    part.meca_material = self.settings.get_mechanical_material(
                        required_type="isotropic"
                    )
        # write
        for part in self.model.parts:
            material = part.meca_material

            if isinstance(material, Mat295):
                # need to write ca2+ curve
                if add_active and not em_couple and material.active is not None:
                    x, y = material.active.ca2_curve.dyna_input

                    cid = self.get_unique_curve_id()
                    curve_kw = create_define_curve_kw(
                        x=x,
                        y=y,
                        curve_name=f"ca2+ of {part.name}",
                        curve_id=cid,
                        lcint=10000,
                    )
                    self.kw_database.material.append(curve_kw)
                    material.active.acid = cid

                material_kw = MaterialHGOMyocardium(
                    id=part.mid, mat=material, ignore_active=not add_active
                )

                self.kw_database.material.append(material_kw)

            elif isinstance(material, NeoHookean):
                material_kw = MaterialNeoHook(
                    mid=part.mid,
                    rho=material.rho,
                    c10=material.c10,
                    nu=material.nu,
                    kappa=material.kappa,
                )
                self.kw_database.material.append(material_kw)

    def _add_cap_bc(self, bc_type: _BoundaryConditionType):
        """Add boundary condition to the cap.

        Parameters
        ----------
        bc_type : BoundaryType
           Boundary condition type.

        """
        # create list of cap names where to add the spring b.c
        constraint_caps = self._get_contraint_caps()

        if bc_type == _BoundaryConditionType.FIX:
            for part in self.model.parts:
                for cap in part.caps:
                    if cap.type in constraint_caps:
                        kw_fix = keywords.BoundarySpcSet()
                        kw_fix.nsid = cap._node_set_id
                        kw_fix.dofx = 1
                        kw_fix.dofy = 1
                        kw_fix.dofz = 1

                        self.kw_database.boundary_conditions.append(kw_fix)

        # if bc type is springs -> add springs
        elif bc_type == _BoundaryConditionType.ROBIN:
            part_id = self.get_unique_part_id()
            section_id = self.get_unique_section_id()
            mat_id = self.get_unique_mat_id()

            # read spring settings
            bc_settings = self.settings.mechanics.boundary_conditions
            spring_stiffness = bc_settings.valve["stiffness"].m
            scale_factor_normal = bc_settings.valve["scale_factor"]["normal"]
            scale_factor_radial = bc_settings.valve["scale_factor"]["radial"]

            part_kw = keywords.Part()
            part_df = pd.DataFrame(
                {
                    "pid": [part_id],
                    "secid": [section_id],
                    "mid": [mat_id],
                    "heading": ["SupportSpring"],
                }
            )
            part_kw.parts = part_df
            section_kw = keywords.SectionDiscrete(secid=section_id, cdl=0, tdl=0)
            mat_kw = keywords.MatSpringElastic(mid=mat_id, k=spring_stiffness)

            self.kw_database.boundary_conditions.append(part_kw)
            self.kw_database.boundary_conditions.append(section_kw)
            self.kw_database.boundary_conditions.append(mat_kw)

            # add springs for each cap
            caps = [cap for part in self.model.parts for cap in part.caps]
            for cap in caps:
                if cap.type in constraint_caps:
                    self.kw_database.boundary_conditions.append(f"$$ spring at {cap.name}$$")
                    self._add_springs_cap_edge(
                        cap,
                        part_id,
                        scale_factor_normal,
                        scale_factor_radial,
                    )

        return

    def _get_contraint_caps(self):
        """Get list of constraint caps depending on models."""
        constraint_caps = []

        if isinstance(self.model, LeftVentricle):
            constraint_caps = [CapType.MITRAL_VALVE, CapType.AORTIC_VALVE]

        elif isinstance(self.model, BiVentricle):
            constraint_caps = [
                CapType.MITRAL_VALVE,
                CapType.AORTIC_VALVE,
                CapType.TRICUSPID_VALVE,
                CapType.PULMONARY_VALVE,
            ]

        elif isinstance(self.model, (FourChamber, FullHeart)):
            constraint_caps = [
                CapType.SUPERIOR_VENA_CAVA,
                CapType.RIGHT_INFERIOR_PULMONARY_VEIN,
                CapType.RIGHT_SUPERIOR_PULMONARY_VEIN,
            ]

            if isinstance(self, ZeroPressureMechanicsDynaWriter):
                # add additional constraint to avoid rotation
                constraint_caps.extend([CapType.PULMONARY_VALVE])

        return constraint_caps

    def _add_springs_cap_edge(
        self,
        cap: Cap,
        part_id: int,
        scale_factor_normal: float,
        scale_factor_radial: float,
    ):
        """Add springs to the cap nodes.

        Notes
        -----
        Appends these to the boundary condition database.
        """
        LOGGER.debug(f"Adding spring b.c. for cap: {cap.name} of type {cap.type}")

        attached_nodes = cap.global_node_ids_edge

        # ? Can we compute this with only the cap mesh?
        #! This computes the nodal areas for all points in the cap mesh, including the central one.
        # compute nodal areas of nodes in cap elements.
        nodal_areas = compute_surface_nodal_area_pyvista(cap._mesh)[cap._local_node_ids_edge]

        # scaled spring stiffness by nodal area
        scale_factor_normal *= nodal_areas
        scale_factor_radial *= nodal_areas

        # add sd_orientiation, element discrete
        # compute the radial components
        # sd_orientations_radial = self.model.mesh.nodes[attached_nodes, :] - cap.centroid
        sd_orientations_radial = cap._mesh.nodes[cap._local_node_ids_edge] - cap.centroid

        # normalize
        norms = np.linalg.norm(sd_orientations_radial, axis=1)
        sd_orientations_radial = sd_orientations_radial / norms[:, None]

        # add sd direction normal to plane
        vector_id_normal = self.id_offset["vector"]
        sd_orientation_normal_kw = create_define_sd_orientation_kw(
            vectors=cap.cap_normal, vector_id_offset=vector_id_normal, iop=0
        )
        vector_id_normal += 1
        self.id_offset["vector"] += 1

        # add sd direction radial to nodes
        sd_orientation_radial_kw = create_define_sd_orientation_kw(
            vectors=sd_orientations_radial,
            vector_id_offset=self.id_offset["vector"],
            iop=0,
        )

        vector_ids_radial = sd_orientation_radial_kw.vectors["vid"].to_numpy()
        self.id_offset["vector"] = vector_ids_radial[-1]

        ## create discrete elements
        nodes_discrete_elements = np.array(
            [attached_nodes + 1, np.zeros(len(attached_nodes))], dtype=int
        ).T
        vector_ids_normal = np.ones(len(attached_nodes), dtype=int) * vector_id_normal

        #  for normal direction
        discrete_element_normal_kw = create_discrete_elements_kw(
            nodes=nodes_discrete_elements,
            part_id=part_id,
            vector_ids=vector_ids_normal,
            scale_factor=scale_factor_normal,
            element_id_offset=self.id_offset["element"]["discrete"],
        )

        self.id_offset["element"]["discrete"] = discrete_element_normal_kw.elements[
            "eid"
        ].to_numpy()[-1]

        #  for radial direction
        discrete_element_radial_kw = create_discrete_elements_kw(
            nodes=nodes_discrete_elements,
            part_id=part_id,
            vector_ids=vector_ids_radial,
            scale_factor=scale_factor_radial,
            element_id_offset=self.id_offset["element"]["discrete"],
        )

        self.id_offset["element"]["discrete"] = discrete_element_radial_kw.elements[
            "eid"
        ].to_numpy()[-1]

        # append to the database
        self.kw_database.boundary_conditions.append(sd_orientation_normal_kw)
        self.kw_database.boundary_conditions.append(sd_orientation_radial_kw)

        self.kw_database.boundary_conditions.append(discrete_element_normal_kw)
        self.kw_database.boundary_conditions.append(discrete_element_radial_kw)

        return

    def _add_pericardium_bc(self, scale=1.0):
        """Add the pericardium."""
        boundary_conditions = copy.deepcopy(self.settings.mechanics.boundary_conditions)
        robin_settings = boundary_conditions.robin

        # collect all pericardium nodes:
        ventricles_epi = self._get_epi_surface(apply=PartType.VENTRICLE)

        #! penalty function is defined on all nodes in the mesh: but just need the epicardial nodes.
        # penalty function
        penalty_function = self._get_longitudinal_penalty(robin_settings["ventricle"])

        ventricles_epi["scale factor"] = penalty_function[
            ventricles_epi.point_data["_global-point-ids"]
        ]
        # remove nodes with scale factor = 0
        ventricles_epi_reduce = ventricles_epi.threshold(
            value=[0.0001, 1], scalars="scale factor"
        ).extract_geometry()  # keep as polydata

        k = scale * robin_settings["ventricle"]["stiffness"].to("MPa/mm").m
        self.kw_database.pericardium.extend(
            self.write_robin_bc("spring", k, ventricles_epi_reduce, normal=None)
        )

        # damper
        dc = robin_settings["ventricle"]["damper"].to("MPa/mm*ms").m
        ventricles_epi.point_data.remove("scale factor")  # remove scale factor for spring
        self.kw_database.pericardium.extend(
            self.write_robin_bc("damper", dc, ventricles_epi, normal=None)
        )

        if isinstance(self.model, FourChamber):
            atrial_epi = self._get_epi_surface(PartType.ATRIUM)

            k = robin_settings["atrial"]["stiffness"].to("MPa/mm").m
            self.kw_database.pericardium.extend(
                self.write_robin_bc("spring", k, atrial_epi, normal=None)
            )

            dc = robin_settings["atrial"]["damper"].to("MPa/mm*ms").m
            self.kw_database.pericardium.extend(
                self.write_robin_bc("damper", dc, atrial_epi, normal=None)
            )
        return

    def _get_epi_surface(
        self, apply: Literal[PartType.VENTRICLE, PartType.ATRIUM] = PartType.VENTRICLE
    ):
        """Get the epicardial surfaces of either the ventricle or atria."""
        LOGGER.debug(f"Collecting epicardium nodesets of {apply}:")

        targets = [part for part in self.model.parts if apply == part.part_type]

        # retrieve combined epicardial surface from the central mesh object:
        # this ensures that we can use the global-point-ids
        epicardium_surface_ids = []
        for part in targets:
            try:
                epicardium_surface_ids.append(part.epicardium.id)
            except AttributeError:
                LOGGER.warning(f"{part.name} has no epicardium surface.")
                # part as "Atrioventricular isolation" may not have epicardium surface
                continue

        epicardium_surface1 = self.model.mesh.get_surface(epicardium_surface_ids)

        return epicardium_surface1

    def _get_longitudinal_penalty(self, pericardium_settings):
        """
        Use the universal ventricular longitudinal coordinate and a sigmoid penalty function.

        Strocchi et al 2020 doi: 10.1016/j.jbiomech.2020.109645.
        """
        penalty_c0 = pericardium_settings["penalty_function"][0]
        penalty_c1 = pericardium_settings["penalty_function"][1]
        self.kw_database.pericardium.append(f"$$ penalty with {penalty_c0}, {penalty_c1} $$")

        def _sigmoid(z):
            """Sigmoid function to scale spring coefficient."""
            return 1 / (1 + np.exp(-z))

        # compute penalty function from longitudinal coordinate
        try:
            uvc_l = self.model.mesh.point_data["apico-basal"]
        except KeyError:
            LOGGER.warning(
                "No apico-basal is found in point data, pericardium spring won't be created."
            )
            uvc_l = np.ones(self.model.mesh.GetNumberOfPoints())
        if np.any(uvc_l < 0):
            LOGGER.warning(
                "Negative normalized longitudinal coordinate detected."
                "Changing {0} negative uvc_l values to 1".format(np.sum((uvc_l < 0))),
            )
        uvc_l[uvc_l < 0] = 1

        penalty_function = -_sigmoid((abs(uvc_l) - penalty_c0) * penalty_c1) + 1
        return penalty_function

    def write_robin_bc(
        self,
        robin_type: Literal["spring", "damper"],
        constant: float,
        surface: pv.PolyData,
        normal: np.ndarray = None,
    ) -> list:
        """Create Robin BC on given surface.

        Parameters
        ----------
        robin_type : Literal[&quot;spring&quot;, &quot;damper&quot;]
            Create spring or damper
        constant : float
            stiffness (MPa/mm) or viscosity (MPa/mm*ms)
        surface : pv.PolyData
            Surface to apply BC, must contain point data '_global-point-ids'.
            Will be scaled by nodal area and point data 'scale factor' if exists
        normal : np.ndarray, optional
            If no normal given, use nodal normals, by default None

        Returns
        -------
        list
            list of dyna input deck
        """
        if surface.n_points == 0:
            LOGGER.error("Surface is empty, no Robin BC is added.")
            return []

        if "_global-point-ids" not in surface.point_data:
            raise ValueError("surface must contain pointdata '_global-point-ids'.")

        # global node ids where to apply the BC
        # NOTE: if we pass in a SurfaceMesh object we could use the
        # .global_node_ids attribute instead.
        nodes = surface["_global-point-ids"]

        # scale factor is nodal area
        # Add area flag in case pyvista defaults change.
        surf2 = surface.compute_cell_sizes(length=False, volume=False, area=True)
        scale_factor = np.array(
            surf2.cell_data_to_point_data().point_data["Area"].copy(), dtype=np.float32
        )
        if "scale factor" in surface.point_data:
            scale_factor *= np.array(surface.point_data["scale factor"], dtype=np.float32)

        # apply direction is nodal normal
        if normal is None:
            directions = surface.compute_normals().point_data["Normals"]
        elif normal.ndim == 1:
            directions = np.tile(normal, (len(nodes), 1))
        else:
            directions = normal

        # define spring orientations
        sd_orientation_kw = create_define_sd_orientation_kw(
            vectors=directions, vector_id_offset=self.id_offset["vector"]
        )
        vector_ids = sd_orientation_kw.vectors["vid"].to_numpy().astype(int)
        # update offset
        self.id_offset["vector"] = sd_orientation_kw.vectors["vid"].to_numpy()[-1]

        # create unique ids for keywords
        part_id = self.get_unique_part_id()
        section_id = self.get_unique_section_id()
        mat_id = self.get_unique_mat_id()

        # define material
        if robin_type == "spring":
            mat_kw = keywords.MatSpringElastic(mid=mat_id, k=constant)
        elif robin_type == "damper":
            mat_kw = keywords.MatDamperViscous(mid=mat_id, dc=constant)

        # define part
        part_kw = keywords.Part()
        part_kw.parts = pd.DataFrame(
            {
                "heading": [f"{robin_type}"],
                "pid": [part_id],
                "secid": [section_id],
                "mid": [mat_id],
            }
        )
        # define section
        section_kw = keywords.SectionDiscrete(secid=section_id, cdl=0, tdl=0)

        # 0: attached to ground
        n1_n2 = np.vstack([nodes + 1, np.zeros(len(nodes))]).T

        # create discrete elements
        discrete_element_kw = create_discrete_elements_kw(
            nodes=n1_n2,
            part_id=part_id,
            vector_ids=vector_ids,
            scale_factor=scale_factor,
            element_id_offset=self.id_offset["element"]["discrete"],
        )
        # add offset
        self.id_offset["element"]["discrete"] = discrete_element_kw.elements["eid"].to_numpy()[-1]

        # add keywords to database
        kw = []
        kw.append(part_kw)
        kw.append(section_kw)
        kw.append(mat_kw)
        kw.append(sd_orientation_kw)
        kw.append(discrete_element_kw)

        return kw

    def _update_cap_elements_db(self, add_mesh=True):
        """Update the database of shell elements.

        Notes
        -----
        Loops over all the defined caps/valves.
        """
        # material
        mat_null_id = self.get_unique_mat_id()
        material_kw = keywords.MatNull(
            mid=mat_null_id,
            ro=0.001,
        )

        # section
        section_id = self.get_unique_section_id()
        section_kw = keywords.SectionShell(
            secid=section_id,
            elform=4,
            shrf=0.8333,
            nip=3,
            t1=1,  # mm
        )

        self.kw_database.cap_elements.append(material_kw)
        self.kw_database.cap_elements.append(section_kw)

        caps = [cap for part in self.model.parts for cap in part.caps]
        # create new part for each cap
        cap_names_used = []
        for cap in caps:
            if cap.name in cap_names_used:
                # avoid to write mitral valve and triscupid valve twice
                LOGGER.debug("Already created material for {}: skipping".format(cap.name))
                continue

            cap.pid = self.get_unique_part_id()

            part_kw = keywords.Part()
            part_kw.parts = pd.DataFrame(
                {
                    "heading": [cap.name],
                    "pid": [cap.pid],
                    "secid": [section_id],
                    "mid": [mat_null_id],
                }
            )
            self.kw_database.cap_elements.append(part_kw)
            cap_names_used.append(cap.name)

            if cap.centroid is not None:
                if cap._node_set_id is None:
                    LOGGER.error("cap node set ID is not yet assigned")
                    exit()

                constraint = keywords.ConstrainedInterpolation(
                    icid=len(cap_names_used) + 1,
                    dnid=cap.global_centroid_id + 1,
                    ddof=123,
                    ityp=1,
                    fgm=0,
                    inid=cap._node_set_id,
                    idof=123,
                )
                self.kw_database.cap_elements.append(constraint)

        # create closing triangles for each cap
        # Note: cap parts already defined in control volume flow area, no mandatory here
        if add_mesh:
            # assumes there are no shells written yet since offset = 0
            # ? Should we use the global cell-index from self.mesh? or start from 0?
            shell_id_offset = 0
            cap_names_used = []
            for cap in caps:
                if cap.name in cap_names_used:
                    continue

                cap_mesh = self.model.mesh.get_surface(cap._mesh.id)

                shell_kw = create_element_shell_keyword(
                    shells=cap_mesh.triangles_global + 1,
                    part_id=cap.pid,
                    id_offset=shell_id_offset,
                )

                self.kw_database.cap_elements.append(shell_kw)

                shell_id_offset = shell_id_offset + cap_mesh.triangles_global.shape[0]
                cap_names_used.append(cap.name)
        return

    def _update_controlvolume_db(self, system_map: list[ControlVolume]):
        """Prepare the keywords for the control volume feature.

        Parameters
        ----------
        system_map : list[ControlVolume]
            list of control volume
        """
        if not self.system_model_name == "ConstantPreloadWindkesselAfterload":
            exit()

        def _create_null_part():
            # material
            mat_id = self.get_unique_mat_id()
            material_kw = keywords.MatNull(
                mid=mat_id,
                ro=0.001,
            )
            # section
            section_id = self.get_unique_section_id()
            section_kw = keywords.SectionShell(
                secid=section_id,
                elform=4,
                shrf=0.8333,
                nip=3,
                t1=1,
            )
            # part
            p_id = self.get_unique_part_id()
            part_kw = keywords.Part()
            part_kw.parts = pd.DataFrame(
                {
                    "heading": ["null flow area"],
                    "pid": [p_id],
                    "secid": [section_id],
                    "mid": [mat_id],
                }
            )

            self.kw_database.control_volume.append(section_kw)
            self.kw_database.control_volume.append(material_kw)
            self.kw_database.control_volume.append(part_kw)

            return p_id

        # create a new null part used in defining flow area
        if self.set_flow_area:
            pid = _create_null_part()

        for control_volume in system_map:
            part = control_volume.part
            cavity = part.cavity

            # DEFINE_CONTROL_VOLUME
            cv_kw = keywords.DefineControlVolume()
            cv_kw.id = control_volume.id
            cv_kw.sid = cavity.surface._seg_set_id
            self.kw_database.control_volume.append(cv_kw)

            if self.set_flow_area:
                # DEFINE_CONTROL_VOLUME_FLOW_AREA
                sid = self.get_unique_segmentset_id()
                sets = []
                for cap in part.caps:
                    sets.append(cap._seg_set_id)
                if len(sets) % 8 == 0:  # dynalib bug when length is 8,16,...
                    sets.append(0)
                self.kw_database.control_volume.append(keywords.SetSegmentAdd(sid=sid, sets=sets))

                # TODO: use dynalib: keywords.DefineControlVolumeFlowArea()
                flow_area_kw = "*DEFINE_CONTROL_VOLUME_FLOW_AREA\n"
                flow_area_kw += "$#    FAID     FCIID     FASID   FASTYPE       PID\n"
                flow_area_kw += "{0:10d}".format(control_volume.id)  # same as CVID
                flow_area_kw += "{0:10d}".format(control_volume.Interactions[0].id)  # first CVI id
                flow_area_kw += "{0:10d}".format(sid)
                flow_area_kw += "{0:10d}".format(2)  # flow area is defined by segment
                flow_area_kw += "{0:10d}".format(pid)
                self.kw_database.control_volume.append(flow_area_kw)

            for interaction in control_volume.Interactions:
                # DEFINE_CONTROL_VOLUME_INTERACTION
                cvi_kw = keywords.DefineControlVolumeInteraction()
                cvi_kw.id = interaction.id
                cvi_kw.cvid1 = interaction.cvid1
                cvi_kw.cvid2 = interaction.cvid2
                cvi_kw.lcid_ = interaction.lcid
                self.kw_database.control_volume.append(cvi_kw)

                # DEFINE FUNCTION
                define_function_wk = _define_function_0d_system(
                    function_id=interaction.lcid,
                    function_name=interaction.name,
                    parameters=interaction.parameters,
                )
                self.kw_database.control_volume.append(define_function_wk)

        return


class ZeroPressureMechanicsDynaWriter(MechanicsDynaWriter):
    """
    Class for preparing the input for a stress-free LS-DYNA simulation.

    Notes
    -----
    Derived from MechanicsDynaWriter and consequently derives all keywords relevant
    for simulations involving mechanics. This class does not write the
    control volume keywords but adds the keyword for computing the stress
    free configuration based on left/right cavity pressures instead.

    """

    def __init__(
        self,
        model: HeartModel,
        settings: SimulationSettings = None,
    ) -> None:
        super().__init__(model=model, settings=settings)

        self.kw_database = MechanicsDecks()
        """Collection of keyword decks relevant for mechanics."""

        return

    def update(self, robin_bcs: list[Callable] = None):
        """Update the keyword database.

        Parameters
        ----------
        robin_bcs : list[Callable], optional
            A list of lambda functions to apply Robin-type BCs, by default None
        """
        # bc_settings = self.settings.mechanics.boundary_conditions

        self._update_main_db()

        self.kw_database.main.title = self.model.__class__.__name__ + " zero-pressure"

        self._update_node_db()
        self._update_parts_db()
        self._update_solid_elements_db(add_fibers=True)
        self._update_segmentsets_db(add_caps=True)
        self._update_nodesets_db()
        self._update_material_db(add_active=False)
        self._update_cap_elements_db()

        # for boundary conditions
        if robin_bcs is None:
            # default BC
            self._add_cap_bc(bc_type=_BoundaryConditionType.FIX)
        else:
            # loop for every Robin BC function
            for robin_bc in robin_bcs:
                self.kw_database.boundary_conditions.extend(robin_bc())

        # Approximate end-diastolic pressures
        self._add_enddiastolic_pressure_bc()

        # zerop key words
        self._add_control_reference_configuration()

        # export dynain file
        save_part_ids = []
        for part in self.model.parts:
            save_part_ids.append(part.pid)

        caps = [cap for part in self.model.parts for cap in part.caps]
        for cap in caps:
            if cap.pid is not None:  # MV,TV for atrial parts get None
                save_part_ids.append(cap.pid)

        partset_id = self.get_unique_partset_id()
        kw = keywords.SetPartList(sid=partset_id)
        # kw.parts._data = save_part_ids
        # NOTE: when len(save_part_ids) = 8/16, dynalib bugs
        str = "\n"
        for i, id in enumerate(save_part_ids):
            str += "{0:10d}".format(id)
            if (i + 1) % 8 == 0:
                str += "\n"
        kw = kw.write() + str

        self.kw_database.main.append(kw)

        self.kw_database.main.append(
            custom_keywords.InterfaceSpringbackLsdyna(
                psid=partset_id,
                nshv=999,
                ftype=3,
                rflag=1,
                optc="OPTCARD",
                ndflag=1,
                cflag=1,
                hflag=1,
            )
        )

        self.kw_database.main.append(
            keywords.InterfaceSpringbackExclude(kwdname="BOUNDARY_SPC_NODE")
        )

        self._get_list_of_includes()
        self._add_includes()

        return

    def _add_export_controls(self, dt_output_d3plot: float = 0.5):
        """Rewrite method for zerop export.

        Parameters
        ----------
        dt_output_d3plot : float, optional
            Writes full D3PLOT results at this time-step spacing, by default 0.5
        """
        # add output control
        self.kw_database.main.append(keywords.ControlOutput(npopt=1, neecho=1, ikedit=0, iflush=0))

        # add export controls
        # self.kw_database.main.append(keywords.DatabaseElout(dt=0.1, binary=2))
        # self.kw_database.main.append(keywords.DatabaseGlstat(dt=0.1, binary=2))
        # self.kw_database.main.append(keywords.DatabaseMatsum(dt=0.1, binary=2))

        # frequency of full results
        self.kw_database.main.append(keywords.DatabaseBinaryD3Plot(dt=dt_output_d3plot))

        # self.kw_database.main.append(keywords.DatabaseExtentBinary(neiph=27, strflg=1, maxint=0))

        # add binout for post-process
        settings = copy.deepcopy(self.settings.stress_free)
        settings._remove_units()

        self.kw_database.main.append(
            keywords.DatabaseNodout(dt=settings.analysis.dt_nodout, binary=2)
        )

        # write for all nodes in nodout
        nodeset_id = self.get_unique_nodeset_id()
        kw = keywords.SetNodeGeneral(option="ALL", sid=nodeset_id)
        self.kw_database.main.append(kw)

        kw = keywords.DatabaseHistoryNodeSet(id1=nodeset_id)
        self.kw_database.main.append(kw)

        return

    def _add_solution_controls(self):
        """Rewrite method for the zerop simulation."""
        settings = copy.deepcopy(self.settings.stress_free)
        settings._remove_units()

        self.kw_database.main.append(keywords.ControlTermination(endtim=settings.analysis.end_time))

        self.kw_database.main.append(keywords.ControlImplicitDynamics(imass=0))

        # add auto step controls
        self.kw_database.main.append(
            keywords.ControlImplicitAuto(
                iauto=1, dtmin=settings.analysis.dtmin, dtmax=settings.analysis.dtmax
            )
        )

        # add general implicit controls
        self.kw_database.main.append(
            keywords.ControlImplicitGeneral(imflag=1, dt0=settings.analysis.dtmax)
        )

        # add implicit solution controls
        self.kw_database.main.append(
            keywords.ControlImplicitSolution(
                # maxref=35,
                dctol=0.01,
                ectol=1e6,
                rctol=1e3,
                abstol=1e-20,
                dnorm=1,
                diverg=2,
                # lsmtd=5,
            )
        )

        # add implicit solver controls
        self.kw_database.main.append(custom_keywords.ControlImplicitSolver(autospc=2))

        # accuracy control
        self.kw_database.main.append(keywords.ControlAccuracy(osu=1, inn=4, iacc=1))

        return

    def _add_control_reference_configuration(self):
        """Add control reference configuration keyword to main."""
        LOGGER.debug("Adding *CONTROL_REFERENCE_CONFIGURATION to main.k")
        settings = self.settings.stress_free.analysis
        kw = keywords.ControlReferenceConfiguration(
            maxiter=settings.max_iters,
            target="nodes.k",
            method=settings.method,
            tol=settings.tolerance,
        )

        self.kw_database.main.append(kw)

        return

    # def _add_enddiastolic_pressure_by_cv(self, pressure_lv: float = 1, pressure_rv: float = 1):
    #     """
    #     Apply end-of-diastolic pressure by control volume.

    #     Notes
    #     -----
    #     LSDYNA stress reference configuration lead to a bug with this load,
    #     it seems due to define function, need to be investigated.
    #     """
    #     cavities = [part.cavity for part in self.model.parts if part.cavity]
    #     for cavity in cavities:
    #         if "atrium" in cavity.name:
    #             continue

    #         # create CV
    #         cv_kw = keywords.DefineControlVolume()
    #         cv_kw.id = cavity.surface.id
    #         cv_kw.sid = cavity.surface._seg_set_id
    #         self.kw_database.main.append(cv_kw)

    #         # define CV interaction
    #         cvi_kw = keywords.DefineControlVolumeInteraction()
    #         cvi_kw.id = cavity.surface.id
    #         cvi_kw.cvid1 = cavity.surface._seg_set_id
    #         cvi_kw.cvid2 = 0  # ambient

    #         if "Left ventricle" in cavity.name:
    #             cvi_kw.lcid_ = 10
    #             pressure = pressure_lv
    #         elif "Right ventricle" in cavity.name:
    #             cvi_kw.lcid_ = 11
    #             pressure = pressure_rv

    #         self.kw_database.main.append(cvi_kw)

    #         # define define function
    #         definefunction_str = _ed_load_template()
    #         self.kw_database.main.append(
    #             definefunction_str.format(
    #                 cvi_kw.lcid_,
    #                 "flow_" + cavity.name.replace(" ", "_"),
    #                 pressure,
    #                 -200,
    #             )
    #         )

    #     self.kw_database.main.append(keywords.DatabaseIcvout(dt=10, binary=2))
    #     return

    def _add_enddiastolic_pressure_bc(self):
        """Add end diastolic pressure boundary condition on the left and right endocardium."""
        bc_settings = self.settings.mechanics.boundary_conditions
        pressure_lv = bc_settings.end_diastolic_cavity_pressure["left_ventricle"].m
        pressure_rv = bc_settings.end_diastolic_cavity_pressure["right_ventricle"].m
        pressure_la = bc_settings.end_diastolic_cavity_pressure["left_atrial"].m
        pressure_ra = bc_settings.end_diastolic_cavity_pressure["right_atrial"].m

        # create unit load curve
        load_curve_id = self.get_unique_curve_id()
        load_curve_kw = create_define_curve_kw(
            [0, 1, 1.001], [0, 1.0, 1.0], "unit load curve", load_curve_id, 100
        )

        load_curve_kw.sfa = 1000

        # append unit curve to main.k
        self.kw_database.main.append(load_curve_kw)

        # create *LOAD_SEGMENT_SETS for each ventricular cavity
        cavities = [part.cavity for part in self.model.parts if part.cavity]
        for cavity in cavities:
            if "Left ventricle" in cavity.name:
                load = keywords.LoadSegmentSet(
                    ssid=cavity.surface._seg_set_id, lcid=load_curve_id, sf=pressure_lv
                )
                self.kw_database.main.append(load)
            elif "Right ventricle" in cavity.name:
                load = keywords.LoadSegmentSet(
                    ssid=cavity.surface._seg_set_id, lcid=load_curve_id, sf=pressure_rv
                )
                self.kw_database.main.append(load)
            elif "Left atrium" in cavity.name:
                load = keywords.LoadSegmentSet(
                    ssid=cavity.surface._seg_set_id, lcid=load_curve_id, sf=pressure_la
                )
                self.kw_database.main.append(load)
            elif "Right atrium" in cavity.name:
                load = keywords.LoadSegmentSet(
                    ssid=cavity.surface._seg_set_id, lcid=load_curve_id, sf=pressure_ra
                )
                self.kw_database.main.append(load)
            else:
                LOGGER.debug(f"No load added to {cavity.name}")
                continue

        return


class FiberGenerationDynaWriter(BaseDynaWriter):
    """Class for preparing the input for a fiber-generation LS-DYNA simulation."""

    def __init__(self, model: HeartModel, settings: SimulationSettings = None) -> None:
        super().__init__(model=model, settings=settings)
        self.kw_database = FiberGenerationDecks()
        """Collection of keywords relevant for fiber generation."""

    def update(self, rotation_angles=None):
        """Update keyword database for Fiber generation: overwrites the inherited function."""
        ##
        self._update_main_db()  # needs updating

        if isinstance(self.model, (FourChamber, FullHeart)):
            LOGGER.warning(
                "Atrium present in the model, these will be removed for ventricle fiber generation."
            )

            parts = [
                part
                for part in self.model.parts
                if part.part_type in [PartType.VENTRICLE, PartType.SEPTUM]
            ]
            #! Note that this only works when tetrahedrons are added at the beginning
            #! of the mesh (file)! E.g. check self.mesh.celltypes to make sure this is the case!
            tet_ids = np.empty((0), dtype=int)
            for part in parts:
                tet_ids = np.append(tet_ids, part.element_ids)
                tets = self.model.mesh.tetrahedrons[tet_ids, :]
            nids = np.unique(tets)

            #  only write nodes attached to ventricle parts
            self._update_node_db(ids=nids)

            # remove parts not belonged to ventricles
            self._keep_ventricles()

            # remove segment which contains atrial nodes
            self._remove_atrial_nodes_from_ventricles_surfaces()

        else:
            self._update_node_db()

        self._update_parts_db()
        self._update_solid_elements_db(add_fibers=False)
        self._update_material_db()

        self._update_segmentsets_db(add_cavities=False)
        self._update_nodesets_db(remove_one_node_from_cell=True)

        # # update ep settings
        self._update_ep_settings()

        if rotation_angles is None:
            # find default settings
            rotation_angles = self.settings.get_ventricle_fiber_rotation(method="LSDYNA")
        self._update_create_fibers(rotation_angles)

        self._get_list_of_includes()
        self._add_includes()

        return

    def _remove_atrial_nodes_from_ventricles_surfaces(self):
        """Remove nodes other than ventricular from ventricular surfaces."""
        parts = [
            part
            for part in self.model.parts
            if part.part_type in [PartType.VENTRICLE, PartType.SEPTUM]
        ]

        tet_ids = np.empty((0), dtype=int)
        for part in parts:
            tet_ids = np.append(tet_ids, part.element_ids)
            tets = self.model.mesh.tetrahedrons[tet_ids, :]
        nids = np.unique(tets)

        for part in parts:
            for surface in part.surfaces:
                nodes_to_remove = surface.node_ids_triangles[
                    np.isin(
                        surface.node_ids_triangles,
                        nids,
                        assume_unique=True,
                        invert=True,
                    )
                ]

                faces = surface.faces.reshape(-1, 4)
                faces_to_remove = np.any(np.isin(faces, nodes_to_remove), axis=1)
                surface.faces = faces[np.invert(faces_to_remove)].ravel()

        return

    def _update_material_db(self):
        """Add simple linear elastic and orthotropic EM material for each defined part."""
        # collect myocardium and septum parts
        ventricles = [part for part in self.model.parts if "ventricle" in part.name]
        if isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
            septum = self.model.get_part("Septum")
            parts = ventricles + [septum]
        else:
            parts = ventricles
        material_settings = self.settings.electrophysiology.material
        for part in parts:
            # element_ids = part.element_ids
            # em_mat_id = self.get_unique_mat_id()
            em_mat_id = part.mid  #! Needs to match material id used in update_parts_db
            self.kw_database.material.extend(
                [
                    keywords.MatElastic(mid=em_mat_id, ro=1e-6, e=1),
                    custom_keywords.EmMat003(
                        mid=em_mat_id,
                        mtype=2,
                        sigma11=material_settings.myocardium["sigma_fiber"].m,
                        sigma22=material_settings.myocardium["sigma_sheet"].m,
                        sigma33=material_settings.myocardium["sigma_sheet_normal"].m,
                        beta=material_settings.myocardium["beta"].m,
                        cm=material_settings.myocardium["cm"].m,
                        aopt=2.0,
                        a1=0,
                        a2=0,
                        a3=1,
                        d1=0,
                        d2=-1,
                        d3=0,
                    ),
                    custom_keywords.EmEpCellmodelTomek(mid=em_mat_id),
                ]
            )

    def _update_ep_settings(self):
        """Add the settings for the electrophysiology solver."""
        self.kw_database.ep_settings.append(
            keywords.EmControl(
                emsol=11, numls=4, macrodt=1, dimtype=None, nperio=None, ncylbem=None
            )
        )

        # use defaults
        self.kw_database.ep_settings.append(custom_keywords.EmControlEp())

        # max iter should be int
        self.kw_database.ep_settings.append(
            keywords.EmSolverFem(reltol=1e-6, maxite=int(1e4), precon=2)
        )

        self.kw_database.ep_settings.append(keywords.EmOutput(mats=1, matf=1, sols=1, solf=1))

        return

    # TODO: Refactor
    def _update_create_fibers(self, rotation_angles):
        """Update the keywords for fiber generation."""
        # collect relevant node and segment sets.
        # node set: apex, base
        # node set: endocardium, epicardium
        # NOTE: could be better if basal nodes are extracted in the preprocessor
        # since that would allow you to robustly extract these nodessets using the
        # input data
        # The below is relevant for all models.
        nodes_base = np.empty(0, dtype=int)
        node_sets_ids_endo = []  # relevant for both models
        node_sets_ids_epi = []  # relevant for both models
        node_set_ids_epi_and_rseptum = []  # only relevant for bv, 4c and full model

        # list of ventricular parts
        ventricles = [part for part in self.model.parts if part.part_type == PartType.VENTRICLE]
        septum = next(
            (part for part in self.model.parts if part.part_type == PartType.SEPTUM),
            None,
        )

        # collect node set ids (already generated previously)
        node_sets_ids_epi = [ventricle.epicardium._node_set_id for ventricle in ventricles]
        node_sets_ids_endo = []
        for ventricle in ventricles:
            for surface in ventricle.surfaces:
                if "endocardium" in surface.name:
                    surf = self.model.mesh.get_surface(surface.id)
                    if surf.n_cells == 0:
                        LOGGER.debug(
                            f"Failed to collect node-set id for {surface.name}. Empty mesh."
                        )
                        continue
                    node_sets_ids_endo.append(surface._node_set_id)

        node_set_id_lv_endo = self.model.get_part("Left ventricle").endocardium._node_set_id
        if isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
            surfaces = [surface for p in self.model.parts for surface in p.surfaces]
            for surface in surfaces:
                #! relies on order of surfaces. Could be tricky.
                if surface.name == "Right ventricle endocardium septum":
                    node_set_ids_epi_and_rseptum = node_sets_ids_epi + [surface._node_set_id]
                    break

        for part in self.model.parts:
            for cap in part.caps:
                nodes_base = np.append(nodes_base, cap.global_node_ids_edge)

        # apex id [0] endocardium, [1] epicardum
        apex_point = self.model.get_part("Left ventricle").apex_points[1]
        if "epicardium" not in apex_point.name:
            raise ValueError("Expecting a point on the epicardium")
        node_apex = apex_point.node_id  #! is this a global node id?

        # validate node set by removing nodes not part of the model without ventricles
        tet_ids_ventricles = np.empty((0), dtype=int)
        if septum:
            parts = ventricles + [septum]
        else:
            parts = ventricles

        for part in parts:
            tet_ids_ventricles = np.append(tet_ids_ventricles, part.element_ids)

        tetra_ventricles = self.model.mesh.tetrahedrons[tet_ids_ventricles, :]

        # remove nodes that occur just in atrial part
        mask = np.isin(nodes_base, tetra_ventricles, invert=True)
        LOGGER.debug("Removing {0} nodes from base nodes".format(np.sum(mask)))
        nodes_base = nodes_base[np.invert(mask)]

        # create set parts for lv and rv myocardium
        myocardium_part_ids = [ventricle.pid for ventricle in ventricles]

        # switch between the various models to generate valid input decks
        if isinstance(self.model, LeftVentricle):
            LOGGER.warning("Model type %s in development " % self.model.__class__.__name__)

            # Define part set for myocardium
            part_list1_kw = keywords.SetPartList(
                sid=1,
            )
            part_list1_kw.parts._data = myocardium_part_ids
            part_list1_kw.options["TITLE"].active = True
            part_list1_kw.title = "myocardium_all"

            self.kw_database.create_fiber.extend([part_list1_kw])

            # combine node sets endocardium uing *SET_NODE_ADD:
            node_set_id_all_endocardium = self.get_unique_nodeset_id()

            set_add_kw = keywords.SetNodeAdd(sid=node_set_id_all_endocardium)
            set_add_kw.options["TITLE"].active = True
            set_add_kw.title = "all_endocardium_segments"
            set_add_kw.nodes._data = node_sets_ids_endo

            self.kw_database.create_fiber.append(set_add_kw)

            # combine node sets epicardium:
            node_set_id_all_epicardium = self.get_unique_nodeset_id()
            set_add_kw = keywords.SetNodeAdd(sid=node_set_id_all_epicardium)
            set_add_kw.options["TITLE"].active = True
            set_add_kw.title = "all_epicardium_segments"
            set_add_kw.nodes._data = node_sets_ids_epi

            self.kw_database.create_fiber.append(set_add_kw)

            node_set_id_base = self.get_unique_nodeset_id()
            node_set_id_apex = self.get_unique_nodeset_id() + 1

            # create node-sets for base and apex
            node_set_base_kw = create_node_set_keyword(
                node_ids=nodes_base + 1,
                node_set_id=node_set_id_base,
                title="base nodes",
            )
            node_set_apex_kw = create_node_set_keyword(
                node_ids=node_apex + 1, node_set_id=node_set_id_apex, title="apex node"
            )

            self.kw_database.create_fiber.extend([node_set_base_kw, node_set_apex_kw])

            # Set up *EM_EP_FIBERINITIAL keyword
            # apex > base
            self.kw_database.create_fiber.append(
                custom_keywords.EmEpFiberinitial(
                    id=1,
                    partid=1,  # set part id 1: myocardium
                    stype=2,  # set type 2 == nodes
                    ssid1=node_set_id_base,
                    ssid2=node_set_id_apex,
                )
            )

            # all epicardium > all endocardium
            self.kw_database.create_fiber.append(
                custom_keywords.EmEpFiberinitial(
                    id=2,
                    partid=1,  # set part id 1: myocardium
                    stype=2,  # set type 1 == segment set, set type 2 == node set
                    ssid1=node_set_id_all_epicardium,
                    ssid2=node_set_id_all_endocardium,
                )
            )

            # add *EM_EP_CREATEFIBERORIENTATION keywords
            self.kw_database.create_fiber.append(
                custom_keywords.EmEpCreatefiberorientation(
                    partsid=1,
                    solvid1=1,
                    solvid2=2,
                    alpha=-101,
                    beta=-102,
                    wfile=1,
                    prerun=1,
                )
            )

            # define functions:
            from ansys.heart.core.writer.define_function_templates import (
                _function_alpha,
                _function_beta,
                _function_beta_septum,
            )

            self.kw_database.create_fiber.append(
                keywords.DefineFunction(
                    fid=101,
                    function=_function_alpha(
                        alpha_endo=rotation_angles["alpha"][0],
                        alpha_epi=rotation_angles["alpha"][1],
                    ),
                )
            )
            self.kw_database.create_fiber.append(
                keywords.DefineFunction(
                    fid=102,
                    function=_function_beta(
                        beta_endo=rotation_angles["beta"][0],
                        beta_epi=rotation_angles["beta"][1],
                    ),
                )
            )

        elif isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
            septum_part_ids = [self.model.get_part("Septum").pid]

            # Define part set for myocardium
            part_list1_kw = keywords.SetPartList(
                sid=1,
            )
            part_list1_kw.parts._data = myocardium_part_ids
            part_list1_kw.options["TITLE"].active = True
            part_list1_kw.title = "myocardium_all"

            # Define part set for septum
            part_list2_kw = keywords.SetPartList(
                sid=2,
            )
            part_list2_kw.options["TITLE"].active = True
            part_list2_kw.title = "septum"
            part_list2_kw.parts._data = septum_part_ids

            self.kw_database.create_fiber.extend([part_list1_kw, part_list2_kw])

            # combine node sets endocardium uing *SET_SEGMENT_ADD:
            node_set_id_all_endocardium = self.get_unique_nodeset_id()
            set_add_kw = keywords.SetNodeAdd(sid=node_set_id_all_endocardium)

            set_add_kw.options["TITLE"].active = True
            set_add_kw.title = "all_endocardium_segments"
            set_add_kw.nodes._data = node_sets_ids_endo

            self.kw_database.create_fiber.append(set_add_kw)

            # combine node sets epicardium:
            node_set_id_all_epicardium = self.get_unique_nodeset_id()
            set_add_kw = keywords.SetNodeAdd(sid=node_set_id_all_epicardium)

            set_add_kw.options["TITLE"].active = True
            set_add_kw.title = "all_epicardium_segments"
            set_add_kw.nodes._data = node_sets_ids_epi

            self.kw_database.create_fiber.append(set_add_kw)

            # combine node sets epicardium and septum:
            node_set_all_but_left_endocardium = self.get_unique_nodeset_id()
            set_add_kw = keywords.SetNodeAdd(sid=node_set_all_but_left_endocardium)

            set_add_kw.options["TITLE"].active = True
            set_add_kw.title = "all_but_left_endocardium"
            set_add_kw.nodes._data = node_set_ids_epi_and_rseptum

            self.kw_database.create_fiber.append(set_add_kw)

            node_set_id_base = self.get_unique_nodeset_id()
            node_set_id_apex = self.get_unique_nodeset_id() + 1
            # create node-sets for base and apex
            node_set_base_kw = create_node_set_keyword(
                node_ids=nodes_base + 1,
                node_set_id=node_set_id_base,
                title="base nodes",
            )
            node_set_apex_kw = create_node_set_keyword(
                node_ids=node_apex + 1, node_set_id=node_set_id_apex, title="apex node"
            )

            self.kw_database.create_fiber.extend([node_set_base_kw, node_set_apex_kw])

            # Set up *EM_EP_FIBERINITIAL keyword
            # apex > base
            self.kw_database.create_fiber.append(
                custom_keywords.EmEpFiberinitial(
                    id=1,
                    partid=1,  # set part id 1: myocardium
                    stype=2,  # set type 2 == nodes
                    ssid1=node_set_id_base,
                    ssid2=node_set_id_apex,
                )
            )

            # all epicardium > all endocardium
            self.kw_database.create_fiber.append(
                custom_keywords.EmEpFiberinitial(
                    id=2,
                    partid=1,  # set part id 1: myocardium
                    stype=2,  # set type 1 == segment set, set type 2 == node set
                    ssid1=node_set_id_all_epicardium,
                    ssid2=node_set_id_all_endocardium,
                )
            )

            # all epicardium > endocardium left ventricle
            self.kw_database.create_fiber.append(
                custom_keywords.EmEpFiberinitial(
                    id=3,
                    partid=2,  # set part id 2: septum
                    stype=2,  # set type 1 == segment set
                    ssid1=node_set_all_but_left_endocardium,
                    ssid2=node_set_id_lv_endo,
                )
            )

            # add *EM_EP_CREATEFIBERORIENTATION keywords
            self.kw_database.create_fiber.append(
                custom_keywords.EmEpCreatefiberorientation(
                    partsid=1,
                    solvid1=1,
                    solvid2=2,
                    alpha=-101,
                    beta=-102,
                    wfile=1,
                    prerun=1,
                )
            )
            # add *EM_EP_CREATEFIBERORIENTATION keywords
            self.kw_database.create_fiber.append(
                custom_keywords.EmEpCreatefiberorientation(
                    partsid=2,
                    solvid1=1,
                    solvid2=3,
                    alpha=-101,
                    beta=-103,
                    wfile=1,
                    prerun=1,
                )
            )

            # define functions:
            from ansys.heart.core.writer.define_function_templates import (
                _function_alpha,
                _function_beta,
                _function_beta_septum,
            )

            self.kw_database.create_fiber.append(
                keywords.DefineFunction(
                    fid=101,
                    function=_function_alpha(
                        alpha_endo=rotation_angles["alpha"][0],
                        alpha_epi=rotation_angles["alpha"][1],
                    ),
                )
            )
            self.kw_database.create_fiber.append(
                keywords.DefineFunction(
                    fid=102,
                    function=_function_beta(
                        beta_endo=rotation_angles["beta"][0],
                        beta_epi=rotation_angles["beta"][1],
                    ),
                )
            )
            self.kw_database.create_fiber.append(
                keywords.DefineFunction(
                    fid=103,
                    function=_function_beta_septum(
                        beta_endo=rotation_angles["beta_septum"][0],
                        beta_epi=rotation_angles["beta_septum"][1],
                    ),
                )
            )

    def _update_main_db(self):
        self.kw_database.main.append(
            keywords.ControlTimeStep(dtinit=1.0, dt2ms=1.0, emscl=None, ihdo=None, rmscl=None)
        )

        self.kw_database.main.append(keywords.ControlTermination(endtim=10))

        self.kw_database.main.append(keywords.DatabaseBinaryD3Plot(dt=1.0))

        return


class PurkinjeGenerationDynaWriter(BaseDynaWriter):
    """Class for preparing the input for a Purkinje LS-DYNA simulation."""

    def __init__(
        self,
        model: HeartModel,
        settings: SimulationSettings = None,
    ) -> None:
        super().__init__(model=model, settings=settings)
        self.kw_database = PurkinjeGenerationDecks()
        """Collection of keywords relevant for Purkinje generation."""

    def update(self):
        """Update keyword database - overwrites the inherited function."""
        ##
        self._update_main_db()  # needs updating

        self._update_node_db()  # can stay the same (could move to base class)
        if isinstance(self.model, (FourChamber, FullHeart)):
            LOGGER.warning(
                "Atrium present in the model, "
                "these will be removed for ventricle Purkinje generation."
            )
            self._keep_ventricles()

        self._update_parts_db()  # can stay the same (could move to base class++++++++++++++++++++)
        self._update_solid_elements_db(add_fibers=False)
        self._update_material_db()

        self._update_segmentsets_db(add_cavities=False)  # can stay the same
        self._update_nodesets_db()  # can stay the same

        # update ep settings
        self._update_ep_settings()
        self._update_create_Purkinje()

        self._get_list_of_includes()
        self._add_includes()

        return

    def _update_material_db(self):
        """Add simple linear elastic material for each defined part."""
        material_settings = self.settings.electrophysiology.material
        for part in self.model.parts:
            em_mat_id = part.pid
            self.kw_database.material.extend(
                [
                    keywords.MatElastic(mid=em_mat_id, ro=1e-6, e=1),
                    custom_keywords.EmMat003(
                        mid=em_mat_id,
                        mtype=2,
                        sigma11=material_settings.myocardium["sigma_fiber"].m,
                        sigma22=material_settings.myocardium["sigma_sheet"].m,
                        sigma33=material_settings.myocardium["sigma_sheet_normal"].m,
                        beta=material_settings.myocardium["beta"].m,
                        cm=material_settings.myocardium["cm"].m,
                        aopt=2.0,
                        a1=0,
                        a2=0,
                        a3=1,
                        d1=0,
                        d2=-1,
                        d3=0,
                    ),
                ]
            )

    def _update_ep_settings(self):
        """Add the settings for the electrophysiology solver."""
        self.kw_database.ep_settings.append(
            keywords.EmControl(
                emsol=11, numls=4, macrodt=1, dimtype=None, nperio=None, ncylbem=None
            )
        )

        self.kw_database.ep_settings.append(keywords.EmOutput(mats=1, matf=1, sols=1, solf=1))

        return

    def _update_create_Purkinje(self):  # noqa N802
        """Update the keywords for Purkinje generation."""
        # collect relevant node and segment sets.
        # node set: apex, base
        # node set: endocardium, epicardium
        # NOTE: could be better if basal nodes are extracted in the preprocessor
        # since that would allow you to robustly extract these nodessets using the
        # input data
        # The below is relevant for all models.

        node_origin_left = np.empty(0, dtype=int)
        node_origin_right = np.empty(0, dtype=int)
        edge_id_start_left = np.empty(0, dtype=int)
        edge_id_start_right = np.empty(0, dtype=int)

        # apex_points[0]: endocardium, apex_points[1]: epicardium
        if isinstance(self.model, (LeftVentricle, BiVentricle, FourChamber, FullHeart)):
            if self.settings.purkinje.node_id_origin_left is None:
                node_origin_left = self.model.left_ventricle.apex_points[0].node_id
            else:
                node_origin_left = self.settings.purkinje.node_id_origin_left

            segment_set_ids_endo_left = self.model.left_ventricle.endocardium._seg_set_id

            # check whether point is on edge of endocardium - otherwise pick another node in
            # the same triangle
            #! Get an up-to-date version of the endocardium.
            endocardium = self.model.mesh.get_surface(self.model.left_ventricle.endocardium.id)
            #! Need to boundary edges to global ids.
            if np.any(
                endocardium.point_data["_global-point-ids"][endocardium.boundary_edges]
                == node_origin_left
            ):
                element_id = np.argwhere(
                    np.any(endocardium.triangles_global == node_origin_left, axis=1)
                )[0][0]

                node_origin_left = endocardium.triangles_global[element_id, :][
                    np.argwhere(
                        np.isin(
                            endocardium.triangles_global[element_id, :],
                            endocardium.point_data["_global-point-ids"][endocardium.boundary_edges],
                            invert=True,
                        )
                    )[0][0]
                ]
                LOGGER.debug(
                    "Node id {0} is on edge of {1}. Picking node id {2}".format(
                        self.model.left_ventricle.apex_points[0].node_id,
                        endocardium.name,
                        node_origin_left,
                    )
                )
                self.model.left_ventricle.apex_points[0].node_id = node_origin_left

            node_set_id_apex_left = self.get_unique_nodeset_id()
            # create node-sets for apex
            node_set_apex_kw = create_node_set_keyword(
                node_ids=[node_origin_left + 1],
                node_set_id=node_set_id_apex_left,
                title="apex node left",
            )

            self.kw_database.node_sets.append(node_set_apex_kw)

            apex_left_coordinates = self.model.mesh.points[node_origin_left, :]

            #! Is this to get unused start node/edge indinces?
            node_id_start_left = self.model.mesh.points.shape[0] + 1

            edge_id_start_left = self.model.mesh.tetrahedrons.shape[0] + 1

            pid = self.get_unique_part_id()
            # Purkinje generation parameters
            self.kw_database.main.append(
                custom_keywords.EmEpPurkinjeNetwork2(
                    purkid=1,
                    buildnet=1,
                    ssid=segment_set_ids_endo_left,
                    mid=pid,
                    pointstx=apex_left_coordinates[0],
                    pointsty=apex_left_coordinates[1],
                    pointstz=apex_left_coordinates[2],
                    edgelen=self.settings.purkinje.edgelen.m,
                    ngen=self.settings.purkinje.ngen.m,
                    nbrinit=self.settings.purkinje.nbrinit.m,
                    nsplit=self.settings.purkinje.nsplit.m,
                    inodeid=node_id_start_left,
                    iedgeid=edge_id_start_left,  # TODO: check if beam elements exist in mesh
                    pmjtype=self.settings.purkinje.pmjtype.m,
                    pmjradius=self.settings.purkinje.pmjradius.m,
                    pmjrestype=self.settings.electrophysiology.material.beam["pmjrestype"].m,
                    pmjres=self.settings.electrophysiology.material.beam["pmjres"].m,
                )
            )

        # Add right purkinje only in biventricular or 4chamber models
        if isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
            if self.settings.purkinje.node_id_origin_right is None:
                node_origin_right = self.model.right_ventricle.apex_points[0].node_id
            else:
                node_origin_right = self.settings.purkinje.node_id_origin_right

            segment_set_ids_endo_right = (
                self.model.right_ventricle.endocardium._seg_set_id
            )  # TODO: Replace

            # check whether point is on edge of endocardium - otherwise pick another node in
            # the same triangle
            #! Make sure endocardium is an updated version (e.g. point/cell data is up to date.)
            endocardium = self.model.mesh.get_surface(self.model.right_ventricle.endocardium.id)
            # endocardium.get_boundary_edges()
            if np.any(endocardium.boundary_edges_global == node_origin_right):
                element_id = np.argwhere(
                    np.any(endocardium.triangles_global == node_origin_right, axis=1)
                )[0][0]

                node_origin_right = endocardium.triangles_global[element_id, :][
                    np.argwhere(
                        np.isin(
                            endocardium.triangles_global[element_id, :],
                            endocardium.boundary_edges_global,
                            invert=True,
                        )
                    )[0][0]
                ]
                LOGGER.debug(
                    "Node id {0} is on edge of {1}. Picking node id {2}".format(
                        self.model.right_ventricle.apex_points[0].node_id,
                        endocardium.name,
                        node_origin_right,
                    )
                )
                self.model.right_ventricle.apex_points[0].node_id = node_origin_right

            node_set_id_apex_right = self.get_unique_nodeset_id()
            # create node-sets for apex
            node_set_apex_kw = create_node_set_keyword(
                node_ids=[node_origin_right + 1],
                node_set_id=node_set_id_apex_right,
                title="apex node right",
            )

            self.kw_database.node_sets.append(node_set_apex_kw)

            apex_right_coordinates = self.model.mesh.points[node_origin_right, :]

            node_id_start_right = (
                2 * self.model.mesh.points.shape[0]
            )  # TODO: find a solution in dyna to better handle id definition

            edge_id_start_right = 2 * self.model.mesh.tetrahedrons.shape[0]
            pid = self.get_unique_part_id() + 1
            # Purkinje generation parameters
            self.kw_database.main.append(
                custom_keywords.EmEpPurkinjeNetwork2(
                    purkid=2,
                    buildnet=1,
                    ssid=segment_set_ids_endo_right,
                    mid=pid,
                    pointstx=apex_right_coordinates[0],
                    pointsty=apex_right_coordinates[1],
                    pointstz=apex_right_coordinates[2],
                    edgelen=self.settings.purkinje.edgelen.m,
                    ngen=self.settings.purkinje.ngen.m,
                    nbrinit=self.settings.purkinje.nbrinit.m,
                    nsplit=self.settings.purkinje.nsplit.m,
                    inodeid=node_id_start_right,  # TODO: check if beam elements exist in mesh
                    iedgeid=edge_id_start_right,
                    pmjtype=self.settings.purkinje.pmjtype.m,
                    pmjradius=self.settings.purkinje.pmjradius.m,
                    pmjrestype=self.settings.electrophysiology.material.beam["pmjrestype"].m,
                    pmjres=self.settings.electrophysiology.material.beam["pmjres"].m,
                )
            )

    def _update_main_db(self):
        return


class ElectrophysiologyDynaWriter(BaseDynaWriter):
    """Class for preparing the input for an Electrophysiology LS-DYNA simulation."""

    def __init__(
        self,
        model: Union[HeartModel, FullHeart, FourChamber, BiVentricle, LeftVentricle],
        settings: SimulationSettings = None,
    ) -> None:
        if isinstance(model, FourChamber):
            model._create_atrioventricular_isolation()
        if model._add_blood_pool:
            model._create_blood_part()

        super().__init__(model=model, settings=settings)
        self.kw_database = ElectrophysiologyDecks()
        """Collection of keywords relevant for Electrophysiology."""

    def update(self):
        """Update keyword database for Electrophysiology."""
        # self._isolate_atria_and_ventricles()

        ##
        self._update_main_db()
        self._update_solution_controls()
        self._update_export_controls()

        self._update_node_db()
        self._update_parts_db()
        self._update_solid_elements_db(add_fibers=True)

        self._update_dummy_material_db()
        self._update_ep_material_db()

        self._update_segmentsets_db(add_cavities=True)

        # TODO: check if no existing node set ids conflict with surface ids
        # For now, new node sets should be created after calling
        # self._update_nodesets_db()
        self._update_nodesets_db()
        self._update_parts_cellmodels()

        if self.model.conduction_system.number_of_cells != 0:
            # with smcoupl=1, mechanical coupling is disabled
            # with thcoupl=1, thermal coupling is disabled
            self.kw_database.ep_settings.append(keywords.EmControlCoupling(thcoupl=1, smcoupl=1))
            self._update_use_Purkinje()

        # update ep settings
        self._update_ep_settings()
        self._update_stimulation()

        if self.model._add_blood_pool:
            self._update_blood_settings()

        if hasattr(self.model, "electrodes") and len(self.model.electrodes) != 0:
            self._update_ECG_coordinates()

        self._get_list_of_includes()
        self._add_includes()

        return

    def _update_dummy_material_db(self):
        """Add simple mechanics material for each defined part."""
        for part in self.model.parts:
            ep_mid = part.pid
            self.kw_database.material.append(
                keywords.MatElastic(mid=ep_mid, ro=1e-6, e=1),
            )

    def _update_ep_material_db(self):
        """Add EP material for each defined part."""
        material_settings = self.settings.electrophysiology.material
        solvertype = self.settings.electrophysiology.analysis.solvertype
        if solvertype == "Monodomain":
            sig1 = material_settings.myocardium["sigma_fiber"].m
            sig2 = material_settings.myocardium["sigma_sheet"].m
            sig3 = material_settings.myocardium["sigma_sheet_normal"].m
        elif solvertype == "Eikonal" or solvertype == "ReactionEikonal":
            sig1 = material_settings.myocardium["velocity_fiber"].m
            sig2 = material_settings.myocardium["velocity_sheet"].m
            sig3 = material_settings.myocardium["velocity_sheet_normal"].m

        for part in self.model.parts:
            if isinstance(part.ep_material, EPMaterial.DummyMaterial):
                LOGGER.info(f"Material of {part.name} will be assigned automatically.")
                if part.active:
                    part.ep_material = EPMaterial.Active(sigma_fiber=sig1)
                else:
                    part.ep_material = EPMaterial.Passive(sigma_fiber=sig1)
                if part.fiber:
                    part.ep_material.sigma_sheet = sig2
                    part.ep_material.sigma_sheet_normal = sig3

            self.kw_database.material.append(f"$$ {part.name} $$")
            ep_mid = part.pid
            kw = self._get_ep_material_kw(ep_mid, part.ep_material)
            self.kw_database.material.append(kw)

        return

    def _update_parts_cellmodels(self):
        """Add cell model for each defined part."""
        for part in self.model.parts:
            if isinstance(part.ep_material, EPMaterial.Active):
                ep_mid = part.pid
                # One cell model for myocardium, default value is epi layer parameters
                self._add_cell_model_keyword(matid=ep_mid, cellmodel=part.ep_material.cell_model)
        # different cell models for endo/mid/epi layer
        # TODO:  this will override previous definition?
        #        what's the situation at setptum? and at atrial?
        if "transmural" in self.model.mesh.point_data.keys():
            (
                endo_id,
                mid_id,
                epi_id,
            ) = self._create_myocardial_nodeset_layers()
            tentusscher_endo = CellModel.TentusscherEndo()
            tentusscher_mid = CellModel.TentusscherMid()
            tentusscher_epi = CellModel.TentusscherEpi()

            self._add_Tentusscher_keyword(matid=-endo_id, params=tentusscher_endo.to_dictionary())
            self._add_Tentusscher_keyword(matid=-mid_id, params=tentusscher_mid.to_dictionary())
            self._add_Tentusscher_keyword(matid=-epi_id, params=tentusscher_epi.to_dictionary())

    def _create_myocardial_nodeset_layers(self):
        percent_endo = self.settings.electrophysiology.material.myocardium["percent_endo"]
        percent_mid = self.settings.electrophysiology.material.myocardium["percent_mid"]
        values = self.model.mesh.point_data["transmural"]
        # Values from experimental data, see:
        # https://www.frontiersin.org/articles/10.3389/fphys.2019.00580/full
        th_endo = percent_endo
        th_mid = percent_endo + percent_mid
        endo_nodes = (np.nonzero(np.logical_and(values >= 0, values < th_endo)))[0]
        mid_nodes = (np.nonzero(np.logical_and(values >= th_endo, values < th_mid)))[0]
        epi_nodes = (np.nonzero(np.logical_and(values >= th_mid, values <= 1)))[0]
        endo_nodeset_id = self.get_unique_nodeset_id()
        node_set_kw = create_node_set_keyword(
            node_ids=endo_nodes + 1,
            node_set_id=endo_nodeset_id,
            title="Layer-Endo",
        )
        self.kw_database.node_sets.append(node_set_kw)
        mid_nodeset_id = self.get_unique_nodeset_id()
        node_set_kw = create_node_set_keyword(
            node_ids=mid_nodes + 1,
            node_set_id=mid_nodeset_id,
            title="Layer-Mid",
        )
        self.kw_database.node_sets.append(node_set_kw)
        epi_nodeset_id = self.get_unique_nodeset_id()
        node_set_kw = create_node_set_keyword(
            node_ids=epi_nodes + 1,
            node_set_id=epi_nodeset_id,
            title="Layer-Epi",
        )
        self.kw_database.node_sets.append(node_set_kw)
        return endo_nodeset_id, mid_nodeset_id, epi_nodeset_id

    def _add_cell_model_keyword(self, matid: int, cellmodel: CellModel):
        """Add cell model keyword to database."""
        if isinstance(cellmodel, CellModel.Tentusscher):
            self._add_Tentusscher_keyword(matid=matid, params=cellmodel.to_dictionary())
        else:
            raise NotImplementedError

    def _add_Tentusscher_keyword(self, matid: int, params: dict):  # noqa N802
        cell_kw = keywords.EmEpCellmodelTentusscher(**{**params})
        cell_kw.mid = matid
        # Note: bug in EmEpCellmodelTentusscher
        # the following 2 parameters can not be assigned by above method
        cell_kw.gas_constant = 8314.472
        cell_kw.faraday_constant = 96485.3415

        self.kw_database.cell_models.append(cell_kw)

    def _update_ep_settings(self):
        """Add the settings for the electrophysiology solver."""
        save_part_ids = []
        for part in self.model.parts:
            save_part_ids.append(part.pid)
        for beams_pid in self.model.conduction_system._line_id_to_pid.values():
            save_part_ids.append(beams_pid)
        partset_id = self.get_unique_partset_id()
        kw = keywords.SetPartList(sid=partset_id)
        # kw.parts._data = save_part_ids
        # NOTE: when len(save_part_ids) = 8/16, dynalib bugs
        str = "\n"
        for i, id in enumerate(save_part_ids):
            str += "{0:10d}".format(id)
            if (i + 1) % 8 == 0:
                str += "\n"
        kw = kw.write() + str

        self.kw_database.ep_settings.append(kw)
        solvertype = self.settings.electrophysiology.analysis.solvertype
        if solvertype == "Monodomain":
            emsol = 11
            self.kw_database.ep_settings.append(custom_keywords.EmControlEp(numsplit=1))
        elif solvertype == "Eikonal":
            emsol = 14
            self.kw_database.ep_settings.append(custom_keywords.EmControlEp(numsplit=1, ionsolvr=0))
            t_end = 500
            dt = 0.1
            # specify simulation time and time step even in the case of a pure
            # Eikonal model (otherwise LS-DYNA crashes)
            self.kw_database.ep_settings.append("$     Tend        dt")
            self.kw_database.ep_settings.append(f"{t_end:>10f}{dt:>10f}")
        elif solvertype == "ReactionEikonal":
            emsol = 15
            self.kw_database.ep_settings.append(custom_keywords.EmControlEp(numsplit=1, ionsolvr=2))
            t_end = 500
            dt = 0.1
            # specify simulation time and time step in case of a spline ionsolver type
            self.kw_database.ep_settings.append("$     Tend        dt")
            self.kw_database.ep_settings.append(f"{t_end:>10f}{dt:>10f}")

        macrodt = self.settings.electrophysiology.analysis.dtmax.m
        if macrodt > self.settings.mechanics.analysis.dtmax.m:
            LOGGER.info(
                "EP Timestep > Mechanics Timestep. Setting EP Timestep to Mechanics Timestep"
            )
            macrodt = self.settings.mechanics.analysis.dtmax.m

        self.kw_database.ep_settings.append(
            keywords.EmControl(
                emsol=emsol,
                numls=4,
                macrodt=macrodt,
                dimtype=None,
                nperio=None,
                ncylbem=None,
            )
        )
        self.kw_database.ep_settings.append(keywords.EmControlTimestep(dtcons=macrodt))

        self.kw_database.ep_settings.append(
            custom_keywords.EmEpIsoch(idisoch=1, idepol=1, dplthr=-20, irepol=1, rplthr=-40)
        )

        self.kw_database.ep_settings.append(
            keywords.EmSolverFem(reltol=1e-6, maxite=int(1e4), precon=2)
        )

        self.kw_database.ep_settings.append(keywords.EmOutput(mats=1, matf=1, sols=1, solf=1))

    def _update_stimulation(self):
        # define stimulation settings
        stimsettings = self.settings.electrophysiology.stimulation
        if not stimsettings:
            stim_nodes = self.get_default_stimulus_nodes()
            stimulation = Stimulation(node_ids=stim_nodes)

            stimsettings = {"stimdefaults": stimulation}

        for stimname in stimsettings.keys():
            stim_nodes = stimsettings[stimname].node_ids
            if stimsettings[stimname].node_ids is None:
                stim_nodes = self.get_default_stimulus_nodes()
            stim = Stimulation(
                node_ids=stim_nodes,
                t_start=stimsettings[stimname].t_start,
                period=stimsettings[stimname].period,
                duration=stimsettings[stimname].duration,
                amplitude=stimsettings[stimname].amplitude,
            )
            node_set_kw, stim_kw = self._add_stimulation_keyword(stim)
            self.kw_database.ep_settings.append(node_set_kw)
            self.kw_database.ep_settings.append(stim_kw)

    def _add_stimulation_keyword(self, stim: Stimulation):
        # create node-sets for stim nodes
        nsid = self.get_unique_nodeset_id()
        node_set_kw = create_node_set_keyword(
            node_ids=np.array(stim.node_ids) + 1,
            node_set_id=nsid,
            title="Stim nodes",
        )

        solvertype = self.settings.electrophysiology.analysis.solvertype
        if solvertype == "Monodomain":
            stim_kw = custom_keywords.EmEpTentusscherStimulus(
                stimid=nsid,
                settype=2,
                setid=nsid,
                stimstrt=stim.t_start.m,
                stimt=stim.period.m,
                stimdur=stim.duration.m,
                stimamp=stim.amplitude.m,
            )

        else:
            # TODO: : add eikonal in custom keywords
            # EM_EP_EIKONAL

            eikonal_stim_content = "*EM_EP_EIKONAL\n"
            eikonal_stim_content += "$    eikId  eikPaSet eikStimNS eikStimDF\n"
            # TODO: get the right part set id
            # setpart_kwds = self.kw_database.ep_settings.get_kwds_by_type()
            # id of the eikonal solver (different eikonal solves
            # can be performed in different parts of the model)
            eikonal_id = 1
            psid = 1
            eikonal_stim_content += f"{eikonal_id:>10d}{psid:>10d}{nsid:>10d}\n"
            if solvertype == "ReactionEikonal":
                eikonal_stim_content += "$ footType     footT     footA  footTauf   footVth\n"
                foot_type = 1
                foot_t = stim.duration.m
                foot_a = stim.amplitude.m
                foot_tauf = 1
                eikonal_stim_content += (
                    f"{foot_type:>10d}{foot_t:>10f}{foot_a:>10f}{foot_tauf:>10f}"
                )
                eikonal_stim_content += "\n$solvetype\n"
                eikonal_stim_content += f"{1:>10d}"  # activate time stepping method by default
            stim_kw = eikonal_stim_content

        return (node_set_kw, stim_kw)

    def get_default_stimulus_nodes(self) -> list[int]:
        """Get default stiumulus nodes.

        1/2 apex point(s) for Left/Bi-ventricle model.

        Sinoatrial node for Fourchamber/Full heart model

        Returns
        -------
        list[int]
            0-based node IDs to sitmulate
        """
        if isinstance(self.model, LeftVentricle):
            stim_nodes = [self.model.left_ventricle.apex_points[0].node_id]

        elif isinstance(self.model, BiVentricle):
            node_apex_left = self.model.left_ventricle.apex_points[0].node_id
            node_apex_right = self.model.right_ventricle.apex_points[0].node_id
            stim_nodes = [node_apex_left, node_apex_right]

        elif isinstance(self.model, (FourChamber, FullHeart)):
            node_apex_left = self.model.left_ventricle.apex_points[0].node_id
            node_apex_right = self.model.right_ventricle.apex_points[0].node_id
            stim_nodes = [node_apex_left, node_apex_right]

            if self.model.right_atrium.get_point("SA_node") is not None:
                # Active SA node (belong to both solid and beam)
                stim_nodes = list(
                    self.model.mesh.find_closest_point(
                        self.model.right_atrium.get_point("SA_node").xyz, n=5
                    )
                )

                #  add more nodes to initiate wave propagation
                if _ConductionType.SAN_AVN.value in list(
                    self.model.conduction_system._line_id_to_name.values()
                ):
                    pointid = self.model.conduction_system.get_lines_by_name(
                        _ConductionType.SAN_AVN.value
                    )["_written-id"][1]
                    stim_nodes.append(pointid)
                if _ConductionType.BACHMANN_BUNDLE.value in list(
                    self.model.conduction_system._line_id_to_name.values()
                ):
                    pointid = self.model.conduction_system.get_lines_by_name(
                        _ConductionType.BACHMANN_BUNDLE.value
                    )["_written-id"][0]
                    stim_nodes.append(pointid)
                    pointid = self.model.conduction_system.get_lines_by_name(
                        _ConductionType.BACHMANN_BUNDLE.value
                    )["_written-id"][1]
                    stim_nodes.append(pointid)

        # stimulate entire elements for Eikonal
        if self.settings.electrophysiology.analysis.solvertype in [
            "Eikonal",
            "ReactionEikonal",
        ]:
            stim_cells = np.where(np.isin(self.model.mesh.tetrahedrons, stim_nodes))[0]
            stim_nodes = np.unique(self.model.mesh.tetrahedrons[stim_cells].ravel())

        return stim_nodes

    def _update_blood_settings(self):
        """Update blood settings."""
        if self.model._add_blood_pool:
            dirichlet_bc_nid = self.get_unique_nodeset_id()
            apex = self.model.left_ventricle.apex_points[0].node_id
            node_set_kw = create_node_set_keyword(
                node_ids=apex + 1,
                node_set_id=dirichlet_bc_nid,
                title="Dirichlet extracellular potential BC",
            )
            self.kw_database.node_sets.append(node_set_kw)
            self.kw_database.ep_settings.append(
                custom_keywords.EmBoundaryPrescribed(
                    bpid=1,
                    bptype=1,
                    settype=2,
                    setid=dirichlet_bc_nid,
                    val=0,
                    sys=0,
                )
            )
            for deckname, deck in vars(self.kw_database).items():
                # lambda_ is the equal anisotropy ratio in the monodomain model.
                # In dyna: lambda_= sigma_i/sigma_e and sigma_i=(1.+lambda)*sigmaElement.
                # when lambda_ is not empty, it activates the computation of extracellular
                # potentials: div((sigma_i+sigma_e) . grad(phi_e)) = div(sigma_i . grad(v))
                # or div(((1.+lambda)*sigmaElement) . grad(phi_e)) = div(sigmaElement . grad(v))
                for kw in deck.keywords:
                    # activate extracellular potential solve
                    if "EM_MAT" in kw.get_title():
                        kw.lambda_ = self.settings.electrophysiology.material.myocardium["lambda"].m

    def _update_ECG_coordinates(self):  # noqa N802
        """Add ECG computation content."""
        # TODO: replace strings by custom dyna keyword
        # TODO: handle dynamic numbering of point set ids "psid'
        psid = 1
        pstype = 0

        # EM_POINT_SET
        em_point_set_content = "*EM_POINT_SET\n"
        em_point_set_content += "$#    psid    pstype        vx        vy        vz\n"
        em_point_set_content += f"{psid:>10d}{pstype:>10d}\n"
        em_point_set_content += "$#     pid         x         y         z       pos"

        self.kw_database.ep_settings.append(em_point_set_content)

        for index, point in enumerate(self.model.electrodes):
            x, y, z = point.xyz
            position_str = (
                f"{index:>10d} {str(f'{x:9.6f}')[:9]} {str(f'{y:9.6f}')[:9]} {str(f'{z:9.6f}')[:9]}"  # noqa
            )

            self.kw_database.ep_settings.append(position_str)

        # EM_EP_EKG
        em_ep_ekg_content = "*EM_EP_EKG\n"
        em_ep_ekg_content += "$#   ekgid      psid\n"
        em_ep_ekg_content += f"{1:>10d}{psid:>10d}\n"

        self.kw_database.ep_settings.append(em_ep_ekg_content)

    def _update_solution_controls(
        self,
    ):
        """Add solution controls and other solver settings as keywords."""
        self.kw_database.main.append(
            keywords.ControlTermination(
                endtim=self.settings.electrophysiology.analysis.end_time.m,
                dtmin=self.settings.electrophysiology.analysis.dtmin.m,
            )
        )
        self.kw_database.main.append(
            keywords.ControlTimeStep(
                dtinit=self.settings.electrophysiology.analysis.dtmax.m,
                dt2ms=self.settings.electrophysiology.analysis.dtmax.m,
            )
        )
        return

    def _update_main_db(self):
        pass

    def _update_use_Purkinje(self, associate_to_segment: bool = True):  # noqa N802
        """Update keywords for Purkinje usage."""
        if not isinstance(self.model, (FullHeart, FourChamber, BiVentricle, LeftVentricle)):
            LOGGER.error("Model type not recognized.")
            return

        sid = self.get_unique_section_id()
        self.kw_database.beam_networks.append(keywords.SectionBeam(secid=sid, elform=3, a=645))

        if type(self) is ElectroMechanicsDynaWriter:
            # id offset due to spring-type elements in mechanical
            beam_elem_id_offset = self.id_offset["element"]["discrete"]
        else:
            beam_elem_id_offset = 0  # no beam elements introduced before

        # write beam nodes
        # Note: the last beam_network saves all beam nodes
        new_nodes = self.model.conduction_system.points[
            (np.where(np.logical_not(self.model.conduction_system["_is-connected"]))[0])
        ]
        ids = (
            np.linspace(
                len(self.model.mesh.points),
                len(self.model.mesh.points) + len(new_nodes) - 1,
                len(new_nodes),
                dtype=int,
            )
            + 1  # dyna start by 1
        )
        nodes_table = np.hstack((ids.reshape(-1, 1), new_nodes))
        kw = add_nodes_to_kw(nodes_table, keywords.Node())
        self.kw_database.beam_networks.append(kw)
        material_settings = self.settings.electrophysiology.material
        solvertype = self.settings.electrophysiology.analysis.solvertype
        default_epmat = EPMaterial.ActiveBeam()
        if solvertype == "Monodomain":
            sig1 = material_settings.beam["sigma"].m
        else:
            sig1 = material_settings.beam["velocity"].m
        default_epmat.sigma_fiber = sig1
        default_epmat.beta = material_settings.beam["beta"].m
        default_epmat.cm = material_settings.beam["cm"].m
        default_epmat.pmjres = material_settings.beam["pmjres"].m
        beam_point_offset_id = 0
        self.model.conduction_system.point_data["_written-id"] = (
            np.zeros(self.model.conduction_system.number_of_points, dtype=int) - 1
        )
        # new for loop
        for netid in self.model.conduction_system._line_id_to_name:
            if isinstance(
                self.model.conduction_system.ep_material[netid], EPMaterial.DummyMaterial
            ):
                epmat = default_epmat
            else:
                epmat = self.model.conduction_system.ep_material[netid]
            pid = self.get_unique_part_id()
            self.model.conduction_system._line_id_to_pid[netid] = pid
            name = self.model.conduction_system._line_id_to_name[netid]
            if name == _ConductionType.LEFT_PURKINJE.value:
                _node_set_id = self.model.left_ventricle.endocardium._seg_set_id
            elif name == _ConductionType.RIGHT_PURKINJE.value:
                _node_set_id = self.model.right_ventricle.endocardium._seg_set_id
            elif name == _ConductionType.SAN_AVN.value:
                _node_set_id = self.model.right_atrium.endocardium._seg_set_id
            elif name == _ConductionType.LEFT_BUNDLE_BRANCH.value:
                _node_set_id = self.model.left_ventricle.cavity.surface._seg_set_id
            elif name == _ConductionType.RIGHT_BUNDLE_BRANCH.value:
                _node_set_id = self.model.right_ventricle.cavity.surface._seg_set_id
            elif name == _ConductionType.HIS.value:
                # His bundle are inside of 3d mesh
                # need to create the segment on which beam elements rely
                surface = self._add_segment_from_surface(name="his_bundle_segment")
                _node_set_id = surface._seg_set_id
            elif name == _ConductionType.BACHMANN_BUNDLE.value:
                # His bundle are inside of 3d mesh
                # need to create the segment on which beam elements rely
                surface = self._add_segment_from_surface(name="Bachman segment")
                _node_set_id = surface._seg_set_id
            else:
                LOGGER.error(f"Unknown network name for {name}.")
                exit()

            # overwrite nsid if beam should not follow the motion of segment
            if not associate_to_segment:
                _node_set_id = -1

            # write
            self.kw_database.beam_networks.append(f"$$ {name} $$")
            origin_coordinates = self.model.conduction_system.get_lines(netid).points[0]
            self.kw_database.beam_networks.append(
                custom_keywords.EmEpPurkinjeNetwork2(
                    purkid=pid,
                    buildnet=0,
                    ssid=_node_set_id,
                    mid=pid,
                    pointstx=origin_coordinates[0],
                    pointsty=origin_coordinates[1],
                    pointstz=origin_coordinates[2],
                    edgelen=self.settings.purkinje.edgelen.m,
                    ngen=self.settings.purkinje.ngen.m,
                    nbrinit=self.settings.purkinje.nbrinit.m,
                    nsplit=self.settings.purkinje.nsplit.m,
                    pmjtype=self.settings.purkinje.pmjtype.m,
                    pmjradius=self.settings.purkinje.pmjradius.m,
                    pmjrestype=self.settings.electrophysiology.material.beam["pmjrestype"].m,
                    pmjres=epmat.pmjres,
                )
            )

            part_df = pd.DataFrame(
                {
                    "heading": [name],
                    "pid": [pid],
                    "secid": [sid],
                    "mid": [pid],
                }
            )
            part_kw = keywords.Part()
            part_kw.parts = part_df
            self.kw_database.beam_networks.append(part_kw)
            self.kw_database.beam_networks.append(keywords.MatNull(mid=pid, ro=1e-11))

            kw = self._get_ep_material_kw(pid, epmat)
            self.kw_database.beam_networks.append(kw)

            # cell model
            self._add_cell_model_keyword(matid=pid, cellmodel=epmat.cell_model)

            # Build connectivity
            # get edges in a 2 column format
            edges = self.model.conduction_system.get_lines(netid).lines.reshape(
                (int(len(self.model.conduction_system.get_lines(netid).lines) / 3), 3)
            )[:, 1:]

            # get info on points to be connected to solid and their coordinates
            connected_point_ids = np.where(
                self.model.conduction_system.get_lines(netid)["_is-connected"]
            )[0]
            mask_nonconnected = np.logical_not(
                self.model.conduction_system.get_lines(netid)["_is-connected"]
            )
            connected_points = self.model.conduction_system.get_lines(netid).points[
                connected_point_ids
            ]

            # got ids in solid mesh of connected points
            kdtree = spatial.cKDTree(self.model.mesh.points)
            _, solid_connected_point_ids = kdtree.query(connected_points)

            # compute writer point ids depending on previously written and connections to solid
            point_ids_to_write = np.zeros(
                self.model.conduction_system.get_lines(netid).number_of_points, dtype=int
            )
            point_ids_to_write[connected_point_ids] = solid_connected_point_ids
            mask_already_written = self.model.conduction_system.get_lines(netid)["_written-id"] >= 0
            point_ids_to_write[mask_already_written] = self.model.conduction_system.get_lines(
                netid
            )["_written-id"][mask_already_written]
            mask_notconnected_notwritten = np.logical_and(
                mask_nonconnected, np.logical_not(mask_already_written)
            )
            point_ids_to_write[mask_notconnected_notwritten] = (
                np.linspace(
                    0,
                    np.sum(mask_notconnected_notwritten) - 1,
                    np.sum(mask_notconnected_notwritten),
                    dtype=int,
                )
                + self.model.mesh.number_of_points
                + beam_point_offset_id
            )

            # replace point id values in edges
            edges = np.vectorize(lambda idvalue: point_ids_to_write[idvalue])(edges)

            # write mesh
            beams_kw = keywords.ElementBeam()
            beams_kw = add_beams_to_kw(
                beams=edges + 1,
                beam_kw=beams_kw,
                pid=pid,
                offset=beam_elem_id_offset,
            )
            # offset beam element id
            beam_elem_id_offset += len(edges)
            # offset beam point id
            beam_point_offset_id += (
                self.model.conduction_system.get_lines(netid).number_of_points
                - len(connected_point_ids)
                - np.sum(mask_already_written)
            )
            # populate the already written ids variable for other beam networks
            point_ids_in_conductionsystem = self.model.conduction_system.get_lines(netid)[
                "_global-point-ids"
            ]
            self.model.conduction_system["_written-id"][point_ids_in_conductionsystem] = (
                point_ids_to_write
            )
            self.kw_database.beam_networks.append(beams_kw)

        self.id_offset["element"]["discrete"] = beam_elem_id_offset

        return

    def _add_segment_from_surface(self, name: str):
        surface = self.model.mesh.get_surface_by_name(name)

        surface._seg_set_id = self.get_unique_segmentset_id()
        surface._node_set_id = self.get_unique_nodeset_id()

        kw = create_segment_set_keyword(
            segments=surface.triangles_global + 1,
            segid=surface._seg_set_id,
            title=surface.name,
        )
        # append this kw to the segment set database
        self.kw_database.segment_sets.append(kw)

        return surface

    def _update_export_controls(self):
        """Add solution controls to the main simulation."""
        self.kw_database.main.append(
            keywords.DatabaseBinaryD3Plot(dt=self.settings.electrophysiology.analysis.dt_d3plot.m)
        )

        return

    def _get_ep_material_kw(self, ep_mid: int, ep_material: EPMaterial):
        if type(ep_material) is EPMaterial.Insulator:
            # insulator mtype
            mtype = 1
            kw = custom_keywords.EmMat001(
                mid=ep_mid,
                mtype=mtype,
                sigma=ep_material.sigma_fiber,
                beta=ep_material.beta,
                cm=ep_material.cm,
            )

        # active myocardium
        elif type(ep_material) is EPMaterial.Active:
            mtype = 2
            # "isotropic" case
            if ep_material.sigma_sheet is None:
                # lSDYNA bug prevents from using isotropic mat (EMMAT001) for active isotropic case
                # Bypass: using EMMAT003 with same sigma value in all directions
                ep_material.sigma_sheet = ep_material.sigma_fiber
                ep_material.sigma_sheet_normal = ep_material.sigma_fiber
            kw = custom_keywords.EmMat003(
                mid=ep_mid,
                mtype=mtype,
                sigma11=ep_material.sigma_fiber,
                sigma22=ep_material.sigma_sheet,
                sigma33=ep_material.sigma_sheet_normal,
                beta=ep_material.beta,
                cm=ep_material.cm,
                aopt=2.0,
                a1=0,
                a2=0,
                a3=1,
                d1=0,
                d2=-1,
                d3=0,
            )

        elif type(ep_material) is EPMaterial.ActiveBeam:
            mtype = 2
            kw = custom_keywords.EmMat001(
                mid=ep_mid,
                mtype=mtype,
                sigma=ep_material.sigma_fiber,
                beta=ep_material.beta,
                cm=ep_material.cm,
            )
        elif type(ep_material) is EPMaterial.Passive:
            mtype = 4
            # isotropic
            if ep_material.sigma_sheet is None:
                kw = custom_keywords.EmMat001(
                    mid=ep_mid,
                    mtype=mtype,
                    sigma=ep_material.sigma_fiber,
                    beta=ep_material.beta,
                    cm=ep_material.cm,
                )
            # Anisotropic
            else:
                kw = custom_keywords.EmMat003(
                    mid=ep_mid,
                    mtype=mtype,
                    sigma11=ep_material.sigma_fiber,
                    sigma22=ep_material.sigma_sheet,
                    sigma33=ep_material.sigma_sheet_normal,
                    beta=ep_material.beta,
                    cm=ep_material.cm,
                    aopt=2.0,
                    a1=0,
                    a2=0,
                    a3=1,
                    d1=0,
                    d2=-1,
                    d3=0,
                )
        return kw


class ElectrophysiologyBeamsDynaWriter(ElectrophysiologyDynaWriter):
    """Class for preparing the input for an Electrophysiology LS-DYNA simulation with beams only."""

    def __init__(self, model: HeartModel, settings: SimulationSettings = None) -> None:
        super().__init__(model=model, settings=settings)
        self.kw_database = ElectrophysiologyDecks()
        """Collection of keywords relevant for Electrophysiology."""

    def update(self):
        """Update keyword database for Electrophysiology."""
        # self._isolate_atria_and_ventricles()

        ##
        self._update_main_db()
        self._update_solution_controls()
        self._update_export_controls()

        self._update_node_db()

        if self.model.conduction_system.number_of_cells != 0:
            # with smcoupl=1, coupling is disabled
            self.kw_database.ep_settings.append(keywords.EmControlCoupling(thcoupl=1, smcoupl=1))
            self._update_use_Purkinje(associate_to_segment=False)

        # update ep settings
        self._update_ep_settings()
        self._update_stimulation()

        self._get_list_of_includes()
        self._add_includes()

        return


class ElectroMechanicsDynaWriter(MechanicsDynaWriter, ElectrophysiologyDynaWriter):
    """Class for preparing the input for LS-DYNA electromechanical simulation."""

    def __init__(
        self,
        model: HeartModel,
        settings: SimulationSettings = None,
    ) -> None:
        if isinstance(model, FourChamber):
            model._create_atrioventricular_isolation()

        BaseDynaWriter.__init__(self, model=model, settings=settings)

        self.kw_database = ElectroMechanicsDecks()
        """Collection of keyword decks relevant for mechanics."""

        self.system_model_name = self.settings.mechanics.system.name
        """Name of system model to use, from MechanicWriter."""

        self.set_flow_area = True
        """from MechanicWriter."""

    def update(self, dynain_name: str = None, robin_bcs=None):
        """Update the keyword database.

        Parameters
        ----------
        dynain_name : str, optional
            dynain file from stress free configuration computation, by default None
        robin_bcs : list[Callable], optional
            A list of lambda functions to apply Robin-type BCs, by default None

        Notes
        -----
        Do not need to write mesh files if dynain file is given.
        """
        if isinstance(self.model, FourChamber):
            self.model.left_atrium.fiber = True
            self.model.left_atrium.active = True
            self.model.right_atrium.fiber = True
            self.model.right_atrium.active = True

        MechanicsDynaWriter.update(self, dynain_name=dynain_name, robin_bcs=robin_bcs)

        if self.model.conduction_system.number_of_cells != 0:
            # Coupling enabled, EP beam nodes follow the motion of surfaces
            self.kw_database.ep_settings.append(keywords.EmControlCoupling(thcoupl=1, smcoupl=0))
            self._update_use_Purkinje()
            self.kw_database.main.append(keywords.Include(filename="beam_networks.k"))

        self._update_parts_cellmodels()
        self.kw_database.main.append(keywords.Include(filename="cell_models.k"))

        self._update_ep_settings()
        self._update_stimulation()

        # coupling parameters
        coupling_str = (
            "*EM_CONTROL_COUPLING\n$    THCPL     SMCPL    THLCID    SMLCID\n         1         0\n"
        )
        self.kw_database.ep_settings.append("$ EM-MECA coupling control")
        self.kw_database.ep_settings.append(coupling_str)
        self.kw_database.main.append(keywords.Include(filename="ep_settings.k"))

        return

    def _update_material_db(self, add_active: bool = True):
        """Update the database of material keywords."""
        MechanicsDynaWriter._update_material_db(self, add_active=add_active, em_couple=True)
        ElectrophysiologyDynaWriter._update_ep_material_db(self)
        return


class LaplaceWriter(BaseDynaWriter):
    """Writer to set Laplace dirichlet problem."""

    # constant node set ID for atrial valves/caps
    _CAP_NODESET_MAP = {
        CapType.RIGHT_INFERIOR_PULMONARY_VEIN: 1,
        CapType.LEFT_ATRIUM_APPENDAGE: 2,
        CapType.RIGHT_SUPERIOR_PULMONARY_VEIN: 3,
        CapType.MITRAL_VALVE_ATRIUM: 4,
        CapType.LEFT_INFERIOR_PULMONARY_VEIN: 5,
        CapType.LEFT_SUPERIOR_PULMONARY_VEIN: 6,
        CapType.TRICUSPID_VALVE_ATRIUM: 7,
        CapType.SUPERIOR_VENA_CAVA: 8,
        CapType.INFERIOR_VENA_CAVA: 9,
    }
    _LANDMARK_RADIUS = 1.5  # mm
    _UVC_APEX_RADIUS = 10.0  # mm

    def __init__(
        self, model: HeartModel, type: Literal["uvc", "la_fiber", "ra_fiber", "D-RBM"], **kwargs
    ):
        """Write thermal input to set up a Laplace dirichlet problem.

        Parameters
        ----------
        model : HeartModel
            Heart model
        type : Literal[&quot;uvc&quot;, &quot;la_fiber&quot;, &quot;ra_fiber&quot;, &quot;D
            simulation type
        """
        super().__init__(model=model)
        self.type = type
        """problem type."""
        self.landmarks = kwargs
        """landmarks can be `laa`, `raa`, `top`."""
        self.target: pv.UnstructuredGrid = None
        """target mesh related to the problem."""

        # remove unnecessary parts
        if self.type == "uvc" or self.type == "D-RBM":
            parts_to_keep = ["Left ventricle", "Right ventricle", "Septum"]
            self._keep_parts(parts_to_keep)
        elif self.type == "la_fiber":
            parts_to_keep = ["Left atrium"]
        elif self.type == "ra_fiber":
            parts_to_keep = ["Right atrium"]

        # remove unnecessary mesh and create target attribute
        if self.type == "uvc" or self.type == "D-RBM":
            elems_to_keep = []
            if isinstance(self.model, LeftVentricle):
                elems_to_keep.extend(model.left_ventricle.element_ids)
            else:
                elems_to_keep.extend(model.left_ventricle.element_ids)
                elems_to_keep.extend(model.right_ventricle.element_ids)
                elems_to_keep.extend(model.septum.element_ids)

            # model.mesh.clear_data()
            model.mesh["cell_ids"] = np.arange(0, model.mesh.n_cells, dtype=int)
            model.mesh["point_ids"] = np.arange(0, model.mesh.n_points, dtype=int)

            self.target = model.mesh.extract_cells(elems_to_keep)

        elif self.type == "la_fiber" or self.type == "ra_fiber":
            self._keep_parts(parts_to_keep)
            # model.mesh.clear_data()
            model.mesh["cell_ids"] = np.arange(0, model.mesh.n_cells, dtype=int)
            model.mesh["point_ids"] = np.arange(0, model.mesh.n_points, dtype=int)

            self.target = model.mesh.extract_cells(model.parts[0].element_ids)

    def _update_ra_top_nodeset(self, atrium: pv.UnstructuredGrid):
        """
        Define right atrium top nodeset with node set id 10.

        Parameters
        ----------
        atrium : pv.UnstructuredGrid
            right atrium pyvista object
        """
        if "top" in self.landmarks.keys():
            top_ids = self._find_top_nodeset_by_geodesic(atrium)
        else:
            top_ids = self._find_top_nodeset_by_cut(atrium)

        # assign top nodeset
        kw = create_node_set_keyword(top_ids + 1, node_set_id=10, title="top")
        self.kw_database.node_sets.append(kw)

    def _find_top_nodeset_by_cut(self, atrium: pv.UnstructuredGrid):
        """
        Define right atrium top nodeset.

        Cut through the center of TV, IVC and SVC, expecting to result in
        3 unconnected regions and the farthest is top.
        This method may fail with varying geometries, then the user
        needs to define the top landmarks.
        """
        cut_center, cut_normal = self._define_ra_cut()

        atrium["cell_ids_tmp"] = np.arange(0, atrium.n_cells, dtype=int)
        atrium["point_ids_tmp"] = np.arange(0, atrium.n_points, dtype=int)
        slice = atrium.slice(origin=cut_center, normal=cut_normal)
        crinkled = atrium.extract_cells(np.unique(slice["cell_ids_tmp"]))

        # After cut, select the top region
        x = crinkled.connectivity()
        if np.max(x.point_data["RegionId"]) != 2:
            # Should only have 3 parts
            LOGGER.error("Cannot find top nodeset...")
            raise ValueError("Please define top start/end points and re-run.")

        # get tricuspid-valve name
        tv_name = CapType.TRICUSPID_VALVE_ATRIUM.value

        # compare closest point with TV nodes, top region should be far with TV node set
        tv_tree = spatial.cKDTree(atrium.points[atrium.point_data[tv_name] == 1])
        min_dst = -1.0
        for i in range(3):
            current_min_dst = np.min(tv_tree.query(x.points[x.point_data["RegionId"] == i])[0])
            if current_min_dst > min_dst:
                min_dst = current_min_dst
                top_region_id = i

        # This region is the top
        mask = x.point_data["RegionId"] == top_region_id

        top_ids = x["point_ids_tmp"][mask]

        atrium.cell_data.remove("cell_ids_tmp")
        atrium.point_data.remove("point_ids_tmp")
        return top_ids

    def _find_top_nodeset_by_geodesic(self, atrium: pv.UnstructuredGrid):
        """Define top nodeset by connecting landmark points with a geodesic path."""
        top_ids = []
        surface: pv.PolyData = atrium.extract_surface()
        for i in range(len(self.landmarks["top"]) - 1):
            p1 = self.landmarks["top"][i]
            p2 = self.landmarks["top"][i + 1]

            path = surface.geodesic(surface.find_closest_point(p1), surface.find_closest_point(p2))
            for point in path.points:
                top_ids.append(atrium.find_closest_point(point))

        return np.unique(np.array(top_ids))

    def _define_ra_cut(self):
        """Define a cut-plane using the three caps of right atrium."""
        for cap in self.model.parts[0].caps:
            if cap.type == CapType.TRICUSPID_VALVE_ATRIUM:
                tv_center = cap.centroid
            elif cap.type == CapType.SUPERIOR_VENA_CAVA:
                svc_center = cap.centroid
            elif cap.type == CapType.INFERIOR_VENA_CAVA:
                ivc_center = cap.centroid
        cut_center = np.vstack((tv_center, svc_center, ivc_center)).mean(axis=0)
        cut_normal = np.cross(svc_center - tv_center, ivc_center - tv_center)

        return cut_center, cut_normal

    def _update_ra_tricuspid_nodeset(self, atrium):
        """Define nodeset for tricuspid_wall and tricuspid_septum."""
        # get tricuspid-valve name
        tv_name = CapType.TRICUSPID_VALVE_ATRIUM.value

        # The cut_normal is determined so 1st part will be septum and 2nd will be free
        cut_center, cut_normal = self._define_ra_cut()

        # need a copied object to do clip, atrium will be corrupted otherwise
        septum, free_wall = copy.deepcopy(atrium).clip(
            origin=cut_center, normal=cut_normal, crinkle=True, return_clipped=True
        )
        # ids in full mesh
        tv_s_ids = septum["point_ids"][np.where(septum[tv_name] == 1)]

        tv_s_ids_sub = np.where(np.isin(atrium["point_ids"], tv_s_ids))[0]
        atrium["tv_s"] = np.zeros(atrium.n_points)
        atrium["tv_s"][tv_s_ids_sub] = 1

        kw = create_node_set_keyword(tv_s_ids_sub + 1, node_set_id=12, title="tv_septum")
        self.kw_database.node_sets.append(kw)

        tv_w_ids = free_wall["point_ids"][np.where(free_wall[tv_name] == 1)]
        tv_w_ids_sub = np.where(np.isin(atrium["point_ids"], tv_w_ids))[0]
        # remove re constraint nodes
        tv_w_ids_sub = np.setdiff1d(tv_w_ids_sub, tv_s_ids_sub)

        atrium["tv_w"] = np.zeros(atrium.n_points)
        atrium["tv_w"][tv_w_ids_sub] = 1

        kw = create_node_set_keyword(tv_w_ids_sub + 1, node_set_id=13, title="tv_wall")
        self.kw_database.node_sets.append(kw)

    def _update_atrial_caps_nodeset(self, atrium: pv.UnstructuredGrid):
        """Define node sets for the caps."""
        for cap in self.model.parts[0].caps:
            # get node IDs for atrium mesh
            cap._mesh = self.model.mesh.get_surface(cap._mesh.id)
            ids_sub = np.where(np.isin(atrium["point_ids"], cap.global_node_ids_edge))[0]
            # create node set
            set_id = self._CAP_NODESET_MAP[cap.type]

            if set_id:  # Can be None for LEFT_ATRIUM_APPENDAGE
                kw = create_node_set_keyword(ids_sub + 1, node_set_id=set_id, title=cap.name)
                self.kw_database.node_sets.append(kw)

                # Add info to pyvista object, necessary for right atrial fibers.
                atrium[cap.type.value] = np.zeros(atrium.n_points, dtype=int)
                atrium[cap.type.value][ids_sub] = 1

        return

    def _update_la_bc(self):
        atrium = self.target

        def get_laa_nodes(atrium, laa: np.ndarray):
            tree = spatial.cKDTree(atrium.points)
            ids = np.array(tree.query_ball_point(laa, self._LANDMARK_RADIUS))
            return ids

        # laa
        if "laa" in self.landmarks.keys():
            # else there should be a LEFT_ATRIUM_APPENDAGE as in Strocchi's data
            laa_ids = get_laa_nodes(atrium, self.landmarks["laa"])

            kw = create_node_set_keyword(
                laa_ids + 1,
                node_set_id=self._CAP_NODESET_MAP[CapType.LEFT_ATRIUM_APPENDAGE],
                title="left atrium appendage",
            )
            self.kw_database.node_sets.append(kw)

        # caps
        self._update_atrial_caps_nodeset(atrium)

        # endo/epi
        endo_nodes = self.model.left_atrium.endocardium.global_node_ids_triangles
        epi_nodes = self.model.left_atrium.epicardium.global_node_ids_triangles
        epi_nodes = np.setdiff1d(epi_nodes, endo_nodes)

        self._add_nodeset(endo_nodes, "endocardium", nodeset_id=100)
        self._add_nodeset(epi_nodes, "epicardium", nodeset_id=200)

        cases = [
            (1, "trans", [100, 200], [0, 1]),
            (2, "ab", [1, 3, 4, 5, 6, 2], [2.0, 2.0, 1.0, 0.0, 0.0, -1.0]),
            (3, "v", [1, 3, 5, 6], [1.0, 1.0, 0.0, 0.0]),
            (4, "r", [4, 1, 2, 3, 5, 6], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        ]
        for case_id, job_name, set_ids, bc_values in cases:
            self.add_case(case_id, job_name, set_ids, bc_values)

    def _update_ra_bc(self):
        atrium = self.target
        # caps
        self._update_atrial_caps_nodeset(atrium)

        # endo/epi
        endo_nodes = self.model.right_atrium.endocardium.global_node_ids_triangles
        epi_nodes = self.model.right_atrium.epicardium.global_node_ids_triangles
        epi_nodes = np.setdiff1d(epi_nodes, endo_nodes)

        self._add_nodeset(endo_nodes, "endocardium", nodeset_id=100)
        self._add_nodeset(epi_nodes, "epicardium", nodeset_id=200)

        # Find appendage apex
        tree = spatial.cKDTree(atrium.points)
        raa_ids = np.array(tree.query_ball_point(self.landmarks["raa"], self._LANDMARK_RADIUS))
        if len(raa_ids) == 0:
            LOGGER.error("No node is identified as right atrium appendage apex.")
            exit()

        kw = create_node_set_keyword(raa_ids + 1, node_set_id=11, title="raa")
        self.kw_database.node_sets.append(kw)
        atrium["raa"] = np.zeros(atrium.n_points)
        atrium["raa"][raa_ids] = 1

        # top nodeset
        self._update_ra_top_nodeset(atrium)
        # tricuspid wall/free nodeset
        self._update_ra_tricuspid_nodeset(atrium)

        cases = [
            (1, "trans", [100, 200], [0, 1]),
            (2, "ab", [9, 7, 8, 11], [2.0, 1.0, 0.0, -1.0]),
            (3, "v", [9, 8, 11], [1.0, 0.0, 0.0]),
            (4, "r", [7, 10], [1.0, 0.0]),
            (5, "w", [12, 13, 10], [1.0, -1.0, 0.0]),
        ]
        for case_id, job_name, set_ids, bc_values in cases:
            self.add_case(case_id, job_name, set_ids, bc_values)

    def update(self):
        """Update keyword database."""
        # nodes
        node_kw = create_node_keyword(self.target.points)
        self.kw_database.nodes.append(node_kw)

        # part and mat
        self._update_parts_materials_db()

        # elems
        kw_elements = create_elemetn_solid_keyword(
            self.target.cells.reshape(-1, 5)[:, 1:] + 1,
            np.arange(1, self.target.n_cells + 1, dtype=int),
            self.model.parts[0].pid,
        )
        self.kw_database.solid_elements.append(kw_elements)

        # main
        self._update_main_db()

        if self.type == "uvc":
            self._update_uvc_bc()
        elif self.type == "la_fiber":
            self._update_la_bc()
        elif self.type == "ra_fiber":
            self._update_ra_bc()
        elif self.type == "D-RBM":
            self._update_drbm_bc()

        self._get_list_of_includes()
        self._add_includes()

    def _update_uvc_bc(self):
        # transmural uvc
        endo_nodes = self.model.left_ventricle.endocardium.global_node_ids_triangles
        epi_nodes = self.model.left_ventricle.epicardium.global_node_ids_triangles

        if not isinstance(self.model, LeftVentricle):
            rv_endo = self.model.right_ventricle.endocardium.global_node_ids_triangles
            septum_endo = self._get_rv_septum_endo_surface().global_node_ids_triangles
            rv_epi = self.model.right_ventricle.epicardium.global_node_ids_triangles

            # septum endocardium is merged into epicardium set, this is
            # consistent with transmural values of LeftVentricle model
            endo_nodes = np.hstack((endo_nodes, rv_endo))
            epi_nodes = np.hstack(
                (
                    epi_nodes,
                    rv_epi,
                    septum_endo,
                )
            )
        epi_nodes = np.setdiff1d(epi_nodes, endo_nodes)

        endo_sid = self._add_nodeset(endo_nodes, "endocardium")
        epi_sid = self._add_nodeset(epi_nodes, "epicardium")

        # base-apical uvc
        # apex is selected only at left ventricle and with a region of 10 mm
        # This avoids mesh sensitivity and seems consistent with Strocchi paper's figure
        apex_nodes = self.model.get_apex_node_set(radius=self._UVC_APEX_RADIUS)
        apex_sid = self._add_nodeset(apex_nodes, "apex")

        # base is with all cap nodes
        (pv_nodes, tv_nodes, av_nodes, mv_nodes), _ = self._update_ventricular_caps_nodes()
        if isinstance(self.model, LeftVentricle):
            base_nodes = np.hstack((mv_nodes, av_nodes))
        else:
            base_nodes = np.hstack((mv_nodes, av_nodes, pv_nodes, tv_nodes))

        base_sid = self._add_nodeset(base_nodes, "base")

        # rotational uvc
        rot_start, rot_end, rot_mid = self._get_uvc_rotation_bc()

        sid_minus_pi = self._add_nodeset(rot_start, title="rotation:-pi")
        sid_plus_pi = self._add_nodeset(rot_end, title="rotation:pi")
        sid_zero = self._add_nodeset(rot_mid, title="rotation:0")

        cases = [
            (1, "transmural", [endo_sid, epi_sid], [0, 1]),
            (2, "apico-basal", [apex_sid, base_sid], [0, 1]),
            (3, "rotational", [sid_minus_pi, sid_plus_pi, sid_zero], [-np.pi, np.pi, 0]),
        ]
        for case_id, job_name, set_ids, bc_values in cases:
            self.add_case(case_id, job_name, set_ids, bc_values)

    def _get_uvc_rotation_bc(self):
        """Select node set on long axis plane."""
        mesh = copy.deepcopy(self.target)
        mesh["cell_ids"] = np.arange(0, mesh.n_cells, dtype=int)
        mesh["point_ids"] = np.arange(0, mesh.n_points, dtype=int)
        slice = mesh.slice(
            origin=self.model.l4cv_axis["center"], normal=self.model.l4cv_axis["normal"]
        )
        crinkled = mesh.extract_cells(np.unique(slice["cell_ids"]))
        free_wall_center, septum_center = crinkled.clip(
            origin=self.model.l2cv_axis["center"],
            normal=-self.model.l2cv_axis["normal"],
            crinkle=True,
            return_clipped=True,
        )

        rotation_mesh = mesh.remove_cells(free_wall_center["cell_ids"])
        LOGGER.info(f"{mesh.n_points - rotation_mesh.n_points} nodes are removed from clip.")

        vn = mesh.points[free_wall_center["point_ids"]] - self.model.l4cv_axis["center"]
        v0 = np.tile(self.model.l4cv_axis["normal"], (len(free_wall_center["point_ids"]), 1))

        dot = np.einsum("ij,ij->i", v0, vn)  # dot product row by row
        set1 = np.unique(free_wall_center["point_ids"][dot >= 0])  # -pi
        set2 = np.unique(free_wall_center["point_ids"][dot < 0])  # pi
        set3 = np.unique(
            np.setdiff1d(septum_center["point_ids"], free_wall_center["point_ids"])
        )  # 0

        return set1, set2, set3

    def _update_parts_materials_db(self):
        """Loop over parts defined in the model and creates keywords."""
        LOGGER.debug("Updating part keywords...")

        # add parts with a dataframe
        section_id = self.get_unique_section_id()

        # get list of cavities from model
        for part in self.model.parts:
            # part.pid = self.get_unique_part_id()
            # material ID = part ID
            part.mid = part.pid

            part_df = pd.DataFrame(
                {
                    "heading": [part.name],
                    "pid": [part.pid],
                    "secid": [section_id],
                    "mid": [0],
                    "tmid": [part.mid],
                }
            )
            part_kw = keywords.Part()
            part_kw.parts = part_df
            self.kw_database.parts.append(part_kw)

            # set up material
            self.kw_database.parts.append(
                keywords.MatThermalIsotropic(tmid=part.mid, tro=1e-9, hc=1, tc=1)
            )

        # set up section solid
        section_kw = keywords.SectionSolid(secid=section_id, elform=10)
        self.kw_database.parts.append(section_kw)

        return

    def _update_main_db(self):
        self.kw_database.main.append(keywords.ControlSolution(soln=1))
        self.kw_database.main.append(keywords.ControlThermalSolver(atype=0, ptype=0, solver=11))
        self.kw_database.main.append(keywords.DatabaseBinaryD3Plot(dt=1.0))
        self.kw_database.main.append(keywords.DatabaseGlstat(dt=1.0))
        self.kw_database.main.append(keywords.DatabaseMatsum(dt=1.0))
        self.kw_database.main.append(keywords.DatabaseTprint(dt=1.0))
        self.kw_database.main.append(keywords.DatabaseExtentBinary(therm=2))  # save heat flux
        self.kw_database.main.append(keywords.ControlTermination(endtim=1, dtmin=1.0))

    def _add_nodeset(self, nodes: np.ndarray, title: str, nodeset_id: int = None) -> int:
        """Convert to local node ID and add to nodeset.

        Parameters
        ----------
        nodes : np.ndarray
            Nodes global ids
        title : str
            nodeset title
        nodeset_id : int, optional
            attribute a nodeset ID if not given, by default None

        Returns
        -------
        int
            nodeset id
        """
        # get node IDs of sub mesh
        nodes = np.where(np.isin(self.target["point_ids"], nodes))[0]
        if nodeset_id is None:
            nodeset_id = self.get_unique_nodeset_id()
        # lsdyna ID start with 1
        kw = create_node_set_keyword(nodes + 1, node_set_id=nodeset_id, title=title)
        self.kw_database.node_sets.append(kw)
        return nodeset_id

    def _update_drbm_bc(self):
        """Update D-RBM boundary conditions."""

        def clean_node_set(nodes: np.ndarray, exclude_nodes: np.ndarray = None):
            """Make sure there are no duplicate or excluded nodes, avoid thermal BC error."""
            nodes = np.unique(nodes)
            if exclude_nodes is not None:
                nodes = np.setdiff1d(nodes, exclude_nodes)
            return nodes

        (pv_nodes, tv_nodes, av_nodes, mv_nodes), combined_av_mv = (
            self._update_ventricular_caps_nodes()
        )

        if isinstance(self.model, LeftVentricle):
            rings_nodes = np.hstack((mv_nodes, av_nodes))
        else:
            rings_nodes = np.hstack((mv_nodes, av_nodes, pv_nodes, tv_nodes))

        # LV endo
        lv_endo_nodes = self.model.left_ventricle.endocardium.global_node_ids_triangles
        lv_endo_nodes = clean_node_set(lv_endo_nodes, rings_nodes)
        # LV epi
        epi_nodes = self.model.left_ventricle.epicardium.global_node_ids_triangles
        epi_nodes = clean_node_set(epi_nodes, np.hstack((lv_endo_nodes, rings_nodes)))
        # LV apex
        la_node = self.model.get_apex_node_set(part="left")

        if not isinstance(self.model, LeftVentricle):
            # Right ventricle endocardium
            septum_endo = self._get_rv_septum_endo_surface()
            rv_endo_nodes = np.hstack(
                (
                    self.model.right_ventricle.endocardium.global_node_ids_triangles,
                    septum_endo.global_node_ids_triangles,
                )
            )
            rv_endo_nodes = clean_node_set(rv_endo_nodes, rings_nodes)

            # append RV epi
            epi_nodes = np.hstack(
                (
                    epi_nodes,
                    self.model.right_ventricle.epicardium.global_node_ids_triangles,
                )
            )
            epi_nodes = clean_node_set(epi_nodes, np.hstack((rv_endo_nodes, rings_nodes)))
            # RV apex
            ra_node = self.model.get_apex_node_set(part="right")

        if isinstance(self.model, LeftVentricle):
            lv_endo_nodeset_id = self._add_nodeset(lv_endo_nodes, "lv endo")
            epi_nodeset_id = self._add_nodeset(epi_nodes, "epi")
            mv_nodeset_id = self._add_nodeset(mv_nodes, "mv")
            av_nodeset_id = self._add_nodeset(av_nodes, "av")
            la_nodeset_id = self._add_nodeset(la_node, "left apex")

            # add case kewyords
            cases = [
                (1, "trans", [lv_endo_nodeset_id, epi_nodeset_id], [1, 0]),
                (2, "ab_l", [mv_nodeset_id, la_nodeset_id], [1, 0]),
                (3, "ot_l", [av_nodeset_id, la_nodeset_id], [1, 0]),
                # If combined MV and AV, mv_nodeset=av_nodeset=combined, solve ab_l = ot_l
                # w_l's has no effect on the result, so set only for structure of code
                (4, "w_l", [mv_nodeset_id, la_nodeset_id], [1, 0])
                if combined_av_mv
                else (4, "w_l", [mv_nodeset_id, la_nodeset_id, av_nodeset_id], [1, 1, 0]),
            ]
        elif isinstance(self.model, (FullHeart, FourChamber, BiVentricle)):
            lv_endo_nodeset_id = self._add_nodeset(lv_endo_nodes, "lv endo")
            rv_endo_nodeset_id = self._add_nodeset(rv_endo_nodes, "rv endo")
            epi_nodeset_id = self._add_nodeset(epi_nodes, "epi")
            mv_nodeset_id = self._add_nodeset(mv_nodes, "mv")
            av_nodeset_id = self._add_nodeset(av_nodes, "av")
            tv_nodeset_id = self._add_nodeset(tv_nodes, "tv")
            pv_nodeset_id = self._add_nodeset(pv_nodes, "pv")
            la_nodeset_id = self._add_nodeset(la_node, "left apex")
            ra_nodeset_id = self._add_nodeset(ra_node, "right apex")

            # add case kewyords
            cases = [
                (1, "trans", [lv_endo_nodeset_id, rv_endo_nodeset_id, epi_nodeset_id], [2, -1, 0]),
                (2, "ab_l", [mv_nodeset_id, la_nodeset_id], [1, 0]),
                (3, "ab_r", [tv_nodeset_id, ra_nodeset_id], [1, 0]),
                (4, "ot_l", [av_nodeset_id, la_nodeset_id], [1, 0]),
                (5, "ot_r", [pv_nodeset_id, ra_nodeset_id], [1, 0]),
                # If combined MV and AV, mv_nodeset=av_nodeset=combined, solve ab_l = ot_l
                # w_l's has no effect on the result, so set only for structure of code
                (6, "w_l", [mv_nodeset_id, la_nodeset_id], [1, 0])
                if combined_av_mv
                else (6, "w_l", [mv_nodeset_id, la_nodeset_id, av_nodeset_id], [1, 1, 0]),
                (7, "w_r", [tv_nodeset_id, ra_nodeset_id, pv_nodeset_id], [1, 1, 0]),
                (8, "lr", [lv_endo_nodeset_id, rv_endo_nodeset_id], [1, -1]),
            ]

        for case_id, job_name, set_ids, bc_values in cases:
            self.add_case(case_id, job_name, set_ids, bc_values)

    def _get_rv_septum_endo_surface(self):
        """Get the right ventricle septum endocardium surface."""
        for surface in self.model.right_ventricle.surfaces:
            if "endocardium" in surface.name and "septum" in surface.name:
                return surface

        raise ValueError("Septum endocardium surface not found in right ventricle.")

    def _update_ventricular_caps_nodes(self):
        combined_av_mv = False  # combined mitral and aortic valve
        mv_nodes = av_nodes = tv_nodes = pv_nodes = None

        for part in self.model.parts:
            for cap in part.caps:
                if cap.type == CapType.MITRAL_VALVE:
                    mv_nodes = cap.global_node_ids_edge
                if cap.type == CapType.AORTIC_VALVE:
                    av_nodes = cap.global_node_ids_edge
                if cap.type == CapType.COMBINED_MITRAL_AORTIC_VALVE:
                    mv_nodes = av_nodes = cap.global_node_ids_edge
                    combined_av_mv = True

                if not isinstance(self.model, LeftVentricle):
                    if cap.type == CapType.TRICUSPID_VALVE:
                        tv_nodes = cap.global_node_ids_edge
                    if cap.type == CapType.PULMONARY_VALVE:
                        pv_nodes = cap.global_node_ids_edge

        return (pv_nodes, tv_nodes, av_nodes, mv_nodes), combined_av_mv

    def add_case(self, case_id: int, case_name: str, set_ids: list[int], bc_values: list[float]):
        """Add case to keyword database.

        Parameters
        ----------
        case_id : int
            case id
        case_name : str
            case name, will be d3plot file name
        set_ids : list[int]
            node set id for boundary condition
        bc_values : list[float]
            boundary condition values
        """
        # declare case
        self.kw_database.main.append(keywords.Case(caseid=case_id, jobid=case_name, scid1=case_id))
        # define BC for this case
        self.kw_database.main.append(f"*CASE_BEGIN_{case_id}")
        for sid, value in zip(set_ids, bc_values):
            self.kw_database.main.append(
                keywords.BoundaryTemperatureSet(
                    nsid=sid,
                    lcid=0,
                    cmult=value,
                ),
            )
        self.kw_database.main.append(f"*CASE_END_{case_id}")


if __name__ == "__main__":
    print("protected")
    pass
