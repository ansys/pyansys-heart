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
"""Base module containing classes for writing LS-DYNA keyword files."""

import os
import shutil
import time
from typing import List

import numpy as np
import pandas as pd

from ansys.dyna.core.keywords import keywords
from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.models import BiVentricle, FourChamber, FullHeart, HeartModel, LeftVentricle
from ansys.health.heart.objects import SurfaceMesh
import ansys.health.heart.parts as anatomy
import ansys.health.heart.settings.settings as sett
from ansys.health.heart.settings.settings import SimulationSettings
from ansys.health.heart.writer import custom_keywords as custom_keywords
from ansys.health.heart.writer.heart_decks import BaseDecks, FiberGenerationDecks
from ansys.health.heart.writer.writer_utils import (
    add_nodes_to_kw,
    create_element_solid_ortho_keyword,
    create_node_set_keyword,
    create_segment_set_keyword,
    fast_element_writer,
    get_list_of_used_ids,
)


class BaseDynaWriter:
    """Base class that contains essential features for all LS-DYNA heart models."""

    def __init__(self, model: HeartModel, settings: SimulationSettings = None) -> None:
        """Initialize writer by loading a heart model and the desired settings.

        Parameters
        ----------
        model : HeartModel
            Object that contains the necessary information for the writer,
            such as nodes, elements, and parts.
        settings : SimulationSettings, default: None
            Simulation settings for creating the LS-DYNA model.
            The dfeault settings are loaded in ``None``is used.

        """
        self.model = model
        """Model information necessary for creating the LS-DYNA K files."""

        self.kw_database = BaseDecks()

        # These are general attributes useful for keeping track of IDs:
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
        """ID offset for several relevant keywords."""

        #! Do we really need the below?
        for part in self.model.parts:
            if not part.pid:
                part.pid = np.max([p.pid for p in self.model.parts if p.pid]) + 1

        self.id_offset["part"] = np.max(self.model.part_ids)

        # ! Removed the below since the part IDs in self.model.parts are already defined.
        # for part in self.model.parts:
        #     id += 1
        #     # cannot use get_unique_part_id() because it checks in Deck()
        #     # part.pid = self.get_unique_part_id()
        #     # part.pid = id
        # """Assign part id for heart parts."""

        if not settings:
            self.settings = SimulationSettings()
            """Simulation settings."""
            LOGGER.warning("No settings provided - loading default values.")
            self.settings.load_defaults()

        else:
            self.settings = settings
            """Simulation settings."""

        self.settings.to_consistent_unit_system()

        return

    def _get_subsettings(self) -> list[sett.Settings]:
        """Get subsettings from the settings object."""
        import ansys.health.heart.settings.settings as sett

        subsettings_classes = [
            getattr(self.settings, attr).__class__
            for attr in self.settings.__dict__
            if isinstance(getattr(self.settings, attr), sett.Settings)
        ]

        return subsettings_classes

    def _update_node_db(self, ids: np.ndarray = None) -> None:
        """Update node database.

        Parameters
        ----------
        ids : np.ndarray, default: None
            0-based IDs of the nodes to write.
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

    def _update_parts_db(self) -> None:
        """Loop over parts defined in the model and create keywords."""
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

    def _update_segmentsets_db(self, add_caps: bool = False, add_cavities: bool = True) -> None:
        """Update the segment set database."""
        # NOTE 0: add all surfaces as segment sets
        # NOTE 1: need to more robustly check segids that are already used?

        # add closed cavity segment sets
        if add_cavities:
            for cavity in self.model.cavities:
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
                    LOGGER.debug(f"Failed to create segment set for {surface.name}.")
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
            caps = self.model.all_caps
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

    def _filter_bc_nodes(self, surface: SurfaceMesh) -> np.ndarray:
        """Remove one or more nodes from tetrahedrons having all nodes in the boundary.

        Notes
        -----
        The removed node must be connected with at least one node outside the boundary. See #656.

        Parameters
        ----------
        surface : SurfaceMesh
            Boundary surface to analyze.

        Returns
        -------
        node_ids : np.ndarray
            Array of boundary nodes after problematic node removal.
        """
        # getting elements in active parts
        element_ids = np.array([], dtype=int)
        node_ids = surface.global_node_ids_triangles

        for part in self.model.parts:
            element_ids = np.append(element_ids, part.get_element_ids(self.model.mesh))

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

        # get tets with 4 nodes in boundary
        issue_tets = np.where(np.sum(tet_mask, axis=0) == 4)[0]

        # get corresponding nodes
        issue_nodes = active_tets[issue_tets, :]

        # count node appearances
        u_active_tets, tet_count_active = np.unique(active_tets, return_counts=True)
        u_issue_nodes, tet_count_issue = np.unique(issue_nodes, return_counts=True)

        # find issue nodes that belong to at least one non-issue tet
        removable_mask = np.array(
            [
                tet_count_active[np.where(u_active_tets == ii)[0][0]]
                != tet_count_issue[np.where(u_issue_nodes == ii)[0][0]]
                for ii in issue_nodes.flatten()
            ]
        ).reshape(-1, 4)

        # remove the first issue node belonging to at least one non-issue tet (for each tet)
        column_idxs = np.argmax(removable_mask, axis=1)
        nodes_toremove = np.unique(
            [issue_nodes[ii, column_idxs[ii]] for ii in range(len(issue_tets))]
        )

        # check that there are no nodes that only belong to non-issue tets
        if not np.all(np.any(removable_mask, axis=1)):
            # remove all such nodes and all their neighbors
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
    ) -> None:
        """Update the nodeset database.

        Parameters
        ----------
        remove_duplicates : bool, default: True
            Whether to remove nodes if they are used in other nodesets.
        remove_one_node_from_cell : bool, default: False
            Whether to remove a node if a cell has all nodes in a nodeset.

        Notes
        -----
        The ``FiberGenerationWriter`` module does  not allow all nodes of the same
        element in one nodeset.
        """
        # formats endo, epi- and septum nodeset keywords, do for all surfaces
        # for each surface in each part add the respective node-set
        # Use same ID as surface
        # TODO: check if database already contains nodesets (there will be duplicates otherwise)
        used_node_ids = np.empty(0, dtype=int)

        # add node-set for each cap
        parts_with_caps = [part for part in self.model.parts if isinstance(part, anatomy.Chamber)]
        for part in parts_with_caps:
            for cap in self.model.all_caps:
                # update cap mesh:
                cap._mesh = self.model.mesh.get_surface(cap._mesh.id)
                if remove_duplicates:
                    node_ids = np.setdiff1d(cap.global_node_ids_edge, used_node_ids)
                else:
                    node_ids = cap.global_node_ids_edge

                if len(node_ids) == 0:
                    LOGGER.debug(
                        "Nodes already used. Skipping nodeset for {0}".format(
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
                    LOGGER.debug(f"Failed to create nodeset for {surface.name}. Empty mesh.")
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
        """Get unique ID of a given keyword.

        Parameters
        ----------
        keyword : str
            Keyword string: valid inputs include:
            ["SECTION", "PART", "MAT", "SET_SEGMENT", "SET_NODE", "CURVE", ...]
        return_used_ids : bool, default: False
            Whether to return used IDs along with the next unique ID.

        Returns
        -------
        int
            Next unique ID.
        """
        used_ids = [0]
        for key in self.kw_database.__dict__.keys():
            db = self.kw_database.__dict__[key]
            used_ids = np.append(used_ids, get_list_of_used_ids(db, keyword))
        used_ids = np.array(used_ids, dtype=int)
        _, counts = np.unique(used_ids, return_counts=True)
        if np.any(counts > 1):
            raise ValueError("{0} Duplicate IDs found for: {1}".format(counts, keyword))

        if return_used_ids:
            return np.max(used_ids) + 1, used_ids
        else:
            return np.max(used_ids) + 1

    def get_unique_part_id(self) -> int:
        """Suggest a unique non-used part ID."""
        return self._get_unique_id("PART")

    def get_unique_mat_id(self) -> int:
        """Suggest a unique non-used material ID."""
        return self._get_unique_id("MAT")

    def get_unique_section_id(self) -> int:
        """Suggest a unique non-used section ID."""
        return self._get_unique_id("SECTION")

    def get_unique_segmentset_id(self) -> int:
        """Suggest a unique non-used segment set ID."""
        return self._get_unique_id("SET_SEGMENT")

    def get_unique_nodeset_id(self) -> int:
        """Suggest a unique non-used nodeset ID."""
        return self._get_unique_id("SET_NODE")

    def get_unique_partset_id(self) -> int:
        """Suggest a unique non-used part ID."""
        return self._get_unique_id("SET_PART")

    def get_unique_curve_id(self) -> int:
        """Suggest a unique curve ID."""
        return self._get_unique_id("DEFINE_CURVE")

    def _get_decknames_of_include(self) -> list[str]:
        """
        Get a list of deck file names in the keyword database.

        Do not get those in the main deck and omit any empty decks.
        """
        include_files = []
        for deckname, deck in vars(self.kw_database).items():
            if deckname == "main":
                continue
            # skip if no keywords are present in the deck
            if len(deck.keywords) == 0:
                LOGGER.debug("No keywords in deck: {0}".format(deckname))
                continue
            include_files.append(deckname + ".k")

        return include_files

    def include_to_main(self, file_list: list[str] | str = []) -> None:
        """Add *INCLUDE keywords into the main decl.

        Parameters
        ----------
        file_list : list[str] | str, default: []
            Files to include.
        """
        if isinstance(file_list, str):
            file_list = [file_list]

        for file in file_list:
            self.kw_database.main.append(keywords.Include(filename=file))

        return

    def export(self, export_directory: str, user_k: list[str] | None = None) -> None:
        """Write the model to files.

        Parameters
        ----------
        export_directory : str
            Export directory.
        user_k : list[str], default: None
            User-provided K files.
        """
        tstart = time.time()
        LOGGER.info("Writing all LS-DYNA K files...")

        if not os.path.isdir(export_directory):
            os.makedirs(export_directory)

        if user_k is not None:
            for k_file in user_k:
                if not os.path.isfile(k_file):
                    error_msg = f"File {k_file} is not found."
                    LOGGER.error(error_msg)
                    raise FileNotFoundError(error_msg)
                else:
                    name = os.path.basename(k_file)
                    shutil.copy(k_file, os.path.join(export_directory, name))
                    self.include_to_main(name)

        # export .k files
        self.export_databases(export_directory)

        # export settings
        self.settings.save(os.path.join(export_directory, "simulation_settings.yml"))

        tend = time.time()
        LOGGER.debug("Time spent writing files: {:.2f} s".format(tend - tstart))

        return

    def export_databases(self, export_directory: str) -> None:
        """Export each non-empty database to a specified directory."""
        if not export_directory:
            export_directory = self.model.info.working_directory

        for deckname, deck in vars(self.kw_database).items():
            # skip empty databases:
            if len(deck.keywords) == 0 and len(deck.string_keywords) == 0:
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

    def _keep_ventricles(self) -> None:
        """Remove any non-ventricular parts."""
        LOGGER.debug("Only keeping ventricular-parts for fiber/Purkinje generation.")
        # Note: we need to use _PartType to check part types.
        # For example, base has _PartType as VENTRICLE,
        # but it is not an instance from Ventricle(Part) so will be missed in writing.
        parts_to_keep = [
            p.name
            for p in self.model.parts
            if p._part_type in [anatomy._PartType.VENTRICLE, anatomy._PartType.SEPTUM]
        ]

        self._keep_parts(parts_to_keep)
        return

    def _keep_parts(self, parts_to_keep: List[str]) -> None:
        """Remove parts by a list of part names."""
        parts_to_remove = [part for part in self.model.part_names if part not in parts_to_keep]
        for part_to_remove in parts_to_remove:
            LOGGER.warning(f"Removing: {part_to_remove}")
            self.model.remove_part(part_to_remove)
        return

    def _update_solid_elements_db(self, add_fibers: bool = True) -> None:
        """
        Create slid (ortho) elements for all parts.

        Parameters
        ----------
        add_fibers: bool, default: True
            Whether to add fibers in general.
        """
        LOGGER.debug("Updating solid element keywords...")

        if add_fibers:
            cell_data_fields = self.model.mesh.cell_data.keys()
            if "fiber" not in cell_data_fields or "sheet" not in cell_data_fields:
                raise KeyError("Mechanics writer requires fiber and sheet fields.")

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
            tetrahedrons = (
                self.model.mesh.tetrahedrons[part.get_element_ids(self.model.mesh), :] + 1
            )
            num_elements = tetrahedrons.shape[0]

            # element_ids = np.arange(1, num_elements + 1, 1) + solid_element_count
            part_ids = np.ones(num_elements, dtype=int) * part.pid

            # format the element keywords
            if not part_add_fibers:
                kw_elements = keywords.ElementSolid()
                elements = pd.DataFrame(
                    {
                        "eid": part.get_element_ids(self.model.mesh) + 1,
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
                fiber = self.volume_mesh.cell_data["fiber"][part.get_element_ids(self.model.mesh)]
                sheet = self.volume_mesh.cell_data["sheet"][part.get_element_ids(self.model.mesh)]

                # normalize fiber and sheet directions:
                # norm = np.linalg.norm(fiber, axis=1)
                # fiber = fiber / norm[:, None]
                # norm = np.linalg.norm(sheet, axis=1)
                # sheet = sheet / norm[:, None]

                kw_elements = create_element_solid_ortho_keyword(
                    elements=tetrahedrons,
                    a_vec=fiber,
                    d_vec=sheet,
                    e_id=part.get_element_ids(self.model.mesh) + 1,
                    part_id=part_ids,
                    element_type="tetra",
                )

            # add elements to database
            self.kw_database.solid_elements.append(kw_elements)
            # solid_element_count = solid_element_count + num_elements

        return


class FiberGenerationDynaWriter(BaseDynaWriter):
    """Class for preparing the input for a fiber-generation LS-DYNA simulation."""

    def __init__(self, model: HeartModel, settings: SimulationSettings = None) -> None:
        super().__init__(model=model, settings=settings)
        self.kw_database = FiberGenerationDecks()
        """Collection of keywords relevant for fiber generation."""

        if sett.Fibers not in self._get_subsettings():
            raise ValueError("Expecting fiber settings.")

    def update(self, rotation_angles: dict[str, list[float]] | None = None) -> None:
        """Update keyword database for fiber generation.

        This method overwrites the inherited function.
        """
        ##
        self._update_main_db()  # needs updating

        if isinstance(self.model, (FourChamber, FullHeart)):
            LOGGER.warning(
                "Atrium are present in the model. These are removed for ventricle fiber generation."
            )

            parts = [
                part
                for part in self.model.parts
                if isinstance(part, (anatomy.Ventricle, anatomy.Septum))
            ]
            #! Note that this only works when tetrahedrons are added at the beginning
            #! of the mesh (file)! E.g. check self.mesh.celltypes to make sure this is the case!
            tet_ids = np.empty((0), dtype=int)
            for part in parts:
                tet_ids = np.append(tet_ids, part.get_element_ids(self.model.mesh))
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

        include_files = self._get_decknames_of_include()
        self.include_to_main(include_files)

        return

    def _remove_atrial_nodes_from_ventricles_surfaces(self) -> None:
        """Remove nodes other than ventricular from ventricular surfaces."""
        parts = [
            part
            for part in self.model.parts
            if isinstance(part, (anatomy.Ventricle, anatomy.Septum))
        ]

        tet_ids = np.empty((0), dtype=int)
        for part in parts:
            tet_ids = np.append(tet_ids, part.get_element_ids(self.model.mesh))
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

    def _update_material_db(self) -> None:
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
            # element_ids = part.get_element_ids(self.model.mesh)
            # em_mat_id = self.get_unique_mat_id()
            em_mat_id = part.mid  # Needs to match material id used in update_parts_db
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

    def _update_ep_settings(self) -> None:
        """Add the settings for the electrophysiology solver."""
        self.kw_database.ep_settings.append(
            keywords.EmControl(
                emsol=11, numls=4, macrodt=1, dimtype=None, nperio=None, ncylbem=None
            )
        )

        self.kw_database.ep_settings.append(keywords.EmControlTimestep(dtcons=1))

        # use defaults
        self.kw_database.ep_settings.append(custom_keywords.EmControlEp())

        # max iter should be int
        self.kw_database.ep_settings.append(
            keywords.EmSolverFem(reltol=1e-6, maxite=int(1e4), precon=2)
        )

        self.kw_database.ep_settings.append(keywords.EmOutput(mats=1, matf=1, sols=1, solf=1))

        return

    # TODO: Refactor
    def _update_create_fibers(self, rotation_angles: dict[str, list[float]]) -> None:
        """Update the keywords for fiber generation."""
        # collect relevant node and segment sets.
        # nodeset: apex, base
        # nodeset: endocardium, epicardium
        # NOTE: could be better if basal nodes are extracted in the preprocessor
        # since that would allow you to robustly extract these nodessets using the
        # input data
        # The below is relevant for all models.
        nodes_base = np.empty(0, dtype=int)
        node_sets_ids_endo = []  # relevant for both models
        node_sets_ids_epi = []  # relevant for both models
        node_set_ids_epi_and_rseptum = []  # only relevant for bv, 4c and full model

        # list of ventricular parts
        ventricles = [part for part in self.model.parts if isinstance(part, anatomy.Ventricle)]
        septum = next(
            (part for part in self.model.parts if isinstance(part, anatomy.Septum)),
            None,
        )

        # collect nodeset IDs (generated previously)
        node_sets_ids_epi = [ventricle.epicardium._node_set_id for ventricle in ventricles]
        node_sets_ids_endo = []
        for ventricle in ventricles:
            for surface in ventricle.surfaces:
                if "endocardium" in surface.name or "septum" in surface.name:
                    surf = self.model.mesh.get_surface(surface.id)
                    if surf.n_cells == 0:
                        LOGGER.debug(
                            f"Failed to collect nodeset ID for {surface.name}. Empty mesh."
                        )
                        continue
                    node_sets_ids_endo.append(surface._node_set_id)

        node_set_id_lv_endo = self.model.get_part("Left ventricle").endocardium._node_set_id
        if isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
            node_set_ids_epi_and_rseptum = node_sets_ids_epi + [
                self.model.right_ventricle.septum._node_set_id
            ]

        for cap in self.model.all_caps:
            nodes_base = np.append(nodes_base, cap.global_node_ids_edge)

        # apex ID [0] endocardium, [1] epicardium
        apex_point = self.model.get_part("Left ventricle").apex_points[1]
        if "epicardium" not in apex_point.name:
            raise ValueError("Expecting a point on the epicardium.")
        node_apex = apex_point.node_id  # is this a global node ID?

        # validate nodeset by removing nodes not part of the model without ventricles
        tet_ids_ventricles = np.empty((0), dtype=int)
        if septum:
            parts = ventricles + [septum]
        else:
            parts = ventricles

        for part in parts:
            tet_ids_ventricles = np.append(
                tet_ids_ventricles, part.get_element_ids(self.model.mesh)
            )

        tetra_ventricles = self.model.mesh.tetrahedrons[tet_ids_ventricles, :]

        # remove nodes that occur just in atrial part
        mask = np.isin(nodes_base, tetra_ventricles, invert=True)
        LOGGER.debug("Removing {0} nodes from base nodes...".format(np.sum(mask)))
        nodes_base = nodes_base[np.invert(mask)]

        # create set parts for lv and rv myocardium
        myocardium_part_ids = [ventricle.pid for ventricle in ventricles]

        # switch between the various models to generate valid input decks
        if isinstance(self.model, LeftVentricle):
            LOGGER.warning("Model type %s is in development. " % self.model.__class__.__name__)

            # Define part set for myocardium
            part_list1_kw = keywords.SetPartList(
                sid=1,
            )
            part_list1_kw.parts._data = myocardium_part_ids
            part_list1_kw.options["TITLE"].active = True
            part_list1_kw.title = "myocardium_all"

            self.kw_database.create_fiber.extend([part_list1_kw])

            # combine nodesets endocardium uing *SET_NODE_ADD:
            node_set_id_all_endocardium = self.get_unique_nodeset_id()

            set_add_kw = keywords.SetNodeAdd(sid=node_set_id_all_endocardium)
            set_add_kw.options["TITLE"].active = True
            set_add_kw.title = "all_endocardium_segments"
            set_add_kw.nodes._data = node_sets_ids_endo

            self.kw_database.create_fiber.append(set_add_kw)

            # combine nodesets epicardium:
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
            from ansys.health.heart.writer.define_function_templates import (
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

            # combine nodesets endocardium uing *SET_SEGMENT_ADD:
            node_set_id_all_endocardium = self.get_unique_nodeset_id()
            set_add_kw = keywords.SetNodeAdd(sid=node_set_id_all_endocardium)

            set_add_kw.options["TITLE"].active = True
            set_add_kw.title = "all_endocardium_segments"
            set_add_kw.nodes._data = node_sets_ids_endo

            self.kw_database.create_fiber.append(set_add_kw)

            # combine nodesets epicardium:
            node_set_id_all_epicardium = self.get_unique_nodeset_id()
            set_add_kw = keywords.SetNodeAdd(sid=node_set_id_all_epicardium)

            set_add_kw.options["TITLE"].active = True
            set_add_kw.title = "all_epicardium_segments"
            set_add_kw.nodes._data = node_sets_ids_epi

            self.kw_database.create_fiber.append(set_add_kw)

            # combine nodesets epicardium and septum:
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
            from ansys.health.heart.writer.define_function_templates import (
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

    def _update_main_db(self) -> None:
        self.kw_database.main.append(
            keywords.ControlTimeStep(dtinit=1.0, dt2ms=1.0, emscl=None, ihdo=None, rmscl=None)
        )

        self.kw_database.main.append(keywords.ControlTermination(endtim=10))

        self.kw_database.main.append(keywords.DatabaseBinaryD3Plot(dt=1.0))

        return
