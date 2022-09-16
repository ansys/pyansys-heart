"""Module contain. classes for writing LS-DYNA keywords based.

Note
----
Uses a HeartModel (from ansys.heart.preprocessor.models).

"""
import json
import os
from pathlib import Path
import time
from typing import List

from ansys.dyna.keywords import keywords
from ansys.heart.custom_logging import LOGGER
from ansys.heart.preprocessor.mesh.objects import Cap
import ansys.heart.preprocessor.mesh.vtkmethods as vtkmethods
from ansys.heart.preprocessor.models import (
    BiVentricle,
    FourChamber,
    FullHeart,
    HeartModel,
    LeftVentricle,
)

# import missing keywords
# import missing keywords
from ansys.heart.writer import custom_dynalib_keywords as custom_keywords
from ansys.heart.writer.heart_decks import (
    BaseDecks,
    ElectrophysiologyDecks,
    FiberGenerationDecks,
    MechanicsDecks,
    PurkinjeGenerationDecks,
)
from ansys.heart.writer.keyword_module import (
    add_nodes_to_kw,
    create_define_curve_kw,
    create_define_sd_orientation_kw,
    create_discrete_elements_kw,
    create_element_shell_keyword,
    create_element_solid_ortho_keyword,
    create_node_set_keyword,
    create_segment_set_keyword,
    fast_element_writer,
    get_list_of_used_ids,
)
from ansys.heart.writer.material_keywords import MaterialAtrium, MaterialHGOMyocardium, active_curve
import numpy as np
import pandas as pd
from vtk.numpy_interface import dataset_adapter as dsa  # noqa

# import commonly used material models


class BaseDynaWriter:
    """Base class that contains essential features for all LS-DYNA heart models."""

    def __init__(self, model: HeartModel) -> None:
        """Initialize writer by loading a HeartModel.

        Parameters
        ----------
        model : HeartModel
            HeartModel object which contains the necessary
            information for the writer, such as nodes and elements.
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

        """List of .k files to include in main. This is derived from the Decks classes."""
        self.include_files = []

        """Load simulation parameters."""
        from ansys.heart.writer.parameters import parameters

        self.parameters = parameters

        if "Improved" in self.model.info.model_type:
            LOGGER.warning(
                "Changing model type from : {0} to {1}".format(
                    self.model.info.model_type, self.model.info.model_type.replace("Improved", "")
                )
            )
            self.model.info.model_type = self.model.info.model_type.replace("Improved", "")

        return

    def _update_node_db(self):
        """Add nodes to the Node database."""
        LOGGER.debug("Updating node keywords...")
        node_kw = keywords.Node()
        node_kw = add_nodes_to_kw(self.model.mesh.nodes, node_kw)

        self.kw_database.nodes.append(node_kw)

        return

    def _update_parts_db(self):
        """Loop over parts defined in the model and creates keywords."""
        LOGGER.debug("Updating part keywords...")
        # add parts with a dataframe

        section_id = self.get_unique_section_id()
        # get list of cavities from model
        for part in self.model.parts:
            mat_id = self.get_unique_mat_id()
            # for element_set in cavity.element_sets:
            part_id = self.get_unique_part_id()
            part_name = part.name
            part_df = pd.DataFrame(
                {"heading": [part_name], "pid": [part_id], "secid": [section_id], "mid": [mat_id]}
            )
            part_kw = keywords.Part()
            part_kw.parts = part_df

            self.kw_database.parts.append(part_kw)

            # store part id for future use
            part.pid = part_id
            part.mid = mat_id

        # set up section solid for cavity myocardium
        section_kw = keywords.SectionSolid(secid=section_id, elform=13)

        self.kw_database.parts.append(section_kw)

        return

    def _update_solid_elements_db(self):
        """Create Solid ortho elements for all cavities.

        Note
        ----
        Each cavity contains one myocardium.

        """
        LOGGER.debug("Updating solid element keywords...")

        # create elements for each part
        solid_element_count = 0  # keeps track of number of solid elements already defined

        for part in self.model.parts:
            tetrahedrons = self.model.mesh.tetrahedrons[part.element_ids, :] + 1
            num_elements = tetrahedrons.shape[0]

            element_ids = np.arange(1, num_elements + 1, 1) + solid_element_count
            part_ids = np.ones(num_elements, dtype=int) * part.pid

            # format the element keywords
            kw_elements = keywords.ElementSolid()
            elements = pd.DataFrame(
                {
                    "eid": element_ids,
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
            # add elements to database
            self.kw_database.solid_elements.append(kw_elements)
            solid_element_count = solid_element_count + num_elements

        return

    def _update_segmentsets_db(self):
        """Update the segment set database."""
        # NOTE 0: add all surfaces as segment sets
        # NOTE 1: need to more robustly check segids that are already used?

        # add closed cavity segment sets
        cavities = [p.cavity for p in self.model.parts if p.cavity]
        for cavity in cavities:
            surface_id = self.get_unique_segmentset_id()
            cavity.surface.id = surface_id
            kw = create_segment_set_keyword(
                segments=cavity.surface.faces + 1,
                segid=cavity.surface.id,
                title=cavity.name,
            )
            # append this kw to the segment set database
            self.kw_database.segment_sets.append(kw)

        # write surfaces as segment sets
        for part in self.model.parts:
            for surface in part.surfaces:
                surface.id = self.get_unique_segmentset_id()
                kw = create_segment_set_keyword(
                    segments=surface.faces + 1,
                    segid=surface.id,
                    title=surface.name,
                )
                # append this kw to the segment set database
                self.kw_database.segment_sets.append(kw)

        return

    def _update_nodesets_db(self, remove_duplicates: bool = True):
        """Update the node set database."""
        # formats endo, epi- and septum nodeset keywords. Do for all surfaces and caps

        surface_ids = [s.id for p in self.model.parts for s in p.surfaces]
        node_set_id = np.max(surface_ids) + 1

        # for each surface in each part add the respective node-set
        # Use same ID as surface
        used_node_ids = np.empty(0, dtype=int)

        for part in self.model.parts:
            kws_surface = []
            for surface in part.surfaces:
                if remove_duplicates:
                    node_ids = np.setdiff1d(surface.node_ids, used_node_ids)
                else:
                    node_ids = surface.node_ids

                kw = create_node_set_keyword(
                    node_ids + 1, node_set_id=surface.id, title=surface.name
                )
                surface.nsid = surface.id
                kws_surface.append(kw)

                used_node_ids = np.append(used_node_ids, node_ids)

            self.kw_database.node_sets.extend(kws_surface)

    def _get_unique_id(self, keyword: str, return_used_ids: bool = False) -> int:
        """Get unique id of given keyword.

        Parameters
        ----------
        keyword : str
            Keyword string: valid inputs include:
            ["SECTION", "PART", "MAT", "SET_SEGMENT", "SET_NODE", "CURVE"]

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
            # use fast element writer for solid ortho elements
            if deckname == "solid_elements":

                element_kws = deck.get_kwds_by_type("ELEMENT")
                if os.path.isfile(filepath):
                    os.remove(filepath)

                for element_kw in element_kws:
                    fast_element_writer(element_kw, filepath)

                fid = open(filepath, "a")
                fid.write("*END")

            elif deckname == "nodes":
                ids = np.arange(0, self.model.mesh.nodes.shape[0], 1) + 1
                content = np.hstack((ids.reshape(-1, 1), self.model.mesh.nodes))
                np.savetxt(
                    os.path.join(export_directory, "nodes.k"),
                    content,
                    fmt="%8d%16.5e%16.5e%16.5e",
                    header="*KEYWORD\n*NODE\n"
                    "$#   nid               x               y               z      tc      rc",
                    footer="*END",
                    comments="",
                )
            else:
                deck.export_file(filepath)
        return

    def _export_cavity_segmentsets(self, export_directory: str):
        """Export the cavity segment sets to separate files."""
        cavities = [part.cavity for part in self.model.parts if part.cavity]
        for cavity in cavities:
            filepath = os.path.join(
                export_directory, "_".join(cavity.name.lower().split()) + ".segment"
            )
            np.savetxt(filepath, cavity.surface.faces + 1, delimiter=",", fmt="%d")

        return

    def _keep_ventricles(self):
        """Remove any non-ventricular parts."""
        # NOTE: Could move "remove part" method to model
        LOGGER.warning("Just keeping ventricular-parts for fiber generation")
        parts_to_keep = ["Left ventricle", "Right ventricle", "Septum"]
        parts_to_remove = [part for part in self.model.part_names if part not in parts_to_keep]
        for part_to_remove in parts_to_remove:
            LOGGER.warning("Removing: {}".format(part_to_remove))
            self.model.remove_part(part_to_remove)
        return

    def get_apex_left(self):
        """Get apex of left ventricle cavity."""
        # collect relevant node and segment sets.
        # node set: apex, base
        # node set: endocardium, epicardium
        # NOTE: could be better if basal nodes are extracted in the preprocessor
        # since that would allow you to robustly extract these nodessets using the
        # input data
        # The below is relevant for all models.
        node_apex_left = np.empty(0, dtype=int)

        # apex_points[0]: endocardium, apex_points[1]: epicardium
        if isinstance(self.model, (LeftVentricle, BiVentricle, FourChamber, FullHeart)):
            node_apex_left = self.model.left_ventricle.apex_points[0].node_id

            # check whether point is on edge of endocardium - otherwise pick another node in
            # the same triangle
            endocardium = self.model.left_ventricle.endocardium
            endocardium.get_boundary_edges()
            if np.any(endocardium.boundary_edges == node_apex_left):
                element_id = np.argwhere(np.any(endocardium.faces == node_apex_left, axis=1))[0][0]

                node_apex_left = endocardium.faces[element_id, :][
                    np.argwhere(
                        np.isin(
                            endocardium.faces[element_id, :],
                            endocardium.boundary_edges,
                            invert=True,
                        )
                    )[0][0]
                ]
                LOGGER.warning(
                    "Node id {0} is on edge of {1}. Picking node id {2}".format(
                        self.model.left_ventricle.apex_points[0].node_id,
                        endocardium.name,
                        node_apex_left,
                    )
                )
                self.model.left_ventricle.apex_points[0].node_id = node_apex_left
        return node_apex_left

    def get_apex_right(self):
        """Get apex of right ventricle cavity."""
        node_apex_right = np.empty(0, dtype=int)
        if isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
            node_apex_right = self.model.right_ventricle.apex_points[0].node_id

            # check whether point is on edge of endocardium - otherwise pick another node in
            # the same triangle
            endocardium = self.model.right_ventricle.endocardium
            endocardium.get_boundary_edges()
            if np.any(endocardium.boundary_edges == node_apex_right):
                element_id = np.argwhere(np.any(endocardium.faces == node_apex_right, axis=1))[0][0]

                node_apex_right = endocardium.faces[element_id, :][
                    np.argwhere(
                        np.isin(
                            endocardium.faces[element_id, :],
                            endocardium.boundary_edges,
                            invert=True,
                        )
                    )[0][0]
                ]
                LOGGER.warning(
                    "Node id {0} is on edge of {1}. Picking node id {2}".format(
                        self.model.right_ventricle.apex_points[0].node_id,
                        endocardium.name,
                        node_apex_right,
                    )
                )
                self.model.right_ventricle.apex_points[0].node_id = node_apex_right
        return node_apex_right


class MechanicsDynaWriter(BaseDynaWriter):
    """Class for preparing the input for a mechanics LS-DYNA simulation."""

    def __init__(self, model: HeartModel, system_model_name: str = "ClosedLoop") -> None:
        super().__init__(model)

        self.kw_database = MechanicsDecks()
        """Collection of keyword decks relevant for mechanics."""

        self.system_model_name = system_model_name
        """Name of system model to use."""

        # Depending on the system model specified give list of parameters

        return

    @property
    def system_model_name(self):
        """System model name.

        Note
        ----
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

    def update(self):
        """Update the keyword database."""
        self._update_main_db()
        self._update_node_db()
        self._update_parts_db()
        self._update_solid_elements_db(add_fibers=True)
        self._update_segmentsets_db()
        self._update_nodesets_db()
        self._update_material_db(add_active=True)

        # for boundary conditions
        self._add_cap_bc(bc_type="springs_caps")
        self._add_pericardium_bc()

        # # for control volume
        self._update_cap_elements_db()
        self._update_controlvolume_db()
        self._update_system_model()

        self._get_list_of_includes()
        self._add_includes()

        return

    def export(self, export_directory: str):
        """Write the model to files."""
        tstart = time.time()
        LOGGER.debug("Writing all LS-DYNA .k files...")

        if not export_directory:
            export_directory = os.path.join(self.model.info.workdir, "mechanics")

        if not os.path.isdir(export_directory):
            os.makedirs(export_directory)

        # export .k files
        self.export_databases(export_directory)

        # add system json in case of closed loop. For open-loop this is already
        # added in the control volume database
        if self.system_model_name == "ClosedLoop":
            # exports system model
            path_system_model_settings = os.path.join(
                export_directory, "system_model_settings.json"
            )
            with open(path_system_model_settings, "w") as outfile:
                json.dump(self.system_model_json, indent=4, fp=outfile)

        # export segment sets to separate file
        self._export_cavity_segmentsets(export_directory)

        tend = time.time()
        LOGGER.debug("Time spent writing files: {:.2f} s".format(tend - tstart))

        return

    def _update_main_db(self, add_damping: bool = True):
        """Update the main .k file.

        Note
        ----
        Consider using a settings (json?) file as input.

        """
        LOGGER.debug("Updating main keywords...")

        self.kw_database.main.title = self.model.model_type

        if self.__class__.__name__ == "ZeroPressureMechanicsDynaWriter":
            scale_time = 1
            if self.parameters["Unit"]["Time"] == "ms":
                scale_time = 1000

            self._add_solution_controls(scale_time=1 * scale_time)
            self._add_export_controls(dt_output_d3plot=0.5 * scale_time)

        elif self.__class__.__name__ == "MechanicsDynaWriter":
            self._add_solution_controls(
                end_time=self.parameters["Time"]["End Time"],
                dtmin=self.parameters["Time"]["dtmin"],
                dtmax=self.parameters["Time"]["dtmax"],
            )
            self._add_export_controls(
                dt_output_d3plot=self.parameters["Time"]["dt_d3plot"],
                dt_output_icvout=self.parameters["Time"]["dt_icvout"],
            )

        if add_damping:
            self._add_damping()

        return

    def _update_node_db(self):
        """Add nodes to the NODE database."""
        LOGGER.debug("Updating node keywords...")
        node_kw = keywords.Node()
        node_kw = add_nodes_to_kw(self.model.mesh.nodes, node_kw)

        self.kw_database.nodes.append(node_kw)

        return

    def _update_parts_db(self):
        """Loop over parts defined in the model and create keywords."""
        LOGGER.debug("Updating part keywords...")
        # add parts with a dataframe

        section_id = self.get_unique_section_id()
        # get list of cavities from model
        for part in self.model.parts:
            mat_id = self.get_unique_mat_id()
            # for element_set in cavity.element_sets:
            part_id = self.get_unique_part_id()
            part_name = part.name
            part_df = pd.DataFrame(
                {"heading": [part_name], "pid": [part_id], "secid": [section_id], "mid": [mat_id]}
            )
            part_kw = keywords.Part()
            part_kw.parts = part_df

            self.kw_database.parts.append(part_kw)

            # store part id for future use
            part.pid = part_id
            part.mid = mat_id

        # set up section solid for cavity myocardium
        section_kw = keywords.SectionSolid(secid=section_id, elform=13)

        self.kw_database.parts.append(section_kw)

        return

    def _update_solid_elements_db(self, add_fibers: bool = True):
        """Create Solid ortho elements for all cavities.

        Note
        ----
        Each cavity contains one myocardium.

        """
        LOGGER.debug("Updating solid element keywords...")

        cell_data_fields = self.model.mesh.cell_data.keys()
        if "fiber" not in cell_data_fields or "sheet" not in cell_data_fields:
            raise KeyError("Mechanics writer requires fiber and sheet fields")
            # logger.warning("Not writing fiber and sheet directions")
            # add_fibers = False

        # create elements for each part
        solid_element_count = 0  # keeps track of number of solid elements already defined

        for part in self.model.parts:
            if type(self) == MechanicsDynaWriter:
                if "ventricle" in part.name.lower() or "septum" in part.name.lower():
                    add_fibers = True
                else:
                    add_fibers = False

            LOGGER.debug(
                "\tAdding elements for {0} | adding fibers: {1}".format(part.name, add_fibers)
            )

            tetrahedrons = self.model.mesh.tetrahedrons[part.element_ids, :] + 1
            num_elements = tetrahedrons.shape[0]

            element_ids = np.arange(1, num_elements + 1, 1) + solid_element_count
            part_ids = np.ones(num_elements, dtype=int) * part.pid

            # format the element keywords
            if not add_fibers:
                kw_elements = keywords.ElementSolid()
                elements = pd.DataFrame(
                    {
                        "eid": element_ids,
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

            elif add_fibers:
                fiber = self.volume_mesh.cell_data["fiber"][part.element_ids]
                sheet = self.volume_mesh.cell_data["sheet"][part.element_ids]

                # normalize fiber and sheet directions:
                norm = np.linalg.norm(fiber, axis=1)
                fiber = fiber / norm[:, None]
                norm = np.linalg.norm(sheet, axis=1)
                sheet = sheet / norm[:, None]
                kw_elements = create_element_solid_ortho_keyword(
                    elements=tetrahedrons,
                    a_vec=fiber,
                    d_vec=sheet,
                    partid=part.pid,
                    id_offset=solid_element_count,
                    element_type="tetra",
                )

            # add elements to database
            self.kw_database.solid_elements.append(kw_elements)
            solid_element_count = solid_element_count + num_elements

        return

    def _add_solution_controls(
        self,
        end_time: float = 5,
        dtmin: float = 0.001,
        dtmax: float = 0.01,
        simulation_type: str = "quasi-static",
    ):
        """Add solution controls, output controls and solver settings."""
        # add termination keywords
        self.kw_database.main.append(keywords.ControlTermination(endtim=end_time, dtmin=dtmin))

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
                "Simulation type not recognized: Please choose " "either quasi-static or static"
            )

        prefill_time = self.parameters["Material"]["Myocardium"]["Active"]["Prefill"]
        self.kw_database.main.append(
            keywords.ControlImplicitDynamics(
                imass=imass,
                gamma=gamma,
                beta=beta,
                # active dynamic process only after prefilling
                tdybir=prefill_time,
            )
        )

        # add auto controls
        lcid = self.get_unique_curve_id()
        # tune time step for better compromise between convergence and performance
        time = [0, prefill_time, prefill_time + dtmax, end_time]
        step = [5 * dtmax, 5 * dtmax, dtmin, dtmax]
        kw_curve = create_define_curve_kw(
            x=time,
            y=step,
            curve_name="time step control",
            curve_id=lcid,
            lcint=0,
        )
        self.kw_database.main.append(kw_curve)
        self.kw_database.main.append(
            keywords.ControlImplicitAuto(iauto=1, dtmin=dtmin, dtmax=-lcid)
        )

        # add general implicit controls
        self.kw_database.main.append(
            keywords.ControlImplicitGeneral(imflag=1, dt0=dtmax)
        )  # imflag=1 means implicit

        # add implicit solution controls: Defaults are OK?
        self.kw_database.main.append(keywords.ControlImplicitSolution())

        # add implicit solver controls
        self.kw_database.main.append(custom_keywords.ControlImplicitSolver())
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

        self.kw_database.main.append(keywords.DatabaseElout(dt=0.1, binary=2))

        self.kw_database.main.append(keywords.DatabaseGlstat(dt=0.1, binary=2))

        self.kw_database.main.append(keywords.DatabaseMatsum(dt=0.1, binary=2))

        # frequency of full results
        lcid = self.get_unique_curve_id()
        time = [
            0,
            self.parameters["Material"]["Myocardium"]["Active"]["Prefill"],
            self.parameters["Material"]["Myocardium"]["Active"]["Prefill"] + dt_output_d3plot,
            self.parameters["Time"]["End Time"],
        ]
        step = [10 * dt_output_d3plot, 10 * dt_output_d3plot, dt_output_d3plot, dt_output_d3plot]
        kw_curve = create_define_curve_kw(
            x=time,
            y=step,
            curve_name="d3plot out control",
            curve_id=lcid,
            lcint=0,
        )

        self.kw_database.main.append(kw_curve)
        self.kw_database.main.append(
            keywords.DatabaseBinaryD3Plot(dt=dt_output_d3plot, lcdt=lcid, ioopt=1)
        )

        self.kw_database.main.append(keywords.DatabaseExtentBinary(neiph=27, strflg=1, maxint=0))

        return

    def _add_damping(self):
        """Add damping to the main file."""
        lcid_damp = self.get_unique_curve_id()

        kw_damp = keywords.DampingGlobal(lcid=lcid_damp)

        kw_damp_curve = create_define_curve_kw(
            x=[0, 10e25],  # to create a constant curve
            y=self.parameters["Global damping"] * np.array([1, 1]),
            curve_name="damping",
            curve_id=lcid_damp,
            lcint=0,
        )
        self.kw_database.main.append(kw_damp)
        self.kw_database.main.append(kw_damp_curve)

        return

    def _update_segmentsets_db(self):
        """Update the segment set database."""
        # NOTE 0: add all surfaces as segment sets
        # NOTE 1: need to more robustly check segids that are already used?

        # add closed cavity segment sets
        cavities = [p.cavity for p in self.model.parts if p.cavity]
        for cavity in cavities:
            surface_id = self.get_unique_segmentset_id()
            cavity.surface.id = surface_id
            kw = create_segment_set_keyword(
                segments=cavity.surface.faces + 1,
                segid=cavity.surface.id,
                title=cavity.name,
            )
            # append this kw to the segment set database
            self.kw_database.segment_sets.append(kw)

        # write surfaces as segment sets
        for part in self.model.parts:
            for surface in part.surfaces:
                surface.id = self.get_unique_segmentset_id()
                kw = create_segment_set_keyword(
                    segments=surface.faces + 1,
                    segid=surface.id,
                    title=surface.name,
                )
                # append this kw to the segment set database
                self.kw_database.segment_sets.append(kw)

        # create corresponding segment sets. Store in new file?
        caps = [cap for part in self.model.parts for cap in part.caps]
        for cap in caps:
            segid = self.get_unique_segmentset_id()
            setattr(cap, "seg_id", segid)
            segset_kw = create_segment_set_keyword(
                segments=cap.triangles + 1,
                segid=cap.seg_id,
                title=cap.name,
            )
            self.kw_database.segment_sets.append(segset_kw)
        return

    def _update_nodesets_db(self, remove_duplicates: bool = True):
        """Update the node set database."""
        # formats endo, epi- and septum nodeset keywords. Do for all surfaces and caps

        surface_ids = [s.id for p in self.model.parts for s in p.surfaces]
        node_set_id = np.max(surface_ids) + 1

        # for each surface in each part add the respective node-set
        # Use same ID as surface
        used_node_ids = np.empty(0, dtype=int)
        for part in self.model.parts:
            # add node-set for each cap
            kws_caps = []
            for cap in part.caps:
                if remove_duplicates:
                    node_ids = np.setdiff1d(cap.node_ids, used_node_ids)
                else:
                    node_ids = cap.node_ids

                cap.nsid = node_set_id
                kw = create_node_set_keyword(node_ids + 1, node_set_id=cap.nsid, title=cap.name)
                kws_caps.append(kw)
                node_set_id = node_set_id + 1

                used_node_ids = np.append(used_node_ids, node_ids)

            self.kw_database.node_sets.extend(kws_caps)

        for part in self.model.parts:
            kws_surface = []
            for surface in part.surfaces:
                if remove_duplicates:
                    node_ids = np.setdiff1d(surface.node_ids, used_node_ids)
                else:
                    node_ids = surface.node_ids

                kw = create_node_set_keyword(
                    node_ids + 1, node_set_id=surface.id, title=surface.name
                )
                surface.nsid = surface.id
                kws_surface.append(kw)

                used_node_ids = np.append(used_node_ids, node_ids)

            self.kw_database.node_sets.extend(kws_surface)

    def _update_material_db(self, add_active: bool = True):
        """Update the database of material keywords."""
        act_curve_id = self.get_unique_curve_id()

        for part in self.model.parts:
            part.mid = part.pid
            mat_id = part.mid

            if "ventricle" in part.name.lower() or "septum" in part.name.lower():
                if not add_active:
                    active_dct = None
                else:
                    active_dct = {
                        "taumax": self.parameters["Material"]["Myocardium"]["Active"]["Tmax"],
                        "ca2ionm": self.parameters["Material"]["Myocardium"]["Active"]["ca2ionm"],
                    }

                myocardium_kw = MaterialHGOMyocardium(
                    mid=part.mid,
                    iso_user=self.parameters["Material"]["Myocardium"]["Isotropic"],
                    anisotropy_user=self.parameters["Material"]["Myocardium"]["Anisotropic"],
                    active_user=active_dct,
                )

                myocardium_kw.acid = act_curve_id

                self.kw_database.material.append(myocardium_kw)

            elif "atrium" in part.name:
                # add atrium material
                # atrium_kw = MaterialAtrium(mid=part.mid)
                atrium_kw = MaterialHGOMyocardium(
                    mid=part.mid, iso_user=self.parameters["Material"]["Atrium"]
                )

                self.kw_database.material.append(atrium_kw)

            else:
                LOGGER.warning("Assuming same material as atrium for: {0}".format(part.name))

                # general_tissue_kw = MaterialAtrium(mid=part.mid)
                general_tissue_kw = MaterialHGOMyocardium(
                    mid=part.mid, iso_user=self.parameters["Material"]["Atrium"]
                )
                self.kw_database.material.append(general_tissue_kw)

        if add_active:
            # write and add active curve to material database
            time_array, active_stress_array = active_curve("Strocchi2020")

            if self.parameters["Unit"]["Time"] == "ms":
                time_array *= 1000

            active_curve_kw = create_define_curve_kw(
                x=time_array,
                y=active_stress_array,
                curve_name="calcium_concentration",
                curve_id=act_curve_id,
                lcint=15000,
            )

            # y scaling
            active_curve_kw.sfo = self.parameters["Material"]["Myocardium"]["Active"]["ca2ionm"]
            # x offset: prefill duration
            active_curve_kw.offa = self.parameters["Material"]["Myocardium"]["Active"]["Prefill"]
            self.kw_database.material.append(active_curve_kw)

        return

    def _add_cap_bc(self, bc_type: str):
        """Add boundary condition to the cap.

        Parameters
        ----------
        bc_type : str
            Boundary condition type. Valid bc's include: ["fix_caps", "springs_caps"].
        """
        valid_bcs = ["fix_caps", "springs_caps"]
        if bc_type not in valid_bcs:
            raise ValueError("Cap/Valve boundary condition must be of type: %r" % valid_bcs)

        # create list of cap names where to add the spring b.c
        caps_to_use = []
        if isinstance(self.model, (LeftVentricle, BiVentricle)):
            # use all caps:
            # for cavity in self.model._mesh._cavities:
            #     for cap in cavity.closing_caps:
            #         caps_to_use.append(cap.name)

            caps_to_use = [
                "mitral-valve",
                "tricuspid-valve",
                "aortic-valve",
                "pulmonary-valve",
            ]

        elif isinstance(self.model, (FourChamber, FullHeart)):
            caps_to_use = [
                "superior-vena-cava",
                "right-inferior-pulmonary-vein",
                "right-superior-pulmonary-vein",
            ]

        if bc_type == "fix_caps":
            for part in self.model.parts:
                for cap in part.caps:
                    if cap.name in caps_to_use:
                        kw_fix = keywords.BoundarySpcSet()
                        kw_fix.nsid = cap.nsid
                        kw_fix.dofx = 1
                        kw_fix.dofy = 1
                        kw_fix.dofz = 1

                        self.kw_database.boundary_conditions.append(kw_fix)

        # if bc type is springs -> add springs
        # NOTE add to boundary condition db or separate spring db?
        elif bc_type == "springs_caps":

            part_id = self.get_unique_part_id()
            section_id = self.get_unique_section_id()
            mat_id = self.get_unique_mat_id()

            if isinstance(self.model, (LeftVentricle, BiVentricle)):
                spring_stiffness = self.parameters["Boundary Condition"]["Valve Spring"][
                    "BV"
                ]  # kPa/mm

            elif isinstance(self.model, (FourChamber, FullHeart)):
                spring_stiffness = self.parameters["Boundary Condition"]["Valve Spring"][
                    "4C"
                ]  # kPa/mm

            scale_factor_normal = self.parameters["Boundary Condition"]["Normal Scale factor"]
            scale_factor_radial = self.parameters["Boundary Condition"]["Radial Scale factor"]

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
                if cap.name in caps_to_use:
                    self._add_springs_cap_edge(
                        cap,
                        part_id,
                        scale_factor_normal,
                        scale_factor_radial,
                    )

        return

    def _add_springs_cap_edge(
        self,
        cap: Cap,
        part_id: int,
        scale_factor_normal: float,
        scale_factor_radial: float,
    ):
        """Add springs to the cap nodes.

        Note
        ----
        Appends these to the boundary condition database.
        """
        # -------------------------------------------------------------------
        LOGGER.debug("Adding spring b.c. for cap: %s" % cap.name)

        # NOTE: may want to extent the node ids to include adjacent nodes
        # num_nodes_edge = len(cap.node_ids)
        mesh = self.model.mesh
        #
        # attached_nodes = cap.node_ids
        #
        for boundary in mesh.boundaries:
            if cap.name.split("-")[0] in boundary.name:
                attached_nodes = boundary.node_ids
                break
        # -------------------------------------------------------------------

        # compute nodal areas:
        # 1. write vtk of volume, 2. read vtk, 3. extract surface, 4. compute nodal areas
        # NOTE: Should do this only once and not for every cap/valve involved
        filename = os.path.join(self.model.info.workdir, "temp_volume_mesh.vtk")
        mesh.write_to_vtk(filename)
        mesh_vtk = vtkmethods.read_vtk_unstructuredgrid_file(filename)
        os.remove(filename)

        surface_vtk = vtkmethods.vtk_surface_filter(mesh_vtk, True)
        nodal_areas = vtkmethods.compute_surface_nodal_area(surface_vtk)
        surface_obj = dsa.WrapDataObject(surface_vtk)
        surface_global_node_ids = surface_obj.PointData["GlobalPointIds"]

        # select only those nodal areas which match the cap node ids
        idx_select = np.nonzero(attached_nodes[:, None] == surface_global_node_ids)[1]
        nodal_areas = nodal_areas[idx_select]

        # scaled spring stiffness by nodal area
        scale_factor_normal *= nodal_areas
        scale_factor_radial *= nodal_areas

        # add part, section discrete, mat spring, sd_orientiation
        # element discrete

        # compute the radial components
        sd_orientations_radial = mesh.nodes[attached_nodes, :] - cap.centroid

        # normalize
        norms = np.linalg.norm(sd_orientations_radial, axis=1)
        sd_orientations_radial = sd_orientations_radial / norms[:, None]

        # add sd direction normal to plane
        vector_id_normal = self.id_offset["vector"]
        sd_orientation_normal_kw = create_define_sd_orientation_kw(
            vectors=cap.normal, vector_id_offset=vector_id_normal, iop=0
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

        # create discrete elements for normal direction
        nodes_discrete_elements = np.array(
            [attached_nodes + 1, np.zeros(len(attached_nodes))], dtype=int
        ).T
        vector_ids_normal = np.ones(len(attached_nodes), dtype=int) * vector_id_normal

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

        # discrete elements for radial direction
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

    def _add_pericardium_bc(self):
        """Add the pericardium.

        Note
        ----
        Uses the universal ventricular longitudinal coordinate
        and a sigmoid penalty function. Strocchi et al 2020 doi: 10.1016/j.jbiomech.2020.109645.
        """

        def _sigmoid(z):
            """Sigmoid function to scale spring coefficient."""
            return 1 / (1 + np.exp(-z))

        # compute penalty function
        uvc_l = self.model.mesh.point_data["uvc_longitudinal"]

        if np.any(uvc_l < 0):
            LOGGER.warning(
                "Negative normalized longitudinal coordinate detected."
                "Changing {0} negative uvc_l values to 1".format(np.sum((uvc_l < 0))),
            )

        uvc_l[uvc_l < 0] = 1
        penalty = (
            -_sigmoid(
                (abs(uvc_l) - self.parameters["Pericardium"]["Penalty function"][0])
                * self.parameters["Pericardium"]["Penalty function"][1]
            )
            + 1
        )

        # collect all pericardium nodes:
        epicardium_nodes = np.empty(0, dtype=int)
        epicardium_faces = np.empty((0, 3), dtype=int)
        LOGGER.debug("Collecting epicardium nodesets of ventricles:")
        ventricles = [part for part in self.model.parts if "ventricle" in part.name]
        epicardium_surfaces = [ventricle.epicardium for ventricle in ventricles]

        for surface in epicardium_surfaces:
            epicardium_nodes = np.append(epicardium_nodes, surface.node_ids)
            epicardium_faces = np.vstack([epicardium_faces, surface.faces])

        # NOTE: some duplicates may exist - fix this in preprocessor
        _, idx, counts = np.unique(epicardium_nodes, return_index=True, return_counts=True)
        if np.any(counts > 1):
            LOGGER.warning("Duplicate nodes found in pericardium")
        epicardium_nodes = epicardium_nodes[np.sort(idx)]

        # select only nodes that are on the epicardium and penalty factor > 0.1
        pericardium_nodes = epicardium_nodes[penalty[epicardium_nodes] > 0.001]

        spring_stiffness = self.parameters["Pericardium"]["Spring Stiffness"]  # kPA/mm
        # compute nodal areas:
        # NOTE: can be simplified
        filename = os.path.join(self.model.info.workdir, "temp_volume_mesh.vtk")
        self.model.mesh.write_to_vtk(filename)
        mesh_vtk = vtkmethods.read_vtk_unstructuredgrid_file(filename)
        os.remove(filename)
        vtk_surface = vtkmethods.vtk_surface_filter(mesh_vtk, True)
        nodal_areas = vtkmethods.compute_surface_nodal_area(vtk_surface)

        surface_obj = dsa.WrapDataObject(vtk_surface)
        surface_global_node_ids = surface_obj.PointData["GlobalPointIds"]

        # select only those nodal areas which match the pericardium node ids
        idx_select = np.nonzero(pericardium_nodes[:, None] == surface_global_node_ids)[1]
        nodal_areas = nodal_areas[idx_select]

        # compute scale factor
        scale_factors = nodal_areas * penalty[pericardium_nodes]

        # write to file
        # np.savetxt(
        #     os.path.join(self.model.info.workdir, "pericardium.txt"),
        #     np.concatenate(
        #         (
        #             self.model.mesh.nodes[pericardium_nodes, :],
        #             penalty[pericardium_nodes].reshape(-1, 1),
        #             scale_factors.reshape(-1, 1),
        #         ),
        #         axis=1,
        #     ),
        #     delimiter=",",
        # )

        # keywords
        # NOTE: Need to be made dynamic
        part_id = self.get_unique_part_id()
        section_id = self.get_unique_section_id()
        mat_id = self.get_unique_mat_id()

        part_kw = keywords.Part()
        part_kw.parts = pd.DataFrame(
            {"heading": ["Pericardium"], "pid": [part_id], "secid": [section_id], "mid": [mat_id]}
        )
        section_kw = keywords.SectionDiscrete(secid=section_id, cdl=0, tdl=0)
        mat_kw = keywords.MatSpringElastic(mid=mat_id, k=spring_stiffness)

        # 1: "omni-directional": equal springs in x,y, and z
        # 2: "apex-mitral-drection": one spring in apex-mitral valve direction
        spring_type = self.parameters["Pericardium"]["Spring Type"]

        if spring_type == "omni-directional":
            # create three unit vectors
            sd_orientation_kw = create_define_sd_orientation_kw(
                vectors=np.eye(3), vector_id_offset=self.id_offset["vector"]
            )

            self.id_offset["vector"] = sd_orientation_kw.vectors["vid"].to_numpy()[-1]
            # create discrete elements for each node and for each direction
            vector_ids = sd_orientation_kw.vectors["vid"].to_numpy()
            num_nodes = len(pericardium_nodes)
            vector_ids = np.tile(vector_ids, num_nodes)
            nodes = np.repeat(pericardium_nodes + 1, 3)
            nodes = np.vstack([nodes, np.zeros(len(nodes))])
            nodes = nodes.T
            scale_factors = np.repeat(scale_factors, 3)

        elif spring_type == "apex-mitral-direction":
            # get center of mitral valve
            left_ventricle = self.model.get_part("Left ventricle")
            apex_node_coordinates = left_ventricle.apex_points[0].xyz
            # midpoint between aortic and mitral valve
            center = np.mean([c.centroid for c in left_ventricle.caps], axis=0)

            # define spring orientation from apex to mitral valve
            orientation = center - apex_node_coordinates
            orientation /= np.linalg.norm(orientation)

            sd_orientation_kw = create_define_sd_orientation_kw(
                vectors=orientation, vector_id_offset=self.id_offset["vector"]
            )
            self.id_offset["vector"] = sd_orientation_kw.vectors["vid"].to_numpy()[-1]

            vector_ids = sd_orientation_kw.vectors["vid"].to_numpy()
            num_nodes = len(pericardium_nodes)
            vector_ids = np.tile(vector_ids, num_nodes)
            nodes = pericardium_nodes + 1
            nodes = np.vstack([nodes, np.zeros(len(nodes))])
            nodes = nodes.T

        # create discrete elements
        discrete_element_kw = create_discrete_elements_kw(
            nodes=nodes,
            part_id=part_id,
            vector_ids=vector_ids,
            scale_factor=scale_factors,
            element_id_offset=self.id_offset["element"]["discrete"],
        )
        # add offset
        self.id_offset["element"]["discrete"] = discrete_element_kw.elements["eid"].to_numpy()[-1]

        # add keywords to database
        self.kw_database.pericardium.append(part_kw)
        self.kw_database.pericardium.append(section_kw)
        self.kw_database.pericardium.append(mat_kw)
        self.kw_database.pericardium.append(sd_orientation_kw)
        self.kw_database.pericardium.append(discrete_element_kw)

        return

    def _update_cap_elements_db(self):
        """Update the database of shell elements.

        Note
        ----
        Loops over all the defined caps/valves.
        """
        # create part for each closing cap
        # used_partids = get_list_of_used_ids(self.kw_database.parts, "PART")
        # used_secids = get_list_of_used_ids(self.kw_database.parts, "SECTION")
        # used_segids = get_list_of_used_ids(self.kw_database.segment_sets, "SET_SEGMENT")

        section_id = self.get_unique_section_id()

        # NOTE should be dynamic
        mat_null_id = self.get_unique_mat_id()

        # material_kw = MaterialCap(mid=mat_null_id)

        material_kw = MaterialAtrium(
            mid=mat_null_id,
            rho=self.parameters["Cap"]["Density"],
            poisson_ratio=self.parameters["Cap"]["nu"],
            c10=self.parameters["Cap"]["c10"],
        )

        section_kw = keywords.SectionShell(
            secid=section_id,
            elform=4,
            shrf=0.8333,
            nip=3,
            t1=self.parameters["Cap"]["Thickness"],
        )

        self.kw_database.cap_elements.append(material_kw)
        self.kw_database.cap_elements.append(section_kw)

        caps = [cap for part in self.model.parts for cap in part.caps]
        # create new part for each cap
        cap_names_used = []
        for cap in caps:
            if cap.name in cap_names_used:
                LOGGER.debug("Already created material for {}: skipping".format(cap.name))
                continue

            cap.pid = self.get_unique_part_id()

            part_kw = keywords.Part()
            part_info = pd.DataFrame(
                {
                    "heading": [cap.name],
                    "pid": [cap.pid],
                    "secid": [section_id],
                    "mid": [mat_null_id],
                }
            )
            part_kw.parts = part_info

            self.kw_database.cap_elements.append(part_kw)
            cap_names_used.append(cap.name)

        # create closing triangles for each cap
        # assumes there are no shells written yet since offset = 0
        shell_id_offset = 0
        cap_names_used = []
        for cap in caps:
            if cap.name in cap_names_used:
                continue

            shell_kw = create_element_shell_keyword(
                shells=cap.triangles + 1,
                part_id=cap.pid,
                id_offset=shell_id_offset,
            )

            self.kw_database.cap_elements.append(shell_kw)

            shell_id_offset = shell_id_offset + cap.triangles.shape[0]
            cap_names_used.append(cap.name)
        return

    def _update_controlvolume_db(self):
        """Prepare the keywords for the control volume feature."""
        # NOTE: Assumes cavity id is reserved for combined
        # segment set

        # set up control volume keywords and interaction of
        # cavity with ambient. Only do for ventricles
        cavities = [part.cavity for part in self.model.parts if part.cavity]
        for cavity in cavities:
            if "atrium" in cavity.name:
                continue

            cv_kw = keywords.DefineControlVolume()
            cv_kw.id = cavity.surface.id
            cv_kw.sid = cavity.surface.id

            self.kw_database.control_volume.append(cv_kw)

        for cavity in cavities:
            if "atrium" in cavity.name:
                continue

            cvi_kw = keywords.DefineControlVolumeInteraction()
            cvi_kw.id = cavity.surface.id
            cvi_kw.cvid1 = cavity.surface.id
            cvi_kw.cvid2 = 0  # ambient

            # NOTE: static for the moment. Maximum of 2 cavities supported
            # but this is valid for the LeftVentricle, BiVentricle and FourChamber models
            if self.system_model_name == "ClosedLoop":
                if "Left ventricle" in cavity.name:
                    cvi_kw.lcid_ = -10
                elif "Right ventricle" in cavity.name:
                    cvi_kw.lcid_ = -11

            elif self.system_model_name == "ConstantPreloadWindkesselAfterload":
                if "Left ventricle" in cavity.name:
                    cvi_kw.lcid_ = 10
                if "Right ventricle" in cavity.name:
                    cvi_kw.lcid_ = 11

            self.kw_database.control_volume.append(cvi_kw)

        return

    def _update_system_model(self):
        """Update json system model settings."""
        model_type = self.model.info.model_type

        # closed loop uses a custom executable
        if self.system_model_name == "ClosedLoop":
            LOGGER.warning(
                "Note that this model type requires a custom executable that "
                "supports the Closed Loop circulation model!"
            )
            if isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
                file_path = os.path.join(
                    Path(__file__).parent.absolute(),
                    "templates",
                    "system_model_settings_bv.json",
                )

            elif isinstance(self.model, LeftVentricle):
                file_path = os.path.join(
                    Path(__file__).parent.absolute(),
                    "templates",
                    "system_model_settings_lv.json",
                )

            fid = open(file_path)
            sys_settings = json.load(fid)

            # update the volumes
            sys_settings["SystemModelInitialValues"]["UnstressedVolumes"][
                "lv"
            ] = self.model.get_part("Left ventricle").cavity.volume

            if isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
                sys_settings["SystemModelInitialValues"]["UnstressedVolumes"][
                    "rv"
                ] = self.model.get_part("Right ventricle").cavity.volume

            self.system_model_json = sys_settings

        # otherwise add the define function
        elif self.system_model_name == "ConstantPreloadWindkesselAfterload":
            from ansys.heart.writer.system_models import define_function_windkessel

            if self.system_model_name != self.parameters["Circulation System"]["Name"]:
                LOGGER.error("Circulation system parameters cannot be rad from Json")

            for cavity in self.model.cavities:
                if "Left ventricle" in cavity.name:

                    define_function_wk = define_function_windkessel(
                        function_id=10,
                        function_name="constant_preload_windkessel_afterload_left",
                        implicit=True,
                        constants=self.parameters["Circulation System"]["Left Ventricle"][
                            "Constant"
                        ],
                        initialvalues=self.parameters["Circulation System"]["Left Ventricle"][
                            "Initial Value"
                        ],
                    )
                    self.kw_database.control_volume.append(define_function_wk)

                elif "Right ventricle" in cavity.name:

                    define_function_wk = define_function_windkessel(
                        function_id=11,
                        function_name="constant_preload_windkessel_afterload_right",
                        implicit=True,
                        constants=self.parameters["Circulation System"]["Right Ventricle"][
                            "Constant"
                        ],
                        initialvalues=self.parameters["Circulation System"]["Right Ventricle"][
                            "Initial Value"
                        ],
                    )
                    self.kw_database.control_volume.append(define_function_wk)

        return

    def _add_enddiastolic_pressure_bc(self, pressure_lv: float = 1, pressure_rv: float = 1):
        """Add end diastolic pressure boundary condition on the left and right endocardium."""
        # create unit load curve
        load_curve_id = self.get_unique_curve_id()
        load_curve_kw = create_define_curve_kw(
            [0, 1, 1.001], [0, 1, 0], "unit load curve", load_curve_id, 100
        )

        if self.parameters["Unit"]["Time"] == "ms":
            load_curve_kw.sfa = 1000  # x scaling: to Millisecond

        # append unit curve to main.k
        self.kw_database.main.append(load_curve_kw)

        # create *LOAD_SEGMENT_SETS for each ventricular cavity
        cavities = [part.cavity for part in self.model.parts if part.cavity]
        for cavity in cavities:
            if "atrium" in cavity.name:
                continue

            if cavity.name == "Left ventricle":
                scale_factor = pressure_lv
                seg_id = cavity.surface.id
            elif cavity.name == "Right ventricle":
                scale_factor = pressure_rv
                seg_id = cavity.surface.id
            load_segset_kw = keywords.LoadSegmentSet(
                ssid=seg_id, lcid=load_curve_id, sf=scale_factor
            )
            self.kw_database.main.append(load_segset_kw)

        return


class ZeroPressureMechanicsDynaWriter(MechanicsDynaWriter):
    """
    Class for preparing the input for a stress-free LS-DYNA simulation.

    Note
    ----
    Derived from MechanicsDynaWriter and consequently derives all keywords relevant
    for simulations involving mechanics. This class does not write the
    control volume keywords but adds the keyword for computing the stress
    free configuration based on left/right cavity pressures instead.

    """

    def __init__(self, model: HeartModel) -> None:
        super().__init__(model)

        self.kw_database = MechanicsDecks()
        """Collection of keyword decks relevant for mechanics."""

        return

    def update(self):
        """Update the keyword database."""
        self._update_main_db(add_damping=False)

        self.kw_database.main.title = self.model.info.model_type + " zero-pressure"

        self._update_node_db()
        self._update_parts_db()
        self._update_solid_elements_db()
        self._update_segmentsets_db()
        self._update_nodesets_db()
        self._update_material_db(add_active=False)

        # for boundary conditions
        # self._update_boundary_conditions_db()
        self._add_cap_bc(bc_type="fix_caps")
        self._add_pericardium_bc()

        self._update_cap_elements_db()

        # # Approximate end-diastolic pressures
        pressure_lv = self.parameters["ED pressure"]["Left Ventricle"]  # kPa
        pressure_rv = self.parameters["ED pressure"]["Right Ventricle"]  # kPa

        self._add_enddiastolic_pressure_bc(pressure_lv=pressure_lv, pressure_rv=pressure_rv)

        # zerop key words
        self._add_control_reference_configuration()

        self._get_list_of_includes()
        self._add_includes()

        return

    def export(self, export_directory: str):
        """Write the model to files."""
        tstart = time.time()
        LOGGER.debug("Writing all LS-DYNA .k files...")

        if not export_directory:
            export_directory = os.path.join(self.model.info.workdir, "zeropressure")

        if not os.path.isdir(export_directory):
            os.makedirs(export_directory)

        # export .k files
        self.export_databases(export_directory)

        # export segment sets to separate file
        self._export_cavity_segmentsets(export_directory)

        tend = time.time()
        LOGGER.debug("Time spent writing files: {:.2f} s".format(tend - tstart))

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
        #
        # self.kw_database.main.append(keywords.DatabaseGlstat(dt=0.1, binary=2))
        #
        # self.kw_database.main.append(keywords.DatabaseMatsum(dt=0.1, binary=2))

        # frequency of full results
        self.kw_database.main.append(keywords.DatabaseBinaryD3Plot(dt=dt_output_d3plot))

        # self.kw_database.main.append(keywords.DatabaseExtentBinary(neiph=27, strflg=1, maxint=0))

        return

    def _add_solution_controls(self, scale_time=1):
        """Rewrite method for the zerop simulation."""
        self.kw_database.main.append(keywords.ControlTermination(endtim=1.0 * scale_time))

        self.kw_database.main.append(keywords.ControlImplicitDynamics(imass=0))
        # self.kw_database.main.append(
        #     keywords.ControlImplicitDynamics(imass=1, gamma=0.6, beta=0.38)
        # )

        # add auto controls
        self.kw_database.main.append(
            keywords.ControlImplicitAuto(iauto=1, dtmin=0.01 * scale_time, dtmax=0.1 * scale_time)
        )

        # add general implicit controls
        self.kw_database.main.append(
            keywords.ControlImplicitGeneral(imflag=1, dt0=0.1 * scale_time)
        )

        # add implicit solution controls: Defaults are OK?
        self.kw_database.main.append(keywords.ControlImplicitSolution())

        # add implicit solver controls
        self.kw_database.main.append(keywords.ControlImplicitSolver())

        # add binout for post-process
        self.kw_database.main.append(keywords.DatabaseNodout(dt=0.2 * scale_time, binary=1))

        # write for all nodes in nodout
        nodeset_id = self.get_unique_nodeset_id()
        kw = keywords.SetNodeGeneral(option="ALL", sid=nodeset_id)

        self.kw_database.main.append(kw)
        kw = keywords.DatabaseHistoryNodeSet(id1=nodeset_id)
        self.kw_database.main.append(kw)
        return

    def _add_control_reference_configuration(self):
        """Add control reference configuration keyword to main."""
        LOGGER.debug("Adding *CONTROL_REFERENCE_CONFIGURATION to main.k")
        kw = keywords.ControlReferenceConfiguration(maxiter=3, target="nodes.k", method=2, tol=5)

        self.kw_database.main.append(kw)

        return

    # def _add_enddiastolic_pressure_bc(self, pressure_lv: float = 1, pressure_rv: float = 1):
    #     """Adds end diastolic pressure boundary condition on the left and right endocardium"""

    #     # create unit load curve
    #     load_curve_id = 2
    #     load_curve_kw = create_define_curve_kw(
    #         [0, 1], [0, 1], "unit load curve", load_curve_id, 100
    #     )

    #     # append unit curve to main.k
    #     self.kw_database.main.append(load_curve_kw)

    #     # create *LOAD_SEGMENT_SETS for each ventricular cavity
    #     for cavity in self.model._mesh._cavities:

    #         if "atrium" in cavity.name:
    #             continue

    #         if cavity.name == "Left ventricle":
    #             scale_factor = pressure_lv
    #         elif cavity.name == "Right ventricle":
    #             scale_factor = pressure_rv

    #         LOGGER.debug(
    #             "Adding end-diastolic pressure of {0} to {1}".format(scale_factor, cavity.name)
    #         )

    #         seg_ids_to_use = []
    #         # find id of endocardium
    #         for segset in cavity.segment_sets:
    #             if "endocardium" in segset["name"]:
    #                 seg_ids_to_use.append(segset["id"])

    #         # create load segment set for each endocardium segment
    #         for seg_id in seg_ids_to_use:
    #             load_segset_kw = keywords.LoadSegmentSet(
    #                 ssid=seg_id, lcid=load_curve_id, sf=scale_factor
    #             )

    #             # append to main.k
    #             self.kw_database.main.append(load_segset_kw)


class FiberGenerationDynaWriter(MechanicsDynaWriter):
    """Class for preparing the input for a fiber-generation LS-DYNA simulation."""

    def __init__(self, model: HeartModel) -> None:
        super().__init__(model)

        self.kw_database = FiberGenerationDecks()
        """Collection of keywords relevant for fiber generation."""

    def update(self):
        """Update keyword database for Fiber generation: overwrites the inherited function."""
        ##
        self._update_main_db()  # needs updating

        self._update_node_db()  # can stay the same (could move to base class)
        if isinstance(self.model, (FourChamber, FullHeart)):
            self._keep_ventricles()

        self._update_parts_db()  # can stay the same (could move to base class++++++++++++++++++++)
        self._update_solid_elements_db(
            add_fibers=False
        )  # can stay the same (could move to base class)
        self._update_material_db()

        self._update_segmentsets_db()  # can stay the same
        self._update_nodesets_db()  # can stay the same

        # # update ep settings
        self._update_ep_settings()
        self._update_create_fibers()

        self._get_list_of_includes()
        self._add_includes()

        return

    def export(self, export_directory: str):
        """Write the model to files."""
        tstart = time.time()
        LOGGER.debug("Writing all LS-DYNA .k files...")

        if not export_directory:
            export_directory = self.model.info.workdir

        if not os.path.isdir(export_directory):
            os.makedirs(export_directory)

        # export .k files
        self.export_databases(export_directory)

        tend = time.time()
        LOGGER.debug("Time spent writing files: {:.2f} s".format(tend - tstart))

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

        for part in parts:
            element_ids = part.element_ids
            em_mat_id = self.get_unique_mat_id()
            self.kw_database.material.extend(
                [
                    keywords.MatElastic(mid=em_mat_id, ro=1e-6, e=1),
                    custom_keywords.EmMat003(
                        mid=em_mat_id,
                        mtype=2,
                        sigma11=5.0e-4,
                        sigma22=1.0e-4,
                        sigma33=1.0e-4,
                        beta=0.14,
                        cm=0.01,
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

    def _update_create_fibers(self):
        """Update the keywords for fiber generation."""
        # collect relevant node and segment sets.
        # node set: apex, base
        # node set: endocardium, epicardium
        # NOTE: could be better if basal nodes are extracted in the preprocessor
        # since that would allow you to robustly extract these nodessets using the
        # input data
        # The below is relevant for all models.
        nodes_base = np.empty(0, dtype=int)
        node_set_ids_endo = []  # relevant for both models
        node_sets_ids_epi = []  # relevant for both models
        node_set_ids_epi_and_rseptum = []  # only relevant for bv, 4c and full model

        # list of ventricular parts
        ventricles = [part for part in self.model.parts if "ventricle" in part.name]
        septum = self.model.get_part("Septum")

        # collect node set ids (already generated previously)
        node_set_ids_endo = [ventricle.endocardium.nsid for ventricle in ventricles]
        node_sets_ids_epi = [ventricle.epicardium.nsid for ventricle in ventricles]
        node_set_id_lv_endo = self.model.get_part("Left ventricle").endocardium.id

        if isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
            surfaces = [surface for p in self.model.parts for surface in p.surfaces]
            for surface in surfaces:
                if "septum" in surface.name and "endocardium" in surface.name:
                    node_set_ids_epi_and_rseptum = node_sets_ids_epi + [surface.id]
                    break

        for part in self.model.parts:
            for cap in part.caps:
                nodes_base = np.append(nodes_base, cap.node_ids)

        # apex id [0] endocardium, [1] epicardum
        apex_point = self.model.get_part("Left ventricle").apex_points[1]
        if "epicardium" not in apex_point.name:
            raise ValueError("Expecting a point on the epicardium")
        node_apex = apex_point.node_id

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
            LOGGER.warning("Model type %s in development " % self.model.info.model_type)

            # Define part set for myocardium
            part_list1_kw = keywords.SetPartList(
                sid=1,
            )
            part_list1_kw.parts._data = myocardium_part_ids
            part_list1_kw.options["TITLE"].active = True
            part_list1_kw.title = "myocardium_all"

            self.kw_database.create_fiber.extend(
                [
                    part_list1_kw,
                ]
            )

            # combine node sets endocardium uing *SET_NODE_ADD:
            node_set_id_all_endocardium = self.get_unique_nodeset_id()

            set_add_kw = keywords.SetNodeAdd(sid=node_set_id_all_endocardium)
            set_add_kw.options["TITLE"].active = True
            set_add_kw.title = "all_endocardium_segments"
            set_add_kw.nodes._data = node_set_ids_endo

            self.kw_database.create_fiber.append(set_add_kw)

            # combine node sets epicardium:
            node_set_id_all_epicardium = self.get_unique_nodeset_id()
            set_add_kw = keywords.SetNodeAdd(sid=node_set_id_all_epicardium)
            set_add_kw.options["TITLE"].active = True
            set_add_kw.title = "all_epicardium_segments"
            set_add_kw.nodes._data = node_sets_ids_epi

            self.kw_database.create_fiber.append(set_add_kw)

            node_set_id_base = self.get_unique_nodeset_id()
            node_set_id_apex = self.get_unique_nodeset_id()

            # create node-sets for base and apex
            node_set_base_kw = create_node_set_keyword(
                node_ids=nodes_base + 1, node_set_id=node_set_id_base, title="base nodes"
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
                    partsid=1, solvid1=1, solvid2=2, alpha=-101, beta=-102, wfile=1, prerun=1
                )
            )

            # define functions:
            from ansys.heart.writer.define_function_strings import function1, function2, function3

            self.kw_database.create_fiber.append(
                keywords.DefineFunction(fid=101, function=function1)
            )
            self.kw_database.create_fiber.append(
                keywords.DefineFunction(fid=102, function=function2)
            )

        elif isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
            LOGGER.warning("Model type %s under development " % self.model.info.model_type)

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
            set_add_kw.nodes._data = node_set_ids_endo

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
                node_ids=nodes_base + 1, node_set_id=node_set_id_base, title="base nodes"
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
                    partsid=1, solvid1=1, solvid2=2, alpha=-101, beta=-102, wfile=1, prerun=1
                )
            )
            # add *EM_EP_CREATEFIBERORIENTATION keywords
            self.kw_database.create_fiber.append(
                custom_keywords.EmEpCreatefiberorientation(
                    partsid=2, solvid1=1, solvid2=3, alpha=-101, beta=-103, wfile=1, prerun=1
                )
            )

            # define functions:
            from ansys.heart.writer.define_function_strings import function1, function2, function3

            self.kw_database.create_fiber.append(
                keywords.DefineFunction(fid=101, function=function1)
            )
            self.kw_database.create_fiber.append(
                keywords.DefineFunction(fid=102, function=function2)
            )
            self.kw_database.create_fiber.append(
                keywords.DefineFunction(fid=103, function=function3)
            )

    def _update_main_db(self):
        self.kw_database.main.append(keywords.ControlTimestep(dtinit=1.0, dt2ms=1.0))

        self.kw_database.main.append(keywords.ControlTermination(endtim=10))

        self.kw_database.main.append(keywords.DatabaseBinaryD3Plot(dt=1.0))

        return


class PurkinjeGenerationDynaWriter(MechanicsDynaWriter):
    """Class for preparing the input for a Purkinje LS-DYNA simulation."""

    def __init__(self, model: HeartModel) -> None:
        super().__init__(model)

        self.kw_database = PurkinjeGenerationDecks()
        """Collection of keywords relevant for Purkinje generation."""

    def update(self):
        """Update keyword database - overwrites the inherited function."""
        ##
        self._update_main_db()  # needs updating

        self._update_node_db()  # can stay the same (could move to base class)
        if isinstance(self.model, (FourChamber, FullHeart)):
            self._keep_ventricles()

        self._update_parts_db()  # can stay the same (could move to base class++++++++++++++++++++)
        self._update_solid_elements_db(
            add_fibers=False
        )  # can stay the same (could move to base class)
        self._update_material_db()

        self._update_segmentsets_db()  # can stay the same
        self._update_nodesets_db()  # can stay the same

        # update ep settings
        self._update_ep_settings()
        self._update_create_Purkinje()

        self._get_list_of_includes()
        self._add_includes()

        return

    def export(self, export_directory: str):
        """Write the model to files."""
        tstart = time.time()
        LOGGER.debug("Writing all LS-DYNA .k files...")

        if not export_directory:
            export_directory = self.model.info.workdir

        if not os.path.isdir(export_directory):
            os.makedirs(export_directory)

        # export .k files
        self.export_databases(export_directory)

        tend = time.time()
        LOGGER.debug("Time spent writing files: {:.2f} s".format(tend - tstart))

        return

    def _update_material_db(self):
        """Add simple linear elastic material for each defined part."""
        for part in self.model.parts:
            part.element_ids
            em_mat_id = part.pid
            self.kw_database.material.extend(
                [
                    keywords.MatElastic(mid=em_mat_id, ro=1e-6, e=1),
                    custom_keywords.EmMat003(
                        mid=em_mat_id,
                        mtype=2,
                        sigma11=5.0e-4,
                        sigma22=1.0e-4,
                        sigma33=1.0e-4,
                        beta=0.14,
                        cm=0.01,
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

    def _update_create_Purkinje(self):
        """Update the keywords for Purkinje generation."""
        # collect relevant node and segment sets.
        # node set: apex, base
        # node set: endocardium, epicardium
        # NOTE: could be better if basal nodes are extracted in the preprocessor
        # since that would allow you to robustly extract these nodessets using the
        # input data
        # The below is relevant for all models.
        node_apex_left = np.empty(0, dtype=int)
        node_apex_right = np.empty(0, dtype=int)
        edge_id_start_left = np.empty(0, dtype=int)
        edge_id_start_right = np.empty(0, dtype=int)

        # apex_points[0]: endocardium, apex_points[1]: epicardium
        if isinstance(self.model, (LeftVentricle, BiVentricle, FourChamber, FullHeart)):
            node_apex_left = self.model.left_ventricle.apex_points[0].node_id
            segment_set_ids_endo_left = self.model.left_ventricle.endocardium.id

            # check whether point is on edge of endocardium - otherwise pick another node in
            # the same triangle
            endocardium = self.model.left_ventricle.endocardium
            endocardium.get_boundary_edges()
            if np.any(endocardium.boundary_edges == node_apex_left):
                element_id = np.argwhere(np.any(endocardium.faces == node_apex_left, axis=1))[0][0]

                node_apex_left = endocardium.faces[element_id, :][
                    np.argwhere(
                        np.isin(
                            endocardium.faces[element_id, :],
                            endocardium.boundary_edges,
                            invert=True,
                        )
                    )[0][0]
                ]
                LOGGER.warning(
                    "Node id {0} is on edge of {1}. Picking node id {2}".format(
                        self.model.left_ventricle.apex_points[0].node_id,
                        endocardium.name,
                        node_apex_left,
                    )
                )
                self.model.left_ventricle.apex_points[0].node_id = node_apex_left

        if isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
            node_apex_right = self.model.right_ventricle.apex_points[0].node_id
            segment_set_ids_endo_right = self.model.right_ventricle.endocardium.id

            # check whether point is on edge of endocardium - otherwise pick another node in
            # the same triangle
            endocardium = self.model.right_ventricle.endocardium
            endocardium.get_boundary_edges()
            if np.any(endocardium.boundary_edges == node_apex_right):
                element_id = np.argwhere(np.any(endocardium.faces == node_apex_right, axis=1))[0][0]

                node_apex_right = endocardium.faces[element_id, :][
                    np.argwhere(
                        np.isin(
                            endocardium.faces[element_id, :],
                            endocardium.boundary_edges,
                            invert=True,
                        )
                    )[0][0]
                ]
                LOGGER.warning(
                    "Node id {0} is on edge of {1}. Picking node id {2}".format(
                        self.model.right_ventricle.apex_points[0].node_id,
                        endocardium.name,
                        node_apex_right,
                    )
                )
                self.model.right_ventricle.apex_points[0].node_id = node_apex_right

        # check whether apical points are on the endge of the endocardium

        # NOTE: is this still relevant applicable with the new structure?
        # NOTE: validate node set by removing any nodes that do not occur in either ventricle
        # tet_ids_ventricles = np.empty((0), dtype=int)
        # for cavity in self.model._mesh._cavities:
        #     for element_set in cavity.element_sets:
        #         if "ventricle" in cavity.name:
        #             tet_ids_ventricles = np.append(tet_ids_ventricles, element_set["set"])
        # tetra_ventricles = self.volume_mesh["tetra"][tet_ids_ventricles, :]

        # # remove nodes that occur just in atrial part
        # mask = np.isin(nodes_base, tetra_ventricles, invert=True)
        # LOGGER.debug("Removing {0} nodes from base nodes".format(np.sum(mask)))
        # nodes_base = nodes_base[np.invert(mask)]

        node_set_id_apex_left = self.get_unique_nodeset_id()
        # create node-sets for apex
        node_set_apex_kw = create_node_set_keyword(
            node_ids=[node_apex_left + 1],
            node_set_id=node_set_id_apex_left,
            title="apex node left",
        )

        self.kw_database.node_sets.append(node_set_apex_kw)

        apex_left_coordinates = self.model.mesh.nodes[node_apex_left, :]

        node_id_start_left = (
            self.model.mesh.nodes.shape[0] + 1
        )  # TODO seek for max id rather than number of rows

        edge_id_start_left = self.model.mesh.tetrahedrons.shape[0] + 1

        # Purkinje generation parameters
        self.kw_database.main_left_ventricle.append(
            custom_keywords.EmEpPurkinjeNetwork(
                purkid=1,
                buildnet=1,
                ssid=segment_set_ids_endo_left,
                mid=25,
                pointstx=apex_left_coordinates[0],
                pointsty=apex_left_coordinates[1],
                pointstz=apex_left_coordinates[2],
                edgelen=2,
                ngen=50,
                nbrinit=8,
                nsplit=2,
                inodeid=node_id_start_left,
                iedgeid=edge_id_start_left,  # TODO check if beam elements exist in mesh
            )
        )

        # Add right purkinje only in biventricular or 4chamber models
        if isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
            LOGGER.warning("Model type %s in development " % self.model.info.model_type)

            node_set_id_apex_right = self.get_unique_nodeset_id()
            # create node-sets for apex
            node_set_apex_kw = create_node_set_keyword(
                node_ids=[node_apex_right + 1],
                node_set_id=node_set_id_apex_right,
                title="apex node right",
            )

            self.kw_database.node_sets.append(node_set_apex_kw)

            apex_right_coordinates = self.model.mesh.nodes[node_apex_right, :]

            node_id_start_right = (
                2 * self.model.mesh.nodes.shape[0]
            )  # TODO find a solution in dyna to better handle id definition

            edge_id_start_right = 2 * self.model.mesh.tetrahedrons.shape[0]

            # Purkinje generation parameters
            self.kw_database.main_right_ventricle.append(
                custom_keywords.EmEpPurkinjeNetwork(
                    purkid=2,
                    buildnet=1,
                    ssid=segment_set_ids_endo_right,
                    mid=26,
                    pointstx=apex_right_coordinates[0],
                    pointsty=apex_right_coordinates[1],
                    pointstz=apex_right_coordinates[2],
                    edgelen=2,
                    ngen=50,
                    nbrinit=8,
                    nsplit=2,
                    inodeid=node_id_start_right,  # TODO check if beam elements exist in mesh
                    iedgeid=edge_id_start_right,
                )
            )

    def _update_main_db(self):

        return

    def _get_list_of_includes(self):
        """Get a list of files to include in main.k. Omit any empty decks."""
        for deckname, deck in vars(self.kw_database).items():
            if deckname == "main_left_ventricle" or deckname == "main_right_ventricle":
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
            self.kw_database.main_left_ventricle.append(
                keywords.Include(filename=filename_to_include)
            )
            self.kw_database.main_right_ventricle.append(
                keywords.Include(filename=filename_to_include)
            )


class ElectrophysiologyDynaWriter(BaseDynaWriter):
    """Class for preparing the input for a Electrophysiology LS-DYNA simulation."""

    def __init__(self, model: HeartModel) -> None:
        super().__init__(model)

        self.kw_database = ElectrophysiologyDecks()
        """Collection of keywords relevant for Electrophysiology."""

    def update(self):
        """Update keyword database for Electrophysiology."""
        ##
        self._update_main_db()

        self._update_solution_controls()
        self._update_export_controls()
        self._update_node_db()

        self._update_parts_db()
        self._update_solid_elements_db()
        self._update_material_db()
        self._update_cellmodels()
        self._update_segmentsets_db()
        self._update_nodesets_db()

        # update ep settings
        self._update_ep_settings()

        self._get_list_of_includes()
        self._add_includes()

        return

    def export(self, export_directory: str):
        """Write the model to files."""
        tstart = time.time()
        LOGGER.debug("Writing all LS-DYNA .k files...")

        if not export_directory:
            export_directory = self.model.info.workdir

        if not os.path.isdir(export_directory):
            os.makedirs(export_directory)

        # export .k files
        self.export_databases(export_directory)

        tend = time.time()
        LOGGER.debug("Time spent writing files: {:.2f} s".format(tend - tstart))

        return

    def _update_material_db(self):
        """Add simple linear elastic material for each defined part."""
        for part in self.model.parts:
            ep_mid = part.pid
            self.kw_database.material.extend(
                [
                    keywords.MatElastic(mid=ep_mid, ro=1e-6, e=1),
                    custom_keywords.EmMat003(
                        mid=ep_mid,
                        mtype=2,
                        sigma11=5.0e-4,
                        sigma22=1.0e-4,
                        sigma33=1.0e-4,
                        beta=0.14,
                        cm=0.01,
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

    def _update_cellmodels(self):
        """Add simple linear elastic material for each defined part."""
        for part in self.model.parts:
            ep_mid = part.pid
            self.kw_database.cell_models.extend(
                [
                    keywords.EmEpCellmodelTentusscher(
                        mid=ep_mid,
                        gas_constant=8314.472,
                        t=310,
                        faraday_constant=96485.3415,
                        cm=0.185,
                        vc=0.016404,
                        vsr=0.001094,
                        vss=0.00005468,
                        pkna=0.03,
                        ko=5.4,
                        nao=140.0,
                        cao=2.0,
                        gk1=5.405,
                        gkr=0.153,
                        gks=0.392,
                        gna=14.838,
                        gbna=0.0002,
                        gcal=0.0000398,
                        gbca=0.000592,
                        gto=0.294,
                        gpca=0.1238,
                        gpk=0.0146,
                        pnak=2.724,
                        km=1.0,
                        kmna=40.0,
                        knaca=1000.0,
                        ksat=0.1,
                        alpha=2.5,
                        gamma=0.35,
                        kmca=1.38,
                        kmnai=87.5,
                        kpca=0.0005,
                        k1=0.15,
                        k2=0.045,
                        k3=0.06,
                        k4=0.005,
                        ec=1.5,
                        maxsr=2.5,
                        minsr=1.0,
                        vrel=0.102,
                        vleak=0.00036,
                        vxfer=0.0038,
                        vmaxup=0.006375,
                        kup=0.00025,
                        bufc=0.2,
                        kbufc=0.001,
                        bufsr=10.0,
                        kbufsf=0.3,
                        bufss=0.4,
                        kbufss=0.00025,
                        v=-85.23,
                        ki=136.89,
                        nai=8.604,
                        cai=0.000126,
                        cass=0.00036,
                        casr=3.64,
                        rpri=0.9073,
                        xr1=0.00621,
                        xr2=0.4712,
                        xs=0.0095,
                        m=0.00172,
                        h=0.7444,
                        j=0.7045,
                        d=3.373e-5,
                        f=0.7888,
                        f2=0.9755,
                        fcass=0.9953,
                        s=0.999998,
                        r=2.42e-8,
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

        # use defaults
        self.kw_database.ep_settings.append(custom_keywords.EmControlEp(numsplit=5))

        # max iter should be int
        self.kw_database.ep_settings.append(
            keywords.EmSolverFem(reltol=1e-6, maxite=int(1e4), precon=2)
        )

        self.kw_database.ep_settings.append(keywords.EmOutput(mats=1, matf=1, sols=1, solf=1))

        node_apex_left = self.get_apex_left()
        node_apex_right = self.get_apex_right()

        node_set_id_apex_left = self.get_unique_nodeset_id()
        # create node-sets for apex left
        node_set_kw = create_node_set_keyword(
            node_ids=[node_apex_left + 1],
            node_set_id=node_set_id_apex_left,
            title="apex node left",
        )
        self.kw_database.node_sets.append(node_set_kw)

        node_set_id_apex_right = self.get_unique_nodeset_id()
        # create node-sets for apex right
        node_set_kw = create_node_set_keyword(
            node_ids=[node_apex_right + 1],
            node_set_id=node_set_id_apex_right,
            title="apex node right",
        )
        self.kw_database.node_sets.append(node_set_kw)
        # TODO add more nodes to initiate wave propagation !!!!
        node_set_id_stimulationnodes = self.get_unique_nodeset_id()
        # create node-sets for apex
        node_set_kw = create_node_set_keyword(
            node_ids=[node_apex_left + 1, node_apex_right + 1],
            node_set_id=node_set_id_stimulationnodes,
            title="Stim nodes",
        )
        self.kw_database.node_sets.append(node_set_kw)

        self.kw_database.ep_settings.append(
            custom_keywords.EmEpTentusscherStimulus(
                stimid=1,
                settype=2,
                setid=node_set_id_stimulationnodes,
                stimstrt=0.0,
                stimt=1000.0,
                stimdur=20.0,
                stimamp=50.0,
            )
        )
        return

    # def _update_use_Purkinje(self):
    #     """Updates the keywords for Purkinje generation."""

    #     # collect relevant node and segment sets.
    #     # node set: apex, base
    #     # node set: endocardium, epicardium
    #     # NOTE: could be better if basal nodes are extracted in the preprocessor
    #     # since that would allow you to robustly extract these nodessets using the
    #     # input data
    #     # The below is relevant for all models.
    #     nodes_base = np.empty(0, dtype=int)
    #     node_apex_left = np.empty(0, dtype=int)
    #     node_apex_right = np.empty(0, dtype=int)
    #     edge_id_start_left = np.empty(0, dtype=int)
    #     edge_id_start_right = np.empty(0, dtype=int)
    #     for cavity in self.model._mesh._cavities:
    #         if cavity.name == "Left ventricle":
    #             node_apex_left = cavity.apex_id["endocardium"]
    #             for segment_set in cavity.segment_sets:
    #                 if "endocardium" in segment_set["name"]:
    #                     segment_set_ids_endo_left = segment_set["id"]
    #         elif cavity.name == "Right ventricle":
    #             node_apex_right = cavity.apex_id["endocardium"]
    #             for segment_set in cavity.segment_sets:
    #                 if "endocardium" in segment_set["name"] and "septum" not
    #                       in segment_set["name"]:
    #                     segment_set_ids_endo_right = segment_set["id"]

    #     # validate node set by removing any nodes that do not occur in either ventricle
    #     # NOTE: can be much more consice
    #     tet_ids_ventricles = np.empty((0), dtype=int)
    #     for cavity in self.model._mesh._cavities:
    #         for element_set in cavity.element_sets:
    #             if "ventricle" in cavity.name:
    #                 tet_ids_ventricles = np.append(tet_ids_ventricles, element_set["set"])
    #     tetra_ventricles = self.volume_mesh["tetra"][tet_ids_ventricles, :]

    #     # remove nodes that occur just in atrial part
    #     mask = np.isin(nodes_base, tetra_ventricles, invert=True)
    #     LOGGER.debug("Removing {0} nodes from base nodes".format(np.sum(mask)))
    #     nodes_base = nodes_base[np.invert(mask)]

    #     node_set_id_apex_left = 201
    #     # create node-sets for apex
    #     node_set_apex_kw = create_node_set_keyword(
    #         node_ids=[node_apex_left + 1],
    #         node_set_id=node_set_id_apex_left,
    #         title="apex node left",
    #     )

    #     self.kw_database.node_sets.extend([node_set_apex_kw])

    #     apex_left_X = self.volume_mesh["nodes"][node_apex_left, 0]
    #     apex_left_Y = self.volume_mesh["nodes"][node_apex_left, 1]
    #     apex_left_Z = self.volume_mesh["nodes"][node_apex_left, 2]
    #     node_id_start_left = (
    #         self.volume_mesh["nodes"].shape[0] + 1
    #     )  # TODO seek for max id rather than number of rows

    #     edge_id_start_left = self.volume_mesh["tetra"].shape[0] + 1

    #     # Purkinje generation parameters
    #     self.kw_database.mainLEFT.append(
    #         custom_keywords.EmEpPurkinjeNetwork(
    #             purkid=1,
    #             buildnet=0,
    #             ssid=segment_set_ids_endo_left,
    #             mid=25,
    #             pointstx=apex_left_X,
    #             pointsty=apex_left_Y,
    #             pointstz=apex_left_Z,
    #             edgelen=2,
    #             ngen=50,
    #             nbrinit=8,
    #             nsplit=2,
    #             inodeid=node_id_start_left,
    #             iedgeid=edge_id_start_left,  # TODO check if beam elements exist in mesh
    #         )
    #     )

    #     # Add right purkinje only in biventricular or 4chamber models
    #     if self.model.info.model_type in ["BiVentricle", "FourChamber"]:
    #         LOGGER.warning("Model type %s in development " % self.model.info.model_type)

    #         node_set_id_apex_right = 202
    #         # create node-sets for apex
    #         node_set_apex_kw = create_node_set_keyword(
    #             node_ids=[node_apex_right + 1],
    #             node_set_id=node_set_id_apex_right,
    #             title="apex node right",
    #         )

    #         self.kw_database.node_sets.extend([node_set_apex_kw])

    #         apex_right_X = self.volume_mesh["nodes"][node_apex_right, 0]
    #         apex_right_Y = self.volume_mesh["nodes"][node_apex_right, 1]
    #         apex_right_Z = self.volume_mesh["nodes"][node_apex_right, 2]
    #         node_id_start_right = (
    #             2 * self.volume_mesh["nodes"].shape[0]
    #         )  # TODO find a solution in dyna to better handle id definition

    #         edge_id_start_right = 2 * self.volume_mesh["tetra"].shape[0]

    #         # Purkinje generation parameters
    #         self.kw_database.mainRIGHT.append(
    #             custom_keywords.EmEpPurkinjeNetwork(
    #                 purkid=2,
    #                 buildnet=0,
    #                 ssid=segment_set_ids_endo_right,
    #                 mid=26,
    #                 pointstx=apex_right_X,
    #                 pointsty=apex_right_Y,
    #                 pointstz=apex_right_Z,
    #                 edgelen=2,
    #                 ngen=50,
    #                 nbrinit=8,
    #                 nsplit=2,
    #                 inodeid=node_id_start_right,  # TODO check if beam elements exist in mesh
    #                 iedgeid=edge_id_start_right,
    #             )
    #         )

    def _update_solution_controls(
        self,
        end_time: float = 800,
    ):
        """Add solution controls, output controls and solver settings."""
        # add termination keywords
        self.kw_database.main.append(keywords.ControlTermination(endtim=end_time, dtmin=0.0))

        self.kw_database.main.append(keywords.ControlTimeStep(dtinit=1.0, dt2ms=1.0))
        return

    def _update_export_controls(self, dt_output_d3plot: float = 1.0):
        """Add solution controls to the main simulation.

        Parameters
        ----------
        dt_output_d3plot : float, optional
            Writes full D3PLOT results at this time-step spacing, by default 0.05
        dt_output_icvout : float, optional
            Writes control volume results at this time-step spacing, by default 0.001

        """
        # frequency of full results
        self.kw_database.main.append(keywords.DatabaseBinaryD3Plot(dt=dt_output_d3plot))

        return

    def _update_main_db(self):

        return

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


if __name__ == "__main__":
    print("protected")
