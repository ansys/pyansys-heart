"""Module contain. classes for writing LS-DYNA keywords based.

Note
----
Uses a HeartModel (from ansys.heart.preprocessor.models).

"""

import copy
import json

# import missing keywords
import logging
import os
import time
from typing import List, Literal

from ansys.dyna.keywords import keywords

LOGGER = logging.getLogger("pyheart_global.writer")
# from importlib.resources import files
from importlib.resources import path as resource_path

from ansys.heart.preprocessor.material.material import MAT295, MechanicalMaterialModel, NeoHookean
from ansys.heart.preprocessor.mesh.objects import Cap
import ansys.heart.preprocessor.mesh.vtkmethods as vtkmethods

global heart_version
heart_version = os.getenv("ANSYS_HEART_MODEL_VERSION")
if not heart_version:
    heart_version = "v0.1"

if heart_version == "v0.2":
    from ansys.heart.preprocessor.models.v0_2.models import (
        BiVentricle,
        FourChamber,
        FullHeart,
        HeartModel,
        LeftVentricle,
    )
elif heart_version == "v0.1":
    from ansys.heart.preprocessor.models.v0_1.models import (
        BiVentricle,
        FourChamber,
        FullHeart,
        HeartModel,
        LeftVentricle,
    )

from ansys.heart.simulator.settings.settings import SimulationSettings
from ansys.heart.writer import custom_dynalib_keywords as custom_keywords
from ansys.heart.writer.heart_decks import (
    BaseDecks,
    ElectroMechanicsDecks,
    ElectrophysiologyDecks,
    FiberGenerationDecks,
    MechanicsDecks,
    PurkinjeGenerationDecks,
)
from ansys.heart.writer.keyword_module import (
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
from ansys.heart.writer.material_keywords import (
    MaterialHGOMyocardium2,
    MaterialNeoHook,
    _Depracated_MaterialHGOMyocardium,
    active_curve,
)
from ansys.heart.writer.system_models import _ed_load_template, define_function_windkessel
import numpy as np
import pandas as pd
import pyvista as pv


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
        <Example to be added>
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

        id = self.id_offset["part"]
        for part in self.model.parts:
            id += 1
            # cannot use get_unique_part_id() because it checks in Deck()
            # part.pid = self.get_unique_part_id()
            part.pid = id
        """Assign part id for heart parts."""

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

        # assign material for part if it's empty
        myocardium, neohookean = self.settings.get_mechanical_material()
        for part in self.model.parts:
            if isinstance(part.meca_material, MechanicalMaterialModel.DummyMaterial):
                LOGGER.info(f"Material of {part.name} will be assigned automatically.")
                if part.has_fiber:
                    if part.is_active:
                        part.meca_material = myocardium
                    else:
                        part.meca_material = copy.deepcopy(myocardium)
                        part.meca_material.active = None
                elif not part.has_fiber:
                    part.meca_material = neohookean

        if "Improved" in self.model.info.model_type:
            LOGGER.warning(
                "Changing model type from : {0} to {1}".format(
                    self.model.info.model_type, self.model.info.model_type.replace("Improved", "")
                )
            )
            self.model.info.model_type = self.model.info.model_type.replace("Improved", "")

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
            if not sett.Mechanics in subsettings_classes:
                raise ValueError("Expecting mechanics settings.")

        elif isinstance(self, FiberGenerationDynaWriter):
            if not sett.Fibers in subsettings_classes:
                raise ValueError("Expecting fiber settings.")

        elif isinstance(self, PurkinjeGenerationDynaWriter):
            if not sett.Purkinje in subsettings_classes:
                raise ValueError("Expecting Purkinje settings.")

        elif isinstance(self, ElectrophysiologyDynaWriter):
            if not sett.Electrophysiology in subsettings_classes:
                raise ValueError("Expecting electrophysiology settings.")
        elif isinstance(self, UHCWriter):
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
            nodes = np.vstack([ids + 1, self.model.mesh.nodes[ids, :].T]).T
            node_kw = add_nodes_to_kw(nodes, node_kw)
        else:
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
                segments=cavity.surface.triangles + 1,
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
                    segments=surface.triangles + 1,
                    segid=surface.id,
                    title=surface.name,
                )
                # append this kw to the segment set database
                self.kw_database.segment_sets.append(kw)

        return

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
        Note
        ----
            In FiberGenerationWriter, we do not allow all nodes of same element in one nodeset.
        """
        # formats endo, epi- and septum nodeset keywords, do for all surfaces

        surface_ids = [s.id for p in self.model.parts for s in p.surfaces]
        node_set_id = np.max(surface_ids) + 1

        # for each surface in each part add the respective node-set
        # Use same ID as surface
        # TODO check if database already contains nodesets (there will be duplicates otherwise)
        used_node_ids = np.empty(0, dtype=int)

        for part in self.model.parts:
            kws_surface = []
            for surface in part.surfaces:
                if remove_duplicates:
                    node_ids = np.setdiff1d(surface.node_ids, used_node_ids)
                else:
                    node_ids = surface.node_ids
                if remove_one_node_from_cell:
                    # make sure not all nodes of the same elements are in the surface
                    node_mask = np.zeros(self.model.mesh.number_of_points, dtype=int)
                    # tag surface nodes with value 1
                    node_mask[node_ids] = 1

                    cell_mask = np.array(
                        [
                            node_mask[self.model.mesh.tetrahedrons[:, 0]],
                            node_mask[self.model.mesh.tetrahedrons[:, 1]],
                            node_mask[self.model.mesh.tetrahedrons[:, 2]],
                            node_mask[self.model.mesh.tetrahedrons[:, 3]],
                        ]
                    )
                    # cells with all nodes in surface are those whose
                    # all nodes are tagged with value 1
                    issue_cells = np.where(np.sum(cell_mask, axis=0) == 4)[0]
                    nodes_toremove = self.model.mesh.tetrahedrons[issue_cells, :][:, 0]
                    node_ids = np.setdiff1d(node_ids, nodes_toremove)

                    for cell in issue_cells:
                        LOGGER.warning(
                            f"All nodes of cell {cell+1} are in nodeset of {surface.name},"
                            " N1 is removed."
                        )

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
        LOGGER.debug("Writing all LS-DYNA .k files...")

        # is this reachable??
        if not export_directory:
            export_directory = os.path.join(
                self.model.info.workdir, self.__class__.__name__.lower().replace("dynawriter", "")
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

    def _export_cavity_segmentsets(self, export_directory: str):
        """Export the cavity segment sets to separate files."""
        cavities = [part.cavity for part in self.model.parts if part.cavity]
        for cavity in cavities:
            filepath = os.path.join(
                export_directory, "_".join(cavity.name.lower().split()) + ".segment"
            )
            np.savetxt(filepath, cavity.surface.triangles + 1, delimiter=",", fmt="%d")

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

    def _keep_atria(self):
        """Remove any non-atrial parts."""
        # NOTE: Could move "remove part" method to model
        LOGGER.warning("Just keeping atrial-parts")
        parts_to_keep = [part for part in self.model.part_names if "atrium" in part.lower()]
        parts_to_remove = [part for part in self.model.part_names if part not in parts_to_keep]
        for part_to_remove in parts_to_remove:
            LOGGER.warning("Removing: {}".format(part_to_remove))
            self.model.remove_part(part_to_remove)
        return

    def _keep_parts(self, parts_to_keep: List[str]):
        """Remove parts by a list of part names."""
        parts_to_remove = [part for part in self.model.part_names if part not in parts_to_keep]
        for part_to_remove in parts_to_remove:
            LOGGER.warning("Removing: {}".format(part_to_remove))
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
            if add_fibers and part.has_fiber:
                part_add_fibers = True
            else:
                part_add_fibers = False

            LOGGER.debug(
                "\tAdding elements for {0} | adding fibers: {1}".format(part.name, part_add_fibers)
            )
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

        # Depending on the system model specified give list of parameters
        self.cap_in_zerop = True
        """
        If include cap (shell) elements in ZeroPressure.
        Experimental feature, please do not change it.
        """
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

    def update(self, with_dynain=False):
        """
        Update the keyword database.

        Parameters
        ----------
        with_dynain: bool, optional
            Use dynain.lsda file from stress free configuration computation.
        """
        self._update_main_db()

        self._add_damping()

        self._update_parts_db()

        if not with_dynain:
            self._update_node_db()
            self._update_solid_elements_db(add_fibers=True)
            # no zerop exists, cap mesh need to be written
            self.cap_in_zerop = False
        else:
            self.kw_database.main.append(keywords.Include(filename="dynain.lsda"))

        self._update_segmentsets_db()
        self._update_nodesets_db()
        self._update_material_db(add_active=True)

        if self.cap_in_zerop:
            # mesh has been defined in Zerop so saved in dynain file
            self._update_cap_elements_db(add_mesh=False)
        else:
            # define new cap element
            self._update_cap_elements_db(add_mesh=True)

        # # for control volume
        self._update_controlvolume_db()
        self._update_system_model()

        # no control volume for atrial, constant pressure instead
        if isinstance(self.model, FourChamber):
            bc_settings = self.settings.mechanics.boundary_conditions
            pressure_lv = bc_settings.end_diastolic_cavity_pressure["left_ventricle"].m
            pressure_rv = bc_settings.end_diastolic_cavity_pressure["right_ventricle"].m
            self._add_constant_atrial_pressure(pressure_lv=pressure_lv, pressure_rv=pressure_rv)

        # for boundary conditions
        self._add_cap_bc(bc_type="springs_caps")
        self._add_pericardium_bc()

        self._get_list_of_includes()
        self._add_includes()

        return

    def export(self, export_directory: str):
        """Write the model to files."""
        super().export(export_directory)

        # write cavity name and volume
        dct = {}
        for cavity in self.model.cavities:
            dct[cavity.name] = cavity.volume
        with open(os.path.join(export_directory, "volumes.json"), "w") as f:
            json.dump(dct, f)

        # todo: Close loop is only available from a customized LSDYNA executable
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
        self.kw_database.main.title = self.model.model_type

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

    # def _update_node_db(self):
    #     """Add nodes to the NODE database."""
    #     LOGGER.debug("Updating node keywords...")
    #     node_kw = keywords.Node()
    #     node_kw = add_nodes_to_kw(self.model.mesh.nodes, node_kw)
    #
    #     self.kw_database.nodes.append(node_kw)
    #
    #     return

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
                "Simulation type not recognized: Please choose " "either quasi-static or static"
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
                maxref=35,
                dctol=0.02,
                ectol=1e6,
                rctol=1e3,
                abstol=-1e-20,
                dnorm=1,
                diverg=2,
                lstol=-0.9,
                lsmtd=5,
                d3itctl=1,
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

        # remove, aha strain is computed from d3plot

        # # control ELOUT file to extract left ventricle's stress/strain
        # if hasattr(self.model, "septum"):
        #     self.kw_database.main.append(
        #         keywords.SetSolidGeneral(
        #             option="PART",
        #             sid=1,
        #             e1=self.model.left_ventricle.pid,
        #             e2=self.model.septum.pid,
        #             user_comment="create left ventricle + septum set for exporting",
        #         )
        #     )
        # else:
        #     self.kw_database.main.append(
        #         keywords.SetSolidGeneral(option="PART", sid=1, e1=self.model.left_ventricle.pid)
        #     )
        # self.kw_database.main.append(keywords.DatabaseHistorySolidSet(id1=1))

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
            kw = keywords.DampingPartStiffness(pid=part.pid, coef=-0.2)
            self.kw_database.main.append(kw)
        return

    def _update_segmentsets_db(self):
        """Update the segment set database."""
        # NOTE 0: add all surfaces as segment sets
        # NOTE 1: need to more robustly check segids that are already used?

        # add closed cavity segment sets
        cavities = [p.cavity for p in self.model.parts if p.cavity]

        # caps = [cap for part in self.model.parts for cap in part.caps]
        # valve_nodes = []
        # for cap in caps:
        #     valve_nodes.extend(cap.node_ids)

        for cavity in cavities:
            segs = cavity.surface.triangles

            # # remove segments related to valve nodes
            # for n in valve_nodes:
            #     index = np.argwhere(n == segs)
            #     segs = np.delete(segs, np.array(index)[:, 0], axis=0)

            surface_id = self.get_unique_segmentset_id()
            cavity.surface.id = surface_id
            kw = create_segment_set_keyword(
                segments=segs + 1,
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
                    segments=surface.triangles + 1,
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
            # # WYE: add a node at center of cap
            # # Note: should not be applied in ZeropWriter, it will impact dynain file
            # nid = len(self.model.mesh.nodes) + segid
            # self.kw_database.segment_sets.append(
            #     "*NODE\n{0:8d}{1:16f}{2:16f}{3:16f}".format(nid + 1, *cap.centroid)
            # )
            # nid_x = cap.triangles[0, 0]
            # cap.triangles[:, 0] = nid
            # cap.triangles = np.insert(
            #     cap.triangles, 0, np.array([nid, nid_x, cap.triangles[0, 1]]), axis=0
            # )
            # cap.triangles = np.insert(
            #     cap.triangles, -1, np.array([nid, cap.triangles[-1, -1], nid_x]), axis=0
            # )
            # # END WYE:

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

                if len(node_ids) == 0:
                    LOGGER.debug(
                        "Nodes already used. Skipping node set for {0}".format(
                            part.name + " " + cap.name
                        )
                    )
                    continue

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

    def _update_material_db(self, add_active: bool = True, em_couple: bool = False):
        for part in self.model.parts:
            material = part.meca_material

            if isinstance(material, MAT295):
                if material.active is not None:
                    # obtain ca2+ curve
                    x, y = material.active.ca2_curve.dyna_input

                    cid = self.get_unique_curve_id()
                    curve_kw = create_define_curve_kw(
                        x=x,
                        y=y,
                        curve_name=f"ca2+ of {part.name}",
                        curve_id=cid,
                        lcint=10000,
                    )

                    # curve is written even if not used
                    self.kw_database.material.append(curve_kw)
                    # assign curve id
                    material.active.acid = cid if not em_couple else None

                material_kw = MaterialHGOMyocardium2(
                    id=part.mid, mat=material, ignore_active=not add_active
                )

                self.kw_database.material.append(material_kw)

            elif isinstance(material, NeoHookean):
                material_kw = MaterialNeoHook(
                    mid=part.mid,
                    rho=material.rho,
                    c10=material.c10,
                )
                self.kw_database.material.append(material_kw)

    def _deprecated_update_material_db(self, add_active: bool = True, em_couple: bool = False):
        """Update the database of material keywords."""
        act_curve_id = self.get_unique_curve_id()

        material_settings = copy.deepcopy(self.settings.mechanics.material)
        # NOTE: since we remove units, we don't have to access quantities by <var_name>.m
        material_settings._remove_units()

        if not add_active:
            active_dict = None
        else:
            if em_couple:
                # TODO: hard coded EM coupling parameters
                active_dict = {
                    "actype": 3,
                    "acthr": 2e-4,
                    "ca2ion50": 1e-3,
                    "n": 2,
                    "sigmax": 0.08,
                    "f": 0,
                    "l": 1.9,
                    "eta": 1.45,
                }
            else:
                active_dict = material_settings.myocardium["active"]

        for part in self.model.parts:
            if part.has_fiber:
                if part.is_active:
                    material_kw = _Depracated_MaterialHGOMyocardium(
                        mid=part.mid,
                        iso_user=material_settings.myocardium["isotropic"],
                        anisotropy_user=material_settings.myocardium["anisotropic"],
                        active_user=active_dict,
                    )

                    if not em_couple:
                        material_kw.acid = act_curve_id

                else:
                    material_kw = _Depracated_MaterialHGOMyocardium(
                        mid=part.mid,
                        iso_user=material_settings.myocardium["isotropic"],
                        anisotropy_user=material_settings.myocardium["anisotropic"],
                        active_user=None,
                    )

            else:
                # add isotropic material
                if material_settings.passive["type"] == "NeoHook":
                    # use MAT77H
                    material_kw = MaterialNeoHook(
                        mid=part.mid,
                        rho=material_settings.passive["rho"],
                        c10=material_settings.passive["mu1"] / 2,
                    )
                else:
                    # use MAT295, should have the same behavior
                    material_kw = _Depracated_MaterialHGOMyocardium(
                        mid=part.mid, iso_user=dict(material_settings.passive)
                    )

            self.kw_database.material.append(material_kw)

        # Add Ca2+ curve if necessary
        if add_active and not em_couple:
            kw = self.add_active_curve(act_curve_id, material_settings)
            self.kw_database.material.append(kw)

        return

    @staticmethod
    def add_active_curve(act_curve_id, material_settings):
        """Add Active curve to material database."""
        if material_settings.myocardium["active"]["actype"] == 1:
            time_array, calcium_array = active_curve("constant")
        elif material_settings.myocardium["active"]["actype"] == 3:
            time_array, calcium_array = active_curve("Strocchi2020", endtime=8000)
            # work around with threshold
            calcium_array[1:] += 1e-6

        curve_kw = create_define_curve_kw(
            x=time_array,
            y=calcium_array,
            curve_name="calcium_concentration",
            curve_id=act_curve_id,
            lcint=10000,
        )

        if material_settings.myocardium["active"]["actype"] == 1:
            # x scaling from beat rate
            curve_kw.sfa = material_settings.myocardium["active"]["beat_time"]
            # y scaling from Ca2
            curve_kw.sfo = 4.35  # same with material ca2ionmax
        return curve_kw

    def _add_cap_bc(self, bc_type: str):
        """Add boundary condition to the cap.

        Parameters
        ----------
        bc_type : str
            Boundary condition type. Valid bc's include: ["fix_caps", "springs_caps"].
        """
        bc_settings = self.settings.mechanics.boundary_conditions

        valid_bcs = ["fix_caps", "springs_caps"]
        if bc_type not in valid_bcs:
            raise ValueError("Cap/Valve boundary condition must be of type: %r" % valid_bcs)

        # create list of cap names where to add the spring b.c
        caps_to_use = []
        if isinstance(self.model, LeftVentricle):
            caps_to_use = [
                "mitral-valve",
                "aortic-valve",
            ]
        elif isinstance(self.model, BiVentricle):
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
            if isinstance(self, ZeroPressureMechanicsDynaWriter):
                # add additional constraint to avoid rotation
                caps_to_use.extend(["pulmonary-valve"])

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

            spring_stiffness = bc_settings.valve["stiffness"].m
            if type(self) == ZeroPressureMechanicsDynaWriter:
                spring_stiffness *= 1e16

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
                if cap.name in caps_to_use:
                    self.kw_database.boundary_conditions.append(f"$$ spring at {cap.name}$$")
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

        attached_nodes = cap.node_ids

        # use pre-computed nodal area
        nodal_areas = self.model.mesh.point_data["nodal_areas"][attached_nodes]

        # scaled spring stiffness by nodal area
        scale_factor_normal *= nodal_areas
        scale_factor_radial *= nodal_areas

        # add sd_orientiation, element discrete

        # compute the radial components
        sd_orientations_radial = self.model.mesh.nodes[attached_nodes, :] - cap.centroid

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
        epicardium_nodes, point_normal = self._get_epicardium_nodes(apply="ventricle")

        # use pre-computed nodal areas
        nodal_areas = self.model.mesh.point_data["nodal_areas"][epicardium_nodes]
        # penalty
        penalty_function = self._get_longitudinal_penalty(robin_settings["ventricle"])
        nodal_penalty = penalty_function[epicardium_nodes]

        # compute scale factor
        scale_factors = nodal_areas * nodal_penalty

        # debug = pv.PointGrid(self.model.mesh.points[epicardium_nodes])
        # debug.point_data["stiff"] = scale_factors
        # debug.point_data["normal"] = point_normal
        # debug.save("pericardium.vtk")

        k = scale * robin_settings["ventricle"]["stiffness"].to("MPa/mm").m
        mask = penalty_function[epicardium_nodes] > 0.0001
        self._write_discret_elements(
            "spring", k, epicardium_nodes[mask], point_normal[mask], scale_factors[mask]
        )
        dc = robin_settings["ventricle"]["damper"].to("MPa/mm*ms").m
        self._write_discret_elements("damper", dc, epicardium_nodes, point_normal, nodal_areas)

        if isinstance(self.model, FourChamber):
            epicardium_nodes, point_normal = self._get_epicardium_nodes(apply="atrial")
            nodal_areas = self.model.mesh.point_data["nodal_areas"][epicardium_nodes]
            k = robin_settings["atrial"]["stiffness"].to("MPa/mm").m
            self._write_discret_elements("spring", k, epicardium_nodes, point_normal, nodal_areas)
            dc = robin_settings["atrial"]["damper"].to("MPa/mm*ms").m
            self._write_discret_elements("damper", dc, epicardium_nodes, point_normal, nodal_areas)
        return

    def _get_epicardium_nodes(self, apply: Literal["ventricle", "atrial"] = "ventricle"):
        """Extract epicardium nodes to apply Robin BC.

        Parameters
        ----------
        apply : Literal[&quot;ventricle&quot;, &quot;atrial&quot;], optional
            Apply to which part, by default "ventricle"

        TODO: move to model

        Returns
        -------
            node list and their normal vectors
        """
        epicardium_faces = np.empty((0, 3), dtype=int)

        LOGGER.debug(f"Collecting epicardium nodesets of {apply}:")
        if apply == "ventricle":
            targets = [part for part in self.model.parts if "ventricle" in part.name]
        elif apply == "atrial":
            targets = [part for part in self.model.parts if "atrium" in part.name]

        epicardium_surfaces = [ventricle.epicardium for ventricle in targets]

        for surface in epicardium_surfaces:
            epicardium_faces = np.vstack([epicardium_faces, surface.triangles])

        # some nodes on the edge must be included
        epicardium_nodes, a = np.unique(epicardium_faces, return_inverse=True)

        # build pericardium polydata
        coord = self.model.mesh.nodes[epicardium_nodes]
        connect = a.reshape(epicardium_faces.shape)
        pericardium_polydata = vtkmethods.create_vtk_surface_triangles(coord, connect, clean=False)
        # vtkmethods.write_vtkdata_to_vtkfile(pericardium_polydata,'pericardium.vtk')

        # compute normal
        _, point_normal = vtkmethods.add_normals_to_polydata(
            pericardium_polydata, return_normals=True
        )
        return epicardium_nodes, point_normal

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
        uvc_l = self.model.mesh.point_data["uvc_longitudinal"]
        if np.any(uvc_l < 0):
            LOGGER.warning(
                "Negative normalized longitudinal coordinate detected."
                "Changing {0} negative uvc_l values to 1".format(np.sum((uvc_l < 0))),
            )
        uvc_l[uvc_l < 0] = 1

        penalty_function = -_sigmoid((abs(uvc_l) - penalty_c0) * penalty_c1) + 1
        return penalty_function

    def _write_discret_elements(
        self, type: Literal["spring", "damper"], global_fact, nodes, directions, local_stiff
    ):
        # create unique ids for keywords
        part_id = self.get_unique_part_id()
        section_id = self.get_unique_section_id()
        mat_id = self.get_unique_mat_id()

        # define material
        if type == "spring":
            mat_kw = keywords.MatSpringElastic(mid=mat_id, k=global_fact)
        elif type == "damper":
            mat_kw = keywords.MatDamperViscous(mid=mat_id, dc=global_fact)

        # define part
        part_kw = keywords.Part()
        part_kw.parts = pd.DataFrame(
            {"heading": [f"{type}"], "pid": [part_id], "secid": [section_id], "mid": [mat_id]}
        )
        # define section
        section_kw = keywords.SectionDiscrete(secid=section_id, cdl=0, tdl=0)

        # define spring orientations
        sd_orientation_kw = create_define_sd_orientation_kw(
            vectors=directions, vector_id_offset=self.id_offset["vector"]
        )
        # add offset
        self.id_offset["vector"] = sd_orientation_kw.vectors["vid"].to_numpy()[-1]
        vector_ids = sd_orientation_kw.vectors["vid"].to_numpy().astype(int)

        # define spring nodes
        nodes_table = np.vstack([nodes + 1, np.zeros(len(nodes))]).T

        # create discrete elements
        discrete_element_kw = create_discrete_elements_kw(
            nodes=nodes_table,
            part_id=part_id,
            vector_ids=vector_ids,
            scale_factor=local_stiff,
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

    def _update_cap_elements_db(self, add_mesh=True):
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
        material_settings = copy.deepcopy(self.settings.mechanics.material)
        material_settings._remove_units()

        def _add_linear_constraint(id: int, slave_id: int, master_ids: List[int]) -> list:
            lin_constraint_kws = []

            for dof in range(1, 4):
                kw = custom_keywords.ConstrainedLinearGlobal(licd=3 * id + dof)
                data = np.empty((0, 3))
                data = np.vstack([data, [slave_id, dof, 1.0]])

                for m_id in master_ids:
                    data = np.vstack([data, [m_id, dof, -1 / len(master_ids)]])

                kw.linear_constraints = pd.DataFrame(columns=["nid", "dof", "coef"], data=data)

                lin_constraint_kws.append(kw)

            return lin_constraint_kws

        # caps are rigid in zerop
        if type(self) == ZeroPressureMechanicsDynaWriter:
            material_kw = keywords.MatRigid(
                mid=mat_null_id,
                ro=material_settings.cap["rho"],
                e=1.0,  # MPa
            )

        else:
            if material_settings.cap["type"] == "stiff":
                material_kw = MaterialNeoHook(
                    mid=mat_null_id,
                    rho=material_settings.cap["rho"],
                    c10=material_settings.cap["mu1"] / 2,
                )

            elif material_settings.cap["type"] == "null":
                material_kw = keywords.MatNull(
                    mid=mat_null_id,
                    ro=material_settings.cap["rho"],
                )
            elif material_settings.cap["type"] == "rigid":
                material_kw = keywords.MatRigid(
                    mid=mat_null_id,
                    ro=material_settings.cap["rho"],
                    e=1.0,  # MPa
                )

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
            part_info = pd.DataFrame(
                {
                    "heading": [cap.name],
                    "pid": [cap.pid],
                    "secid": [section_id],
                    "mid": [mat_null_id],
                }
            )
            part_kw.parts = part_info

            if cap.centroid is not None:
                # cap centroids already added to mesh for v0.2
                if heart_version == "v0.1":
                    if add_mesh:
                        # Add center node
                        node_kw = keywords.Node()
                        df = pd.DataFrame(
                            data=np.insert(cap.centroid, 0, cap.centroid_id + 1).reshape(1, -1),
                            columns=node_kw.nodes.columns[0:4],
                        )
                        node_kw.nodes = df
                        # comment the line '*NODE' so nodes.k can be parsed by zerop solver
                        # correctly otherwise, these nodes will not be updated in iterations
                        s = "$" + node_kw.write()
                        self.kw_database.nodes.append(s)

                if type(self) == MechanicsDynaWriter:
                    # center node constraint: average of edge nodes
                    n = len(cap.node_ids) // 7  # select n+1 node for interpolation
                    constraint_list = _add_linear_constraint(
                        len(cap_names_used), cap.centroid_id + 1, cap.node_ids[::n] + 1
                    )
                    self.kw_database.cap_elements.extend(constraint_list)

                # # # do not work with mpp
                # # constraint = keywords.ConstrainedInterpolation(
                # #     icid=len(cap_names_used) + 1,
                # #     dnid=cap.centroid_id + 1,
                # #     ddof=123,
                # #     ityp=1,
                # #     fgm=1,
                # #     inid=cap.nsid,
                # #     idof=123,
                # # )
                # # self.kw_database.cap_elements.append(constraint)

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
            if add_mesh:
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

        system_settings = copy.deepcopy(self.settings.mechanics.system)
        system_settings._remove_units()

        # closed loop uses a custom executable
        if system_settings.name == "ClosedLoop":
            raise NotImplementedError("Closed loop circulation not yet supported.")
            LOGGER.warning(
                "Note that this model type requires a custom executable that "
                "supports the Closed Loop circulation model!"
            )
            if isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
                file_path = resource_path(
                    "ansys.heart.writer", "templates/system_model_settings_bv.json"
                ).__enter__()
            elif isinstance(self.model, LeftVentricle):
                file_path = resource_path(
                    "ansys.heart.writer", "templates/system_model_settings_lv.json"
                ).__enter__()

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
        elif system_settings.name == "ConstantPreloadWindkesselAfterload":
            if self.system_model_name != system_settings.name:
                LOGGER.error("Circulation system parameters cannot be rad from Json")

            for cavity in self.model.cavities:
                if "Left ventricle" in cavity.name:
                    define_function_wk = define_function_windkessel(
                        function_id=10,
                        function_name="constant_preload_windkessel_afterload_left",
                        implicit=True,
                        constants=dict(system_settings.left_ventricle["constants"]),
                        initialvalues=system_settings.left_ventricle["initial_value"]["part"],
                    )
                    self.kw_database.control_volume.append(define_function_wk)

                elif "Right ventricle" in cavity.name:
                    define_function_wk = define_function_windkessel(
                        function_id=11,
                        function_name="constant_preload_windkessel_afterload_right",
                        implicit=True,
                        constants=dict(system_settings.right_ventricle["constants"]),
                        initialvalues=system_settings.right_ventricle["initial_value"]["part"],
                    )
                    self.kw_database.control_volume.append(define_function_wk)

        return

    def _add_enddiastolic_pressure_bc2(self, pressure_lv: float = 1, pressure_rv: float = 1):
        """
        Apply ED pressure by control volume.

        Notes
        -----
        LSDYNA stress reference configuration bug with this load due to define function.
        """
        cavities = [part.cavity for part in self.model.parts if part.cavity]
        for cavity in cavities:
            if "atrium" in cavity.name:
                continue

            # create CV
            cv_kw = keywords.DefineControlVolume()
            cv_kw.id = cavity.surface.id
            cv_kw.sid = cavity.surface.id
            self.kw_database.main.append(cv_kw)

            # define CV interaction
            cvi_kw = keywords.DefineControlVolumeInteraction()
            cvi_kw.id = cavity.surface.id
            cvi_kw.cvid1 = cavity.surface.id
            cvi_kw.cvid2 = 0  # ambient

            if "Left ventricle" in cavity.name:
                cvi_kw.lcid_ = 10
                pressure = pressure_lv
            elif "Right ventricle" in cavity.name:
                cvi_kw.lcid_ = 11
                pressure = pressure_rv

            self.kw_database.main.append(cvi_kw)

            # define define function
            definefunction_str = _ed_load_template()
            self.kw_database.main.append(
                definefunction_str.format(
                    cvi_kw.lcid_, "flow_" + cavity.name.replace(" ", "_"), pressure, -200
                )
            )

        self.kw_database.main.append(keywords.DatabaseIcvout(dt=10, binary=2))
        return

    def _add_enddiastolic_pressure_bc(self, pressure_lv: float = 1, pressure_rv: float = 1):
        """Add end diastolic pressure boundary condition on the left and right endocardium."""
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
            if cavity.name == "Left ventricle":
                load = keywords.LoadSegmentSet(
                    ssid=cavity.surface.id, lcid=load_curve_id, sf=pressure_lv
                )
                self.kw_database.main.append(load)
            elif cavity.name == "Right ventricle":
                load = keywords.LoadSegmentSet(
                    ssid=cavity.surface.id, lcid=load_curve_id, sf=pressure_rv
                )
                self.kw_database.main.append(load)
            elif cavity.name == "Left atrium":
                load = keywords.LoadSegmentSet(
                    ssid=cavity.surface.id, lcid=load_curve_id, sf=pressure_lv
                )
                self.kw_database.main.append(load)
            elif cavity.name == "Right atrium":
                load = keywords.LoadSegmentSet(
                    ssid=cavity.surface.id, lcid=load_curve_id, sf=pressure_rv
                )
                self.kw_database.main.append(load)
            else:
                continue

        # # load only endocardium segment (exclude cap shells)
        # for part in self.model.parts:
        #     for surface in part.surfaces:
        #         if surface.name == "Left ventricle endocardium":
        #             scale_factor = pressure_lv
        #             seg_id = surface.id
        #             load = keywords.LoadSegmentSet(
        #                 ssid=seg_id, lcid=load_curve_id, sf=scale_factor
        #             )
        #             self.kw_database.main.append(load)
        #         elif (
        #             surface.name == "Right ventricle endocardium"
        #             or surface.name == "Right ventricle endocardium septum"
        #         ):
        #             scale_factor = pressure_rv
        #             seg_id = surface.id
        #             load = keywords.LoadSegmentSet(
        #                 ssid=seg_id, lcid=load_curve_id, sf=scale_factor
        #             )
        #             self.kw_database.main.append(load)
        return

    def _add_constant_atrial_pressure(self, pressure_lv: float = 1, pressure_rv: float = 1):
        """Missing circulation model for atrial cavity, apply constant ED pressure."""
        # create unit load curve
        load_curve_id = self.get_unique_curve_id()
        load_curve_kw = create_define_curve_kw(
            [0, 1e20], [1.0, 1.0], "constant load curve", load_curve_id, 100
        )

        # append unit curve to main.k
        self.kw_database.main.append(load_curve_kw)

        # create *LOAD_SEGMENT_SETS for each ventricular cavity
        cavities = [part.cavity for part in self.model.parts if part.cavity]
        for cavity in cavities:
            if cavity.name == "Left atrium":
                load = keywords.LoadSegmentSet(
                    ssid=cavity.surface.id, lcid=load_curve_id, sf=pressure_lv
                )
                self.kw_database.main.append(load)
            elif cavity.name == "Right atrium":
                load = keywords.LoadSegmentSet(
                    ssid=cavity.surface.id, lcid=load_curve_id, sf=pressure_rv
                )
                self.kw_database.main.append(load)


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

    def __init__(
        self,
        model: HeartModel,
        settings: SimulationSettings = None,
    ) -> None:
        super().__init__(model=model, settings=settings)

        self.kw_database = MechanicsDecks()
        """Collection of keyword decks relevant for mechanics."""

        return

    def update(self):
        """Update the keyword database."""
        bc_settings = self.settings.mechanics.boundary_conditions

        self._update_main_db()

        self.kw_database.main.title = self.model.info.model_type + " zero-pressure"

        self._update_node_db()
        self._update_parts_db()
        self._update_solid_elements_db(add_fibers=True)
        self._update_segmentsets_db()
        self._update_nodesets_db()
        self._update_material_db(add_active=False)
        if self.cap_in_zerop:
            # define cap element
            self._update_cap_elements_db()

        # TODO: it should be after cap creation, or it will be written in dynain
        # for boundary conditions
        self._add_cap_bc(bc_type="springs_caps")
        if isinstance(self.model, FourChamber):
            # add a small constraint to avoid rotation
            self._add_pericardium_bc(scale=0.01)

        # # Approximate end-diastolic pressures
        pressure_lv = bc_settings.end_diastolic_cavity_pressure["left_ventricle"].m
        pressure_rv = bc_settings.end_diastolic_cavity_pressure["right_ventricle"].m
        self._add_enddiastolic_pressure_bc(pressure_lv=pressure_lv, pressure_rv=pressure_rv)

        # zerop key words
        self._add_control_reference_configuration()

        # export dynain file
        save_part_ids = []
        for part in self.model.parts:
            save_part_ids.append(part.pid)

        caps = [cap for part in self.model.parts for cap in part.caps]
        for cap in caps:
            if cap.pid != None:  # MV,TV for atrial parts get None
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

    def export(self, export_directory: str):
        """Write the model to files."""
        super().export(export_directory)

        # export segment sets to separate file
        self._export_cavity_segmentsets(export_directory)

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
            keywords.ControlImplicitGeneral(imflag=1, dt0=settings.analysis.dtmin)
        )

        # add implicit solution controls
        self.kw_database.main.append(
            keywords.ControlImplicitSolution(
                maxref=35,
                dctol=0.01,
                ectol=1e6,
                rctol=1e3,
                abstol=-1e-20,
                dnorm=1,
                diverg=2,
                lsmtd=5,
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


class FiberGenerationDynaWriter(BaseDynaWriter):
    """Class for preparing the input for a fiber-generation LS-DYNA simulation."""

    def __init__(self, model: HeartModel, settings: SimulationSettings = None) -> None:
        super().__init__(model=model, settings=settings)
        self.kw_database = FiberGenerationDecks()
        """Collection of keywords relevant for fiber generation."""

    def update(self):
        """Update keyword database for Fiber generation: overwrites the inherited function."""
        ##
        self._update_main_db()  # needs updating

        if isinstance(self.model, (FourChamber, FullHeart)):
            LOGGER.warning(
                "Atrium present in the model, they will be removed for ventricle fiber generation."
            )

            parts = [
                part
                for part in self.model.parts
                if part.part_type == "ventricle" or part.part_type == "septum"
            ]
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

        self._update_segmentsets_db()
        self._update_nodesets_db(remove_one_node_from_cell=True)

        # # update ep settings
        self._update_ep_settings()
        self._update_create_fibers()

        self._get_list_of_includes()
        self._add_includes()

        return

    def _remove_atrial_nodes_from_ventricles_surfaces(self):
        """Remove nodes other than ventricular from ventricular surfaces."""
        parts = [
            part
            for part in self.model.parts
            if part.part_type == "ventricle" or part.part_type == "septum"
        ]

        tet_ids = np.empty((0), dtype=int)
        for part in parts:
            tet_ids = np.append(tet_ids, part.element_ids)
            tets = self.model.mesh.tetrahedrons[tet_ids, :]
        nids = np.unique(tets)

        for part in parts:
            for surface in part.surfaces:
                nodes_to_remove = surface.node_ids[
                    np.isin(surface.node_ids, nids, assume_unique=True, invert=True)
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
            element_ids = part.element_ids
            em_mat_id = self.get_unique_mat_id()
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
        node_sets_ids_endo = []  # relevant for both models
        node_sets_ids_epi = []  # relevant for both models
        node_set_ids_epi_and_rseptum = []  # only relevant for bv, 4c and full model

        # list of ventricular parts
        ventricles = [part for part in self.model.parts if "ventricle" in part.name]
        septum = self.model.get_part("Septum")

        # collect node set ids (already generated previously)
        node_sets_ids_epi = [ventricle.epicardium.nsid for ventricle in ventricles]
        node_sets_ids_endo = []
        for ventricle in ventricles:
            for surface in ventricle.surfaces:
                if "endocardium" in surface.name:
                    node_sets_ids_endo.append(surface.nsid)

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
            from ansys.heart.writer.define_function_strings import (
                function_alpha,
                function_beta,
                function_beta_septum,
            )

            self.kw_database.create_fiber.append(
                keywords.DefineFunction(
                    fid=101,
                    function=function_alpha(
                        alpha_endo=self.settings.fibers.alpha_endo.m,
                        alpha_epi=self.settings.fibers.alpha_epi.m,
                    ),
                )
            )
            self.kw_database.create_fiber.append(
                keywords.DefineFunction(
                    fid=102,
                    function=function_beta(
                        beta_endo=self.settings.fibers.beta_endo.m,
                        beta_epi=self.settings.fibers.beta_epi.m,
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
            from ansys.heart.writer.define_function_strings import (
                function_alpha,
                function_beta,
                function_beta_septum,
            )

            self.kw_database.create_fiber.append(
                keywords.DefineFunction(
                    fid=101,
                    function=function_alpha(
                        alpha_endo=self.settings.fibers.alpha_endo.m,
                        alpha_epi=self.settings.fibers.alpha_epi.m,
                    ),
                )
            )
            self.kw_database.create_fiber.append(
                keywords.DefineFunction(
                    fid=102,
                    function=function_beta(
                        beta_endo=self.settings.fibers.beta_endo.m,
                        beta_epi=self.settings.fibers.beta_epi.m,
                    ),
                )
            )
            self.kw_database.create_fiber.append(
                keywords.DefineFunction(
                    fid=103,
                    function=function_beta_septum(
                        beta_endo=self.settings.fibers.beta_endo_septum.m,
                        beta_epi=self.settings.fibers.beta_epi_septum.m,
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


# todo: why it's from MechanicsDynaWriter not BaseWriter?
class PurkinjeGenerationDynaWriter(MechanicsDynaWriter):
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
                "they will be removed for ventricle Purkinje generation."
            )
            self._keep_ventricles()

        self._update_parts_db()  # can stay the same (could move to base class++++++++++++++++++++)
        self._update_solid_elements_db(add_fibers=False)
        self._update_material_db()

        self._update_segmentsets_db()  # can stay the same
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
                element_id = np.argwhere(np.any(endocardium.triangles == node_apex_left, axis=1))[
                    0
                ][0]

                node_apex_left = endocardium.triangles[element_id, :][
                    np.argwhere(
                        np.isin(
                            endocardium.triangles[element_id, :],
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

            node_set_id_apex_left = self.get_unique_nodeset_id()
            # create node-sets for apex
            node_set_apex_kw = create_node_set_keyword(
                node_ids=[node_apex_left + 1],
                node_set_id=node_set_id_apex_left,
                title="apex node left",
            )

            self.kw_database.node_sets.append(node_set_apex_kw)

            apex_left_coordinates = self.model.mesh.nodes[node_apex_left, :]

            node_id_start_left = self.model.mesh.nodes.shape[0] + 1

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
                    iedgeid=edge_id_start_left,  # TODO check if beam elements exist in mesh
                    pmjtype=self.settings.purkinje.pmjtype.m,
                    pmjradius=self.settings.purkinje.pmjradius.m,
                    pmjrestype=self.settings.electrophysiology.material.beam["pmjrestype"].m,
                    pmjres=self.settings.electrophysiology.material.beam["pmjres"].m,
                )
            )

        # Add right purkinje only in biventricular or 4chamber models
        if isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
            node_apex_right = self.model.right_ventricle.apex_points[0].node_id
            segment_set_ids_endo_right = self.model.right_ventricle.endocardium.id

            # check whether point is on edge of endocardium - otherwise pick another node in
            # the same triangle
            endocardium = self.model.right_ventricle.endocardium
            endocardium.get_boundary_edges()
            if np.any(endocardium.boundary_edges == node_apex_right):
                element_id = np.argwhere(np.any(endocardium.triangles == node_apex_right, axis=1))[
                    0
                ][0]

                node_apex_right = endocardium.triangles[element_id, :][
                    np.argwhere(
                        np.isin(
                            endocardium.triangles[element_id, :],
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
                    inodeid=node_id_start_right,  # TODO check if beam elements exist in mesh
                    iedgeid=edge_id_start_right,
                    pmjtype=self.settings.purkinje.pmjtype.m,
                    pmjradius=self.settings.purkinje.pmjradius.m,
                    pmjrestype=self.settings.electrophysiology.material.beam["pmjrestype"].m,
                    pmjres=self.settings.electrophysiology.material.beam["pmjres"].m,
                )
            )

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


class ElectrophysiologyDynaWriter(BaseDynaWriter):
    """Class for preparing the input for an Electrophysiology LS-DYNA simulation."""

    def __init__(self, model: HeartModel, settings: SimulationSettings = None) -> None:
        if isinstance(model, FourChamber):
            model._create_isolation_part()
        if model.info.add_blood_pool == True:
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

        self._update_segmentsets_db()

        # TODO check if no existing node set ids conflict with surface ids
        # For now, new node sets should be created after calling
        # self._update_nodesets_db()
        self._update_nodesets_db()
        self._update_cellmodels()

        if self.model.beam_network:
            # with smcoupl=1, coupling is disabled
            self.kw_database.ep_settings.append(keywords.EmControlCoupling(smcoupl=1))
            self._update_use_Purkinje()

        # update ep settings
        self._update_ep_settings()
        self._update_stimulation()

        if self.model.info.add_blood_pool == True:
            self._update_blood_settings()

        if hasattr(self.model, "electrodes") and len(self.model.electrodes) != 0:
            self._update_ECG_coordinates()

        self._get_list_of_includes()
        self._add_includes()

        return

    def _isolate_atria_and_ventricles(self):
        """Add duplicate nodes between atria and ventricles."""
        if isinstance(self.model, (FourChamber, FullHeart)):
            # 1,3 - 2,4
            self.model.mesh.establish_connectivity()
            left_ventricle_left_atrium = []
            right_ventricle_right_atrium = []
            left_ventricle_right_atrium = []
            right_ventricle_left_atrium = []
            left_ventricle_left_atrium_name = "left-ventricle_left-atrium"
            right_ventricle_right_atrium_name = "right-ventricle_right-atrium"
            left_ventricle_right_atrium_name = "left-ventricle_right-atrium"
            right_ventricle_left_atrium_name = "right-ventricle_left-atrium"

            # build atrio-ventricular tag-id pairs
            # labels_to_ids stores the mapping between tag-ids and the corresponding label.
            labels_to_tag_ids = self.model.info.labels_to_ids
            left_ventricle_left_atrium = [
                labels_to_tag_ids["Left ventricle myocardium"],
                labels_to_tag_ids["Left atrium myocardium"],
            ]
            right_ventricle_right_atrium = [
                labels_to_tag_ids["Right ventricle myocardium"],
                labels_to_tag_ids["Right atrium myocardium"],
            ]
            left_ventricle_right_atrium = [
                labels_to_tag_ids["Left ventricle myocardium"],
                labels_to_tag_ids["Right atrium myocardium"],
            ]
            right_ventricle_left_atrium = [
                labels_to_tag_ids["Right ventricle myocardium"],
                labels_to_tag_ids["Left atrium myocardium"],
            ]

            # build atrioventricular tag_id pairs
            left_ventricle_left_atrium = np.unique(left_ventricle_left_atrium)
            right_ventricle_right_atrium = np.unique(right_ventricle_right_atrium)
            left_ventricle_right_atrium = np.unique(left_ventricle_right_atrium)
            right_ventricle_left_atrium = np.unique(right_ventricle_left_atrium)
            # find atrioventricular shared nodes/interfaces
            self.model.mesh.add_interfaces(
                [
                    left_ventricle_left_atrium,
                    right_ventricle_right_atrium,
                    left_ventricle_right_atrium,
                    right_ventricle_left_atrium,
                ],
                [
                    left_ventricle_left_atrium_name,
                    right_ventricle_right_atrium_name,
                    left_ventricle_right_atrium_name,
                    right_ventricle_left_atrium_name,
                ],
            )

            # duplicate nodes of each interface in atrium side
            old_nodes = []
            new_nodes = []
            for interface in self.model.mesh.interfaces:
                if interface.name != None and interface.name == left_ventricle_left_atrium_name:
                    interface_nids = interface.node_ids
                    tets_atrium = self.model.mesh.tetrahedrons[
                        self.model.left_atrium.element_ids, :
                    ]

                    nids_tobe_replaced = tets_atrium[np.isin(tets_atrium, interface_nids)]
                    new_nids = np.array(
                        list(
                            map(
                                lambda id: (np.where(interface_nids == id))[0][0],
                                nids_tobe_replaced,
                            )
                        )
                    ) + len(self.model.mesh.nodes)
                    tets_atrium[np.isin(tets_atrium, interface_nids)] = new_nids
                    old_nodes.extend(nids_tobe_replaced)
                    new_nodes.extend(new_nids)

                    # TODO refactor this and avoid big ndarray copies
                    tets: np.ndarray = self.model.mesh.tetrahedrons
                    tets[self.model.left_atrium.element_ids, :] = tets_atrium

                    self.model.mesh.tetrahedrons = tets
                    self.model.mesh.nodes = np.append(
                        self.model.mesh.nodes, self.model.mesh.nodes[interface_nids, :], axis=0
                    )

                elif interface.name != None and interface.name == right_ventricle_right_atrium_name:
                    interface_nids = interface.node_ids
                    tets_atrium = self.model.mesh.tetrahedrons[
                        self.model.right_atrium.element_ids, :
                    ]

                    nids_tobe_replaced = tets_atrium[np.isin(tets_atrium, interface_nids)]
                    new_nids = np.array(
                        list(
                            map(
                                lambda id: (np.where(interface_nids == id))[0][0],
                                nids_tobe_replaced,
                            )
                        )
                    ) + len(self.model.mesh.nodes)
                    tets_atrium[np.isin(tets_atrium, interface_nids)] = new_nids
                    old_nodes.extend(nids_tobe_replaced)
                    new_nodes.extend(new_nids)
                    tets: np.ndarray = self.model.mesh.tetrahedrons
                    tets[self.model.right_atrium.element_ids, :] = tets_atrium

                    self.model.mesh.tetrahedrons = tets

                    self.model.mesh.nodes = np.append(
                        self.model.mesh.nodes, self.model.mesh.nodes[interface_nids, :], axis=0
                    )
                elif interface.name != None and interface.name == left_ventricle_right_atrium_name:
                    interface_nids = interface.node_ids
                    tets_atrium = self.model.mesh.tetrahedrons[
                        self.model.right_atrium.element_ids, :
                    ]

                    nids_tobe_replaced = tets_atrium[np.isin(tets_atrium, interface_nids)]
                    new_nids = np.array(
                        list(
                            map(
                                lambda id: (np.where(interface_nids == id))[0][0],
                                nids_tobe_replaced,
                            )
                        )
                    ) + len(self.model.mesh.nodes)
                    tets_atrium[np.isin(tets_atrium, interface_nids)] = new_nids
                    old_nodes.extend(nids_tobe_replaced)
                    new_nodes.extend(new_nids)

                    tets: np.ndarray = self.model.mesh.tetrahedrons
                    tets[self.model.right_atrium.element_ids, :] = tets_atrium

                    self.model.mesh.tetrahedrons = tets

                    self.model.mesh.nodes = np.append(
                        self.model.mesh.nodes, self.model.mesh.nodes[interface_nids, :], axis=0
                    )
                elif interface.name != None and interface.name == right_ventricle_left_atrium_name:
                    interface_nids = interface.node_ids
                    tets_atrium = self.model.mesh.tetrahedrons[
                        self.model.left_atrium.element_ids, :
                    ]

                    nids_tobe_replaced = tets_atrium[np.isin(tets_atrium, interface_nids)]
                    new_nids = np.array(
                        list(
                            map(
                                lambda id: (np.where(interface_nids == id))[0][0],
                                nids_tobe_replaced,
                            )
                        )
                    ) + len(self.model.mesh.nodes)
                    tets_atrium[np.isin(tets_atrium, interface_nids)] = new_nids
                    old_nodes.extend(nids_tobe_replaced)
                    new_nodes.extend(new_nids)

                    tets: np.ndarray = self.model.mesh.tetrahedrons
                    tets[self.model.left_atrium.element_ids, :] = tets_atrium

                    self.model.mesh.tetrahedrons = tets

                    self.model.mesh.nodes = np.append(
                        self.model.mesh.nodes, self.model.mesh.nodes[interface_nids, :], axis=0
                    )

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
        for part in self.model.parts:
            self.kw_database.material.append(f"$$ {part.name} $$")
            partname = part.name.lower()
            if ("atrium" in partname) or ("ventricle" in partname) or ("septum" in partname):
                # Electrically "active" tissue (mtype=1)
                ep_mid = part.pid
                self.kw_database.material.append(
                    custom_keywords.EmMat003(
                        mid=ep_mid,
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
                )
            elif "isolation" in partname:
                # assign insulator material to isolation layer.
                ep_mid = part.pid
                self.kw_database.material.append(
                    custom_keywords.EmMat001(mid=ep_mid, mtype=1, sigma=1),
                )
            else:
                # Electrically non-active tissue (mtype=4)
                # These bodies are still conductive bodies
                # in the extra-cellular space
                ep_mid = part.pid
                self.kw_database.material.append(
                    custom_keywords.EmMat001(
                        mid=ep_mid,
                        mtype=4,
                        sigma=self.settings.electrophysiology.material.myocardium[
                            "sigma_passive"
                        ].m,
                    ),
                )

        return

    def _update_cellmodels(self):
        """Add cell model for each defined part."""
        for part in self.model.parts:
            partname = part.name.lower()
            if ("atrium" in partname) or ("ventricle" in partname) or ("septum" in partname):
                ep_mid = part.pid
                # One cell model for myocardium, default value is epi layer parameters
                self._add_Tentusscher_keyword_epi(matid=ep_mid)

        # different cell models for endo/mid/epi layer
        # TODO:  this will override previous definition?
        #        what's the situation at setptum? and at atrial?
        if "transmural" in self.model.mesh.point_data.keys():
            (
                endo_id,
                mid_id,
                epi_id,
            ) = self._create_myocardial_nodeset_layers()
            self._add_Tentusscher_keyword_endo(matid=-endo_id)
            self._add_Tentusscher_keyword_mid(matid=-mid_id)
            self._add_Tentusscher_keyword_epi(matid=-epi_id)

    def _create_myocardial_nodeset_layers(
        self, percent_endo=0.17, percent_mid=0.41, percent_epi=0.42
    ):
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

    def _add_Tentusscher_keyword_epi(self, matid):
        cell_kw = keywords.EmEpCellmodelTentusscher(
            mid=matid,
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
        )
        cell_kw.gas_constant = 8314.472
        cell_kw.faraday_constant = 96485.3415
        self.kw_database.cell_models.append(cell_kw)
        return

    def _add_Tentusscher_keyword_endo(self, matid):
        cell_kw = keywords.EmEpCellmodelTentusscher(
            mid=matid,
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
            gto=0.073,
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
            v=-86.709,
            ki=138.4,
            nai=10.355,
            cai=0.00013,
            cass=0.00036,
            casr=3.715,
            rpri=0.9068,
            xr1=0.00448,
            xr2=0.476,
            xs=0.0087,
            m=0.00155,
            h=0.7573,
            j=0.7225,
            d=3.164e-5,
            f=0.8009,
            f2=0.9778,
            fcass=0.9953,
            s=0.3212,
            r=2.235e-8,
        )
        cell_kw.gas_constant = 8314.472
        cell_kw.faraday_constant = 96485.3415
        self.kw_database.cell_models.append(cell_kw)
        return

    def _add_Tentusscher_keyword_mid(self, matid):
        cell_kw = keywords.EmEpCellmodelTentusscher(
            mid=matid,
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
            gks=0.098,
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
            vleak=0.00042,
            vxfer=0.0038,
            vmaxup=0.006375,
            kup=0.00025,
            bufc=0.2,
            kbufc=0.001,
            bufsr=10.0,
            kbufsf=0.3,
            bufss=0.4,
            kbufss=0.00025,
            v=-85.423,
            ki=138.52,
            nai=10.132,
            cai=0.000153,
            cass=0.00036,
            casr=4.272,
            rpri=0.8978,
            xr1=0.0165,
            xr2=0.473,
            xs=0.0174,
            m=0.00165,
            h=0.749,
            j=0.6788,
            d=3.288e-5,
            f=0.7026,
            f2=0.9526,
            fcass=0.9942,
            s=0.999998,
            r=2.347e-8,
        )
        cell_kw.gas_constant = 8314.472
        cell_kw.faraday_constant = 96485.3415
        self.kw_database.cell_models.append(cell_kw)
        return

    def _update_ep_settings(self):
        """Add the settings for the electrophysiology solver."""
        self.kw_database.ep_settings.append(
            keywords.EmControl(
                emsol=11, numls=4, macrodt=1, dimtype=None, nperio=None, ncylbem=None
            )
        )

        self.kw_database.ep_settings.append(
            custom_keywords.EmEpIsoch(idisoch=1, idepol=1, dplthr=-20, irepol=1, rplthr=-40)
        )

        # use defaults
        self.kw_database.ep_settings.append(custom_keywords.EmControlEp(numsplit=1))

        self.kw_database.ep_settings.append(
            keywords.EmSolverFem(reltol=1e-6, maxite=int(1e4), precon=2)
        )

        self.kw_database.ep_settings.append(keywords.EmOutput(mats=1, matf=1, sols=1, solf=1))

    def _update_stimulation(self):
        # # define apex node set
        # node_apex_left = self.model.left_ventricle.apex_points[0].node_id
        # node_set_id_apex_left = self.get_unique_nodeset_id()

        # # create node-sets for apex left
        # node_set_kw = create_node_set_keyword(
        #     node_ids=[node_apex_left + 1],
        #     node_set_id=node_set_id_apex_left,
        #     title="apex node left",
        # )
        # self.kw_database.node_sets.append(node_set_kw)

        # # right ventricle apex
        # if not isinstance(self.model, LeftVentricle):
        #     node_apex_right = self.model.right_ventricle.apex_points[0].node_id
        #     node_set_id_apex_right = self.get_unique_nodeset_id()

        #     # create node-sets for apex right
        #     node_set_kw = create_node_set_keyword(
        #         node_ids=[node_apex_right + 1],
        #         node_set_id=node_set_id_apex_right,
        #         title="apex node right",
        #     )

        # define stimulation node set
        if isinstance(self.model, LeftVentricle):
            stim_nodes = np.array(self.model.left_ventricle.apex_points[0].node_id)

        elif isinstance(self.model, BiVentricle):
            node_apex_left = self.model.left_ventricle.apex_points[0].node_id
            node_apex_right = self.model.right_ventricle.apex_points[0].node_id
            stim_nodes = np.array([node_apex_left, node_apex_right])

        elif isinstance(self.model, (FourChamber, FullHeart)):
            node_apex_left = self.model.left_ventricle.apex_points[0].node_id
            node_apex_right = self.model.right_ventricle.apex_points[0].node_id
            stim_nodes = np.array([node_apex_left, node_apex_right])

            if self.model.right_atrium.get_point("SA_node") != None:
                # Active SA node (belong to both solid and beam)
                stim_nodes = [self.model.right_atrium.get_point("SA_node").node_id]

                #  add more nodes to initiate wave propagation
                # id offset due to cap center nodes TODO do once
                if type(self) == ElectroMechanicsDynaWriter:
                    beam_node_id_offset = len(self.model.cap_centroids)
                else:
                    beam_node_id_offset = 0

                for network in self.model.beam_network:
                    if network.name == "SAN_to_AVN":
                        stim_nodes.append(network.edges[1, 0] + beam_node_id_offset)
                        # stim_nodes.append(network.edges[2, 0] + beam_node_id_offset)
                        # stim_nodes.append(network.edges[3, 0] + beam_node_id_offset)

        # create node-sets for stim nodes
        node_set_id_stimulationnodes = self.get_unique_nodeset_id()
        node_set_kw = create_node_set_keyword(
            node_ids=np.array(stim_nodes) + 1,
            node_set_id=node_set_id_stimulationnodes,
            title="Stim nodes",
        )
        self.kw_database.ep_settings.append(node_set_kw)

        # stimulation
        self.kw_database.ep_settings.append(
            custom_keywords.EmEpTentusscherStimulus(
                stimid=1,
                settype=2,
                setid=node_set_id_stimulationnodes,
                stimstrt=0.0,
                stimt=800.0,
                stimdur=20.0,
                stimamp=50.0,
            )
        )

        # TODO: His bundle is removed for EPMECA model due to unfinished development in LSDYNA
        # instead we directly stimulate the apical points with a delay.
        if type(self) == ElectroMechanicsDynaWriter and isinstance(self.model, FourChamber):
            second_stim_nodes = self.get_unique_nodeset_id()
            beam_node_id_offset = len(self.model.cap_centroids)
            # get apical points
            stim_nodes = [
                apex.node_id
                for p in self.model.parts
                if "ventricle" in p.name
                for apex in p.apex_points
                if "endocardium" in apex.name
            ]

            self.kw_database.ep_settings.append("$$ second stimulation at Left/Right bundle. $$")
            node_set_kw = create_node_set_keyword(
                node_ids=np.array(stim_nodes) + 1,
                node_set_id=second_stim_nodes,
                title="origin of left/right bundle",
            )
            self.kw_database.ep_settings.append(node_set_kw)

            self.kw_database.ep_settings.append(
                custom_keywords.EmEpTentusscherStimulus(
                    stimid=2,
                    settype=2,
                    setid=second_stim_nodes,
                    stimstrt=90.0,  # 90ms delay
                    stimt=800.0,
                    stimdur=20.0,
                    stimamp=50.0,
                )
            )

    def _update_blood_settings(self):
        """Update blood settings."""
        if self.model.info.add_blood_pool == True:
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

    def _update_ECG_coordinates(self):
        """Add ECG computation content."""
        # TODO replace strings by custom dyna keyword
        # TODO handle dynamic numbering of point set ids "psid'
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
                f"{index:>10d} {str(f'{x:9.6f}')[:9]} {str(f'{y:9.6f}')[:9]} {str(f'{z:9.6f}')[:9]}"
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

    def _update_use_Purkinje(self):
        """Update keywords for Purkinje usage."""
        sid = self.get_unique_section_id()
        self.kw_database.beam_networks.append(keywords.SectionBeam(secid=sid, elform=3, a=645))

        if self.__class__.__name__ == "ElectroMechanicsDynaWriter":
            # id offset due to cap center nodes
            beam_node_id_offset = len(self.model.cap_centroids)
            # id offset due to spring type elements
            beam_elem_id_offset = self.id_offset["element"]["discrete"]
        else:
            beam_node_id_offset = 0
            beam_elem_id_offset = 0  # no beam elements introduced before

        # write beam nodes

        # Note: the las beam_network saves all bam nodes
        new_nodes = self.model.beam_network[-1]._all_beam_nodes
        ids = (
            np.linspace(
                len(self.model.mesh.nodes),
                len(self.model.mesh.nodes) + len(new_nodes) - 1,
                len(new_nodes),
                dtype=int,
            )
            + 1  # dyna start by 1
            + beam_node_id_offset  # apply node offset
        )
        nodes_table = np.hstack((ids.reshape(-1, 1), new_nodes))
        kw = add_nodes_to_kw(nodes_table, keywords.Node())
        self.kw_database.beam_networks.append(kw)

        for network in self.model.beam_network:
            ## TODO do not write His Bundle when coupling, it leads to crash
            if self.__class__.__name__ == "ElectroMechanicsDynaWriter" and network.name == "His":
                continue

            # It is previously defined from purkinje generation step
            # but needs to reassign part ID here
            # to make sure no conflict with 4C/full heart case.
            network.pid = self.get_unique_part_id()

            if network.name == "Left-purkinje":
                network.nsid = self.model.left_ventricle.endocardium.id
            elif network.name == "Right-purkinje":
                network.nsid = self.model.right_ventricle.endocardium.id
            elif network.name == "SAN_to_AVN":
                network.nsid = self.model.right_atrium.endocardium.id
            elif network.name == "Left bundle branch":
                network.nsid = self.model.left_ventricle.cavity.surface.id
            elif network.name == "Right bundle branch":
                network.nsid = self.model.right_ventricle.cavity.surface.id
            # His bundle is inside of surface, no segment will associated
            elif network.name == "His":
                network.nsid = -1
            else:
                LOGGER.error(f"Unknown network name for {network.name}.")
                exit()

            # write
            self.kw_database.beam_networks.append(f"$$ {network.name} $$")

            origin_coordinates = network.nodes[network.edges[0, 0]]
            self.kw_database.beam_networks.append(
                custom_keywords.EmEpPurkinjeNetwork2(
                    purkid=network.pid,
                    buildnet=0,
                    ssid=network.nsid,
                    mid=network.pid,
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
                    pmjres=self.settings.electrophysiology.material.beam["pmjres"].m,
                )
            )

            part_df = pd.DataFrame(
                {
                    "heading": [network.name],
                    "pid": [network.pid],
                    "secid": [sid],
                    "mid": [network.pid],
                }
            )
            part_kw = keywords.Part()
            part_kw.parts = part_df
            self.kw_database.beam_networks.append(part_kw)
            self.kw_database.beam_networks.append(keywords.MatNull(mid=network.pid, ro=1e-11))
            beam_material = self.settings.electrophysiology.material.beam
            self.kw_database.beam_networks.append(
                custom_keywords.EmMat001(
                    mid=network.pid,
                    mtype=2,
                    sigma=beam_material["sigma"].m,
                    beta=beam_material["beta"].m,
                    cm=beam_material["cm"].m,
                )
            )

            # cell model
            # use endo property
            self._add_Tentusscher_keyword_endo(matid=network.pid)

            # mesh
            # apply offset for beam connectivity
            connect = network.edges
            connect[network.beam_nodes_mask] += beam_node_id_offset

            beams_kw = keywords.ElementBeam()
            beams_kw = add_beams_to_kw(
                beams=connect + 1,
                beam_kw=beams_kw,
                pid=network.pid,
                offset=beam_elem_id_offset,
            )
            beam_elem_id_offset += len(network.edges)
            self.kw_database.beam_networks.append(beams_kw)

    def _update_export_controls(self):
        """Add solution controls to the main simulation."""
        self.kw_database.main.append(
            keywords.DatabaseBinaryD3Plot(dt=self.settings.electrophysiology.analysis.dt_d3plot.m)
        )

        return


class ElectroMechanicsDynaWriter(MechanicsDynaWriter, ElectrophysiologyDynaWriter):
    """Class for preparing the input for LS-DYNA electromechanical simulation."""

    def __init__(
        self,
        model: HeartModel,
        settings: SimulationSettings = None,
    ) -> None:
        if isinstance(model, FourChamber):
            model._create_isolation_part()

        BaseDynaWriter.__init__(self, model=model, settings=settings)

        self.kw_database = ElectroMechanicsDecks()
        """Collection of keyword decks relevant for mechanics."""

        self.system_model_name = self.settings.mechanics.system.name
        """Name of system model to use."""

        # Depending on the system model specified give list of parameters
        self.cap_in_zerop = True
        """
        If include cap (shell) elements in ZeroPressure.
        Experimental feature, please do not change it.
        """

    def __duplicate_ventricle_atrial_nodes_tie(self):
        """Test feature, not working with mpp < DEV104400."""
        # find interface nodes between ventricles and atrial
        v_ele = np.hstack(
            (self.model.left_ventricle.element_ids, self.model.right_ventricle.element_ids)
        )
        a_ele = np.hstack((self.model.left_atrium.element_ids, self.model.right_atrium.element_ids))

        ventricles = self.model.mesh.extract_cells(v_ele)
        atrial = self.model.mesh.extract_cells(a_ele)

        interface_nids = np.intersect1d(
            ventricles["vtkOriginalPointIds"], atrial["vtkOriginalPointIds"]
        )

        # duplicate these nodes and update mesh
        sets = []
        new_coords = np.array([])
        tets: np.ndarray = self.model.mesh.tetrahedrons
        nid_offset = len(self.model.mesh.nodes)
        cnt = 0

        for old_nid in interface_nids:
            set = np.array([old_nid])
            ele_ids = np.where(np.any(np.isin(tets, old_nid), axis=1))[
                0
            ]  # Ids of elements on this interface node
            for ele_id in ele_ids[1:]:
                # first element will keep the original node,
                # other elements will take duplicated nodes (update mesh)
                new_id = nid_offset + cnt
                new_coords = np.append(new_coords, self.model.mesh.nodes[old_nid, :])
                tets[ele_id][tets[ele_id] == old_nid] = new_id
                set = np.append(set, new_id)
                cnt += 1
            sets.append(set)

        self.model.mesh.nodes = np.vstack((self.model.mesh.nodes, new_coords.reshape(-1, 3)))
        self.model.mesh.tetrahedrons = tets

        # tie duplicated nodes
        sid_offset = self.get_unique_nodeset_id()  # slow and move out of loop
        sid_offset = 100
        count = 0
        for set in sets:
            sid = sid_offset + count
            count += 1
            kw = create_node_set_keyword(set + 1, node_set_id=sid, title="tied_" + str(set[0] + 1))
            self.kw_database.duplicate_nodes.append(kw)
            kw = keywords.ConstrainedTiedNodesFailure(nsid=sid, eppf=1.0e25, etype=1)
            self.kw_database.duplicate_nodes.append(kw)

        # self.kw_database.main.append(keywords.Include(filename="duplicate_nodes.k"))
        return

    def update(self, with_dynain=False):
        """Update the keyword database."""
        if isinstance(self.model, FourChamber):
            self.model.left_atrium.has_fiber = True
            self.model.left_atrium.is_active = True
            self.model.right_atrium.has_fiber = True
            self.model.right_atrium.is_active = True

        MechanicsDynaWriter.update(self, with_dynain=with_dynain)

        if self.model.beam_network:
            # Coupling enabled, EP beam nodes follow the motion of surfaces
            self.kw_database.ep_settings.append(keywords.EmControlCoupling(smcoupl=0))
            self._update_use_Purkinje()
            self.kw_database.main.append(keywords.Include(filename="beam_networks.k"))

        self._update_cellmodels()
        self.kw_database.main.append(keywords.Include(filename="cell_models.k"))

        self._update_ep_settings()
        self._update_stimulation()

        # coupling parameters
        coupling_str = (
            "*EM_CONTROL_TIMESTEP\n"
            "$   TSTYPE   DTCONST      LCID    FACTOR     DTMIN     DTMAX\n"
            "         1       1.0\n"
            "*EM_CONTROL_COUPLING\n"
            "$    THCPL     SMCPL    THLCID    SMLCID\n"
            "                   0\n"
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


class UHCWriter(BaseDynaWriter):
    """Universal Heart Coordinate Writer."""

    def __init__(self, model, type: Literal["uvc", "la_fiber", "ra_fiber"], **kwargs):
        """
        Write thermal input to set up a Laplace dirichlet problem.

        Parameters
        ----------
        model: Heart Model
            Heart model to simulate.
        type : Literal[]
            Type of simulation to set up.
        """
        super().__init__(model=model)
        self.type = type

        # remove unnecessary parts
        if self.type == "uvc":
            parts_to_keep = ["Left ventricle", "Right ventricle", "Septum"]
            self._keep_parts(parts_to_keep)
        elif self.type == "la_fiber":
            parts_to_keep = ["Left atrium"]
            #  A manual point for LA fiber
            for key, value in kwargs.items():
                if key == "laa":
                    self.left_appendage_apex = value
        elif self.type == "ra_fiber":
            parts_to_keep = ["Right atrium"]
            #  A manual point for RA fiber
            for key, value in kwargs.items():
                if key == "raa":
                    self.right_appendage_apex = value

        # remove unnecessary mesh
        if self.type == "uvc":
            elems_to_keep = []
            if isinstance(self.model, LeftVentricle):
                elems_to_keep.extend(model.parts[0].element_ids)
            else:
                elems_to_keep.extend(model.parts[0].element_ids)
                elems_to_keep.extend(model.parts[1].element_ids)
                elems_to_keep.extend(model.parts[2].element_ids)

            model.mesh.clear_data()
            model.mesh["cell_ids"] = np.arange(0, model.mesh.n_cells, dtype=int)
            model.mesh["point_ids"] = np.arange(0, model.mesh.n_points, dtype=int)

            self.target = model.mesh.extract_cells(elems_to_keep)

        elif self.type == "la_fiber" or self.type == "ra_fiber":
            # In original model, mitral/tricuspid valves are assigned with ventricle parts
            # so we need to update caps information at first
            for part in model.parts:
                part.caps = []
                for surface in part.surfaces:
                    surface.edge_groups = []
            model.cap_centroids = []
            model._assign_surfaces_to_parts()
            model._assign_caps_to_parts(unique_mitral_tricuspid_valve=False)

            self._keep_parts(parts_to_keep)
            model.mesh.clear_data()
            model.mesh["cell_ids"] = np.arange(0, model.mesh.n_cells, dtype=int)
            model.mesh["point_ids"] = np.arange(0, model.mesh.n_points, dtype=int)

            self.target = model.mesh.extract_cells(model.parts[0].element_ids)

        # if self.type == "la_fiber":
        #     if 6 != len(self.model.parts[0].caps):
        #         LOGGER.error("Input left atrium is not suitable for set up BC.")
        #         exit(-1)
        # elif self.type == "ra_fiber":
        #     if 3 != len(self.model.parts[0].caps):
        #         LOGGER.error("Input left atrium is not suitable for set up BC.")
        #         exit(-1)

    def additional_right_atrium_bc(self, atrium: pv.UnstructuredGrid):
        """
        Find additional node sets for right atrium.

        Find appendage, top, tricuspid wall and septum node set.

        Parameters
        ----------
        atrium : pv.UnstructuredGrid
            right atrium pyvista object
        """
        # Find appendage apex
        import scipy.spatial as spatial

        tree = spatial.cKDTree(atrium.points)
        # radius = 1.5 mm
        raa_ids = np.array(tree.query_ball_point(self.right_appendage_apex, 1.5))
        if len(raa_ids) == 0:
            LOGGER.error("No node is identified as right atrium appendage apex.")
            exit()

        kw = create_node_set_keyword(raa_ids + 1, node_set_id=11, title="raa")
        self.kw_database.node_sets.append(kw)
        atrium["raa"] = np.zeros(atrium.n_points)
        atrium["raa"][raa_ids] = 1

        # Find top
        for cap in self.model.parts[0].caps:
            if "tricuspid" in cap.name:
                tv_center = cap.centroid
            elif "superior" in cap.name:
                svc_center = cap.centroid
            elif "inferior" in cap.name:
                ivc_center = cap.centroid
        cut_center = np.vstack((tv_center, svc_center, ivc_center)).mean(axis=0)
        cut_normal = np.cross(svc_center - tv_center, ivc_center - tv_center)

        atrium["cell_ids_tmp"] = np.arange(0, atrium.n_cells, dtype=int)
        atrium["point_ids_tmp"] = np.arange(0, atrium.n_points, dtype=int)
        slice = atrium.slice(origin=cut_center, normal=cut_normal)
        crinkled = atrium.extract_cells(np.unique(slice["cell_ids_tmp"]))

        # After cut, select the top region
        x = crinkled.connectivity()
        if np.max(x.point_data["RegionId"]) != 2:
            # Should only have 3 parts
            LOGGER.error("Cannot find top node set for right atrium.")
            exit()

        # compare closest point with TV nodes, top region should be far with TV node set
        tv_tree = spatial.cKDTree(atrium.points[atrium.point_data["tricuspid-valve"] == 1])
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

        # assign
        kw = create_node_set_keyword(top_ids + 1, node_set_id=10, title="top")
        self.kw_database.node_sets.append(kw)
        atrium["top"] = np.zeros(atrium.n_points)
        atrium["top"][top_ids] = 1

        # Find tricuspid_wall and tricuspid_septum
        # need a copied object to do clip, atrium will be corrupted otherwise
        septum, free_wall = copy.deepcopy(atrium).clip(
            origin=cut_center, normal=cut_normal, crinkle=True, return_clipped=True
        )
        # ids in full mesh
        tv_s_ids = septum["point_ids"][np.where(septum["tricuspid-valve"] == 1)]

        tv_s_ids_sub = np.where(np.isin(atrium["point_ids"], tv_s_ids))[0]
        atrium["tv_s"] = np.zeros(atrium.n_points)
        atrium["tv_s"][tv_s_ids_sub] = 1

        kw = create_node_set_keyword(tv_s_ids_sub + 1, node_set_id=12, title="tv_septum")
        self.kw_database.node_sets.append(kw)

        tv_w_ids = free_wall["point_ids"][np.where(free_wall["tricuspid-valve"] == 1)]
        tv_w_ids_sub = np.where(np.isin(atrium["point_ids"], tv_w_ids))[0]
        # remove re constraint nodes
        tv_w_ids_sub = np.setdiff1d(tv_w_ids_sub, tv_s_ids_sub)

        atrium["tv_w"] = np.zeros(atrium.n_points)
        atrium["tv_w"][tv_w_ids_sub] = 1

        kw = create_node_set_keyword(tv_w_ids_sub + 1, node_set_id=13, title="tv_wall")
        self.kw_database.node_sets.append(kw)

    def update_atrium_fiber_bc(self, atrium: pv.UnstructuredGrid):
        """Define boundary condition."""

        def get_nodeset_id_by_cap_name(cap):
            # ID map:
            # RIP 1 LAP 2 RSP 3 MV 4 LIP 5 LSP 6 TV 7 SVC 8 IVC 9
            if "right" in cap.name:
                if "inferior" in cap.name:
                    set_id = 1
                elif "superior" in cap.name:
                    set_id = 3
            elif "left" in cap.name:
                if "appendage" in cap.name:
                    set_id = 2
                elif "inferior" in cap.name:
                    set_id = 5
                elif "superior" in cap.name:
                    set_id = 6
            elif "mitral" in cap.name:
                set_id = 4
            elif "tricuspid" in cap.name:
                set_id = 7
            elif "vena" in cap.name:
                if "superior" in cap.name:
                    set_id = 8
                elif "inferior" in cap.name:
                    set_id = 9
            else:
                set_id = 99
            return set_id

        id_sorter = np.argsort(atrium["point_ids"])
        ids_edges = []  # all nodes belong to valves
        for cap in self.model.parts[0].caps:
            # get node IDs for atrium mesh
            ids_sub = np.where(np.isin(atrium["point_ids"], cap.node_ids))[0]
            # create node set
            set_id = get_nodeset_id_by_cap_name(cap)
            kw = create_node_set_keyword(ids_sub + 1, node_set_id=set_id, title=cap.name)
            self.kw_database.node_sets.append(kw)

            ids_edges.extend(ids_sub)

            # Add info to pyvista object (RA fiber use this)
            atrium[cap.name] = np.zeros(atrium.n_points, dtype=int)
            atrium[cap.name][ids_sub] = 1

        if self.type == "la_fiber" and hasattr(self, "left_appendage_apex"):
            import scipy.spatial as spatial

            tree = spatial.cKDTree(atrium.points)
            # radius = 1.5 mm
            laa_ids = np.array(tree.query_ball_point(self.left_appendage_apex, 1.5))
            kw = create_node_set_keyword(laa_ids + 1, node_set_id=2, title="left atrium appendage")
            self.kw_database.node_sets.append(kw)

        # endo nodes ID
        ids_endo = np.where(np.isin(atrium["point_ids"], self.model.parts[0].endocardium.node_ids))[
            0
        ]

        atrium["endo"] = np.zeros(atrium.n_points, dtype=int)
        atrium["endo"][ids_endo] = 1
        kw = create_node_set_keyword(ids_endo + 1, node_set_id=100, title="endo")
        self.kw_database.node_sets.append(kw)

        # epi node ID
        # epi cannot use directly Surface because new free surface exposed
        ids_surface = atrium.extract_surface()["vtkOriginalPointIds"]
        ids_epi = np.setdiff1d(ids_surface, ids_endo)
        ids_epi = np.setdiff1d(ids_epi, ids_edges)

        atrium["epi"] = np.zeros(atrium.n_points, dtype=int)
        atrium["epi"][ids_epi] = 1
        kw = create_node_set_keyword(ids_epi + 1, node_set_id=200, title="epi")
        self.kw_database.node_sets.append(kw)

        # set BC in DYNA case
        if self.type == "la_fiber":
            self.kw_database.main.append(keywords.Case(caseid=1, jobid="trans", scid1=1))
            self.kw_database.main.append(keywords.Case(caseid=2, jobid="ab", scid1=2))
            self.kw_database.main.append(keywords.Case(caseid=3, jobid="v", scid1=3))
            self.kw_database.main.append(keywords.Case(caseid=4, jobid="r", scid1=4))

            self.kw_database.main.append("*CASE_BEGIN_1")
            self._define_Laplace_Dirichlet_bc(set_ids=[100, 200], bc_values=[0, 1])
            self.kw_database.main.append("*CASE_END_1")

            self.kw_database.main.append("*CASE_BEGIN_2")
            self._define_Laplace_Dirichlet_bc(
                set_ids=[1, 3, 4, 5, 6, 2], bc_values=[2.0, 2.0, 1.0, 0.0, 0.0, -1.0]
            )
            self.kw_database.main.append("*CASE_END_2")

            self.kw_database.main.append("*CASE_BEGIN_3")
            self._define_Laplace_Dirichlet_bc(set_ids=[1, 3, 5, 6], bc_values=[1.0, 1.0, 0.0, 0.0])
            self.kw_database.main.append("*CASE_END_3")

            self.kw_database.main.append("*CASE_BEGIN_4")
            self._define_Laplace_Dirichlet_bc(
                set_ids=[4, 1, 2, 3, 5, 6], bc_values=[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            )
            self.kw_database.main.append("*CASE_END_4")

        elif self.type == "ra_fiber":
            self.additional_right_atrium_bc(atrium)

            self.kw_database.main.append(keywords.Case(caseid=1, jobid="trans", scid1=1))
            self.kw_database.main.append(keywords.Case(caseid=2, jobid="ab", scid1=2))
            self.kw_database.main.append(keywords.Case(caseid=3, jobid="v", scid1=3))
            self.kw_database.main.append(keywords.Case(caseid=4, jobid="r", scid1=4))
            self.kw_database.main.append(keywords.Case(caseid=5, jobid="w", scid1=5))

            self.kw_database.main.append("*CASE_BEGIN_1")
            self._define_Laplace_Dirichlet_bc(set_ids=[100, 200], bc_values=[0, 1])
            self.kw_database.main.append("*CASE_END_1")

            self.kw_database.main.append("*CASE_BEGIN_2")
            self._define_Laplace_Dirichlet_bc(
                set_ids=[9, 7, 8, 11], bc_values=[2.0, 1.0, 0.0, -1.0]
            )
            self.kw_database.main.append("*CASE_END_2")

            self.kw_database.main.append("*CASE_BEGIN_3")
            self._define_Laplace_Dirichlet_bc(set_ids=[9, 8, 11], bc_values=[1.0, 0.0, 0.0])
            self.kw_database.main.append("*CASE_END_3")

            self.kw_database.main.append("*CASE_BEGIN_4")
            self._define_Laplace_Dirichlet_bc(set_ids=[7, 10], bc_values=[1.0, 0.0])
            self.kw_database.main.append("*CASE_END_4")

            # Differently with article, we add Gamma_top = 0 to enforce BC
            self.kw_database.main.append("*CASE_BEGIN_5")
            self._define_Laplace_Dirichlet_bc(set_ids=[12, 13, 10], bc_values=[1.0, -1.0, 0.0])
            self.kw_database.main.append("*CASE_END_5")

        return atrium

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
        elif self.type == "la_fiber" or self.type == "ra_fiber":
            self.update_atrium_fiber_bc(self.target)

        self._get_list_of_includes()
        self._add_includes()

    def _update_uvc_bc(self):
        self.kw_database.main.append(keywords.Case(caseid=1, jobid="transmural", scid1=1))
        self.kw_database.main.append(keywords.Case(caseid=2, jobid="apico-basal", scid1=2))
        self.kw_database.main.append(keywords.Case(caseid=3, jobid="rotational", scid1=3))

        # transmural uvc
        endo_set = []
        epi_set = []
        for part in self.model.parts:
            for surf in part.surfaces:
                if "endocardium" in surf.name:
                    endo_set.extend(surf.node_ids)
                # elif "epicardium" in surf.name:
                #     epi_set.extend(surf.node_ids)

        # map IDs to sub mesh
        endo_set_new = np.where(np.isin(self.target["point_ids"], endo_set))[0]

        endo_sid = self.get_unique_nodeset_id()
        kw = create_node_set_keyword(endo_set_new + 1, node_set_id=endo_sid, title="endo")
        self.kw_database.node_sets.append(kw)

        # epi_set_new = id_sorter[
        #     np.searchsorted(self.target["point_ids"], epi_set, sorter=id_sorter)
        # ]
        # epi_set_new = np.unique(epi_set_new)
        # epi_set_new = np.setdiff1d(epi_set_new, endo_set_new)

        # epi cannot use directly Surface because new free surface exposed
        ids_surface = self.target.extract_surface()["vtkOriginalPointIds"]
        epi_set_new = np.setdiff1d(ids_surface, endo_set_new)

        epi_sid = self.get_unique_nodeset_id()
        kw = create_node_set_keyword(epi_set_new + 1, node_set_id=epi_sid, title="epi")
        self.kw_database.node_sets.append(kw)

        self.kw_database.main.append("*CASE_BEGIN_1")
        self._define_Laplace_Dirichlet_bc(set_ids=[endo_sid, epi_sid], bc_values=[0, 1])
        self.kw_database.main.append("*CASE_END_1")

        # apicobasal uvc
        apex_sid = self._create_apex_nodeset()
        base_sid = self._create_base_nodeset()

        self.kw_database.main.append("*CASE_BEGIN_2")
        self._define_Laplace_Dirichlet_bc(set_ids=[apex_sid, base_sid], bc_values=[0, 1])
        self.kw_database.main.append("*CASE_END_2")

        # rotational uc
        [sid_minus_pi, sid_plus_pi, sid_zero] = self._create_rotational_nodesets()

        self.kw_database.main.append("*CASE_BEGIN_3")
        self._define_Laplace_Dirichlet_bc(
            set_ids=[sid_minus_pi, sid_plus_pi, sid_zero], bc_values=[-np.pi, np.pi, 0]
        )
        self.kw_database.main.append("*CASE_END_3")

    def _create_apex_nodeset(self):
        # apex
        # select a region within 1 cm, this seems consistent with Strocchi database
        apex_set = self.model.get_apex_node_set(radius=10)
        # get local ID
        ids_submesh = np.where(np.isin(self.target["point_ids"], apex_set))[0]

        sid = self.get_unique_nodeset_id()
        kw = create_node_set_keyword(ids_submesh + 1, node_set_id=sid, title="apex")
        self.kw_database.node_sets.append(kw)
        return sid

    def _create_base_nodeset(self):
        # base
        base_set = np.array([])
        for part in self.model.parts:
            for cap in part.caps:
                #  Strocchi database use only mv and tv
                # if ("mitral" in cap.name) or ("tricuspid" in cap.name):
                base_set = np.append(base_set, cap.node_ids)
        # get local ID
        ids_submesh = np.where(np.isin(self.target["point_ids"], base_set))[0]

        sid = self.get_unique_nodeset_id()
        kw = create_node_set_keyword(ids_submesh + 1, node_set_id=sid, title="base")
        self.kw_database.node_sets.append(kw)
        return sid

    def _create_surface_nodeset(self, surftype: str, cavity_type: str):
        nodeset = np.array([])
        for part in self.model.parts:
            if cavity_type in part.name:
                for surf in part.surfaces:
                    if surftype in surf.name:
                        nodeset = np.append(nodeset, surf.node_ids)
        nodeset = np.unique(nodeset.astype(int))

        # map IDs to sub mesh
        ids_submesh = np.where(np.isin(self.target["point_ids"], nodeset))[0]

        sid = self.get_unique_nodeset_id()
        kw = create_node_set_keyword(
            ids_submesh + 1, node_set_id=sid, title=cavity_type + " " + surftype + " all"
        )
        self.kw_database.node_sets.append(kw)

        return sid

    def _create_rotational_nodesets(self):
        # Find nodes on target mesh
        rot_start, rot_end, septum = self.model._compute_uvc_rotation_bc(copy.deepcopy(self.target))

        sid_minus_pi = self.get_unique_nodeset_id()
        kw = create_node_set_keyword(rot_start + 1, node_set_id=sid_minus_pi, title="rotation:-pi")
        self.kw_database.node_sets.append(kw)
        sid_plus_pi = self.get_unique_nodeset_id()
        kw = create_node_set_keyword(rot_end + 1, node_set_id=sid_plus_pi, title="rotation:pi")
        self.kw_database.node_sets.append(kw)
        sid_zero = self.get_unique_nodeset_id()
        kw = create_node_set_keyword(septum + 1, node_set_id=sid_zero, title="rotation:0")
        self.kw_database.node_sets.append(kw)
        return [sid_minus_pi, sid_plus_pi, sid_zero]

    def _define_Laplace_Dirichlet_bc(
        self,
        set_ids: List[int],
        bc_values: List[float],
    ):
        for sid, value in zip(set_ids, bc_values):
            self.kw_database.main.append(
                keywords.BoundaryTemperatureSet(
                    nsid=sid,
                    lcid=0,
                    cmult=value,
                ),
            )
        return

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
        self.kw_database.main.append(keywords.ControlTermination(endtim=1, dtmin=1.0))


if __name__ == "__main__":
    print("protected")
