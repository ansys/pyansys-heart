"""Module contain. classes for writing LS-DYNA keywords based.

Note
----
Uses a HeartModel (from ansys.heart.preprocessor.models).

"""
import copy
import json
import os
import time
from typing import List, Literal

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
from ansys.heart.simulator.settings.settings import SimulationSettings

# import missing keywords
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
    MaterialHGOMyocardium,
    MaterialNeoHook,
    active_curve,
)
from ansys.heart.writer.system_models import _ed_load_template, define_function_windkessel
import numpy as np
import pandas as pd
import pkg_resources
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

    def _update_node_db(self, ids=None):
        """Add nodes to the Node database."""
        LOGGER.debug("Updating node keywords...")
        node_kw = keywords.Node()
        if ids is not None:
            nodes = np.vstack([ids, self.model.mesh.nodes.T]).T
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
            part.pid = self.get_unique_part_id()
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
            # use fast element writer for solid ortho elements
            if deckname == "solid_elements":
                element_kws = deck.get_kwds_by_type("ELEMENT")
                if os.path.isfile(filepath):
                    os.remove(filepath)

                for element_kw in element_kws:
                    fast_element_writer(element_kw, filepath)

                fid = open(filepath, "a")
                fid.write("*END")

            # elif deckname == "nodes":
            #     ids = np.arange(0, self.model.mesh.nodes.shape[0], 1) + 1
            #     content = np.hstack((ids.reshape(-1, 1), self.model.mesh.nodes))
            #     np.savetxt(
            #         os.path.join(export_directory, "nodes.k"),
            #         content,
            #         fmt="%8d%16.5e%16.5e%16.5e",
            #         header="*KEYWORD\n*NODE\n"
            #         "$#   nid               x               y               z      tc      rc",
            #         footer="*END",
            #         comments="",
            #     )
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
        Create Solid ortho elements for all cavities.

        Note
        ----
        Each cavity contains one myocardium.

        Parameters
        ----------
        add_fibers: bool, True
            if add fiber information into solid element.
        """
        LOGGER.debug("Updating solid element keywords...")

        if add_fibers:
            cell_data_fields = self.model.mesh.cell_data.keys()
            if "fiber" not in cell_data_fields or "sheet" not in cell_data_fields:
                raise KeyError("Mechanics writer requires fiber and sheet fields")

        # create elements for each part
        for part in self.model.parts:
            # Atrium do not contain fiber information in any way.
            if add_fibers:
                if "ventricle" in part.name.lower() or "septum" in part.name.lower():
                    part_add_fibers = True
                else:
                    part_add_fibers = False
            else:
                part_add_fibers = add_fibers

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

        # for boundary conditions
        if isinstance(self.model, (FourChamber, FullHeart)):
            self._add_cap_bc(bc_type="fix_caps")
        else:
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

    def _update_main_db(self, add_damping: bool = True):
        """Update the main .k file.

        Note
        ----
        Consider using a settings (json?) file as input.

        """
        LOGGER.debug("Updating main keywords...")

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

        if add_damping:
            self._add_damping()

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

        # # add auto controls
        # lcid = self.get_unique_curve_id()
        # # tune time step for better compromise between convergence and performance
        # time = [0, prefill_time, prefill_time + dtmax, end_time]
        # step = [5 * dtmax, 5 * dtmax, dtmin, dtmax]
        # kw_curve = create_define_curve_kw(
        #     x=time,
        #     y=step,
        #     curve_name="time step control",
        #     curve_id=lcid,
        #     lcint=0,
        # )
        # self.kw_database.main.append(kw_curve)
        self.kw_database.main.append(
            keywords.ControlImplicitAuto(iauto=1, dtmin=dtmin, dtmax=dtmax)
        )

        # add general implicit controls
        self.kw_database.main.append(
            keywords.ControlImplicitGeneral(imflag=1, dt0=dtmax)
        )  # imflag=1 means implicit

        # add implicit solution controls
        # Nil's suggestion
        self.kw_database.main.append(
            keywords.ControlImplicitSolution(
                dctol=0.01, ectol=1e6, rctol=1e3, abstol=-1e-20, dnorm=1, nlnorm=4, lsmtd=5
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

        self.kw_database.main.append(keywords.DatabaseExtentBinary(neiph=27, strflg=1, maxint=0))

        # control ELOUT file to extract left ventricle's stress/strain
        if hasattr(self.model, "septum"):
            self.kw_database.main.append(
                keywords.SetSolidGeneral(
                    option="PART",
                    sid=1,
                    e1=self.model.left_ventricle.pid,
                    e2=self.model.septum.pid,
                    user_comment="create left ventricle + septum set for exporting",
                )
            )
        else:
            self.kw_database.main.append(
                keywords.SetSolidGeneral(option="PART", sid=1, e1=self.model.left_ventricle.pid)
            )
        self.kw_database.main.append(keywords.DatabaseHistorySolidSet(id1=1))

        # lcid = self.get_unique_curve_id()
        # time = [
        #     0,
        #     self.parameters["Time"]["End Time"] * 0.8 * 0.99,
        #     self.parameters["Time"]["End Time"] * 0.8,
        #     self.parameters["Time"]["End Time"],
        # ]
        # step = [100 * dt_output_d3plot, 100 * dt_output_d3plot, dt_output_d3plot,
        #         dt_output_d3plot]
        # kw_curve = create_define_curve_kw(
        #     x=time,
        #     y=step,
        #     curve_name="elout control, only save during the last 20% ",
        #     curve_id=lcid,
        #     lcint=0,
        # )
        # self.kw_database.main.append(kw_curve)

        # self.kw_database.main.append(
        #     keywords.DatabaseElout(dt=dt_output_d3plot, binary=2, option1=27)
        # )

        return

    def _add_damping(self):
        """Add damping to the main file."""
        lcid_damp = self.get_unique_curve_id()

        kw_damp = keywords.DampingGlobal(lcid=lcid_damp)

        kw_damp_curve = create_define_curve_kw(
            x=[0, 10e25],  # to create a constant curve
            y=self.settings.mechanics.analysis.global_damping.m * np.array([1, 1]),
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

    def _update_material_db(self, add_active: bool = True):
        """Update the database of material keywords."""
        act_curve_id = self.get_unique_curve_id()

        material_settings = copy.deepcopy(self.settings.mechanics.material)
        # NOTE: since we remove units, we don't have to access quantities by <var_name>.m
        material_settings._remove_units()

        for part in self.model.parts:
            if "ventricle" in part.name.lower() or "septum" in part.name.lower():
                if not add_active:
                    active_dict = None
                else:
                    active_dict = material_settings.myocardium["active"]

                myocardium_kw = MaterialHGOMyocardium(
                    mid=part.mid,
                    iso_user=material_settings.myocardium["isotropic"],
                    anisotropy_user=material_settings.myocardium["anisotropic"],
                    active_user=active_dict,
                )

                myocardium_kw.acid = act_curve_id

                self.kw_database.material.append(myocardium_kw)

            elif "atrium" in part.name:
                # add atrium material
                if material_settings.atrium["type"] == "NeoHook":
                    # use MAT77H
                    atrium_kw = MaterialNeoHook(
                        mid=part.mid,
                        rho=material_settings.atrium["rho"],
                        c10=material_settings.atrium["mu1"] / 2,
                    )
                else:
                    # use MAT295
                    atrium_kw = MaterialHGOMyocardium(
                        mid=part.mid, iso_user=dict(material_settings.atrium)
                    )

                self.kw_database.material.append(atrium_kw)

            else:
                LOGGER.warning("Assuming same material as atrium for: {0}".format(part.name))

                general_tissue_kw = MaterialNeoHook(
                    mid=part.mid,
                    rho=material_settings.atrium["rho"],
                    c10=material_settings.atrium["mu1"] / 2,
                )
                # general_tissue_kw = MaterialHGOMyocardium(
                #     mid=part.mid,
                #     iso_user=dict(material_settings.atrium),
                # )
                self.kw_database.material.append(general_tissue_kw)

        if add_active:
            # write and add active curve to material database
            if material_settings.myocardium["active"]["actype"] == 1:
                time_array, calcium_array = active_curve("constant")
            elif material_settings.myocardium["active"]["actype"] == 2:
                time_array, calcium_array = active_curve("Strocchi2020")

            active_curve_kw = create_define_curve_kw(
                x=time_array,
                y=calcium_array,
                curve_name="calcium_concentration",
                curve_id=act_curve_id,
                lcint=10000,
            )

            # x scaling from beat rate
            active_curve_kw.sfa = 1 / material_settings.myocardium["active"]["heart rate"]
            # y scaling from Ca2
            active_curve_kw.sfo = material_settings.myocardium["active"]["ca2ionm"]

            self.kw_database.material.append(active_curve_kw)

        return

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
            ]
            if isinstance(self, ZeroPressureMechanicsDynaWriter):
                caps_to_use.extend(
                    [
                        "aortic-valve",
                        "pulmonary-valve",
                    ]
                )

        elif isinstance(self.model, (FourChamber, FullHeart)):
            caps_to_use = [
                "superior-vena-cava",
                "right-inferior-pulmonary-vein",
                "right-superior-pulmonary-vein",
            ]
            if isinstance(self, ZeroPressureMechanicsDynaWriter):
                caps_to_use.extend(
                    [
                        "left-superior-pulmonary-vein",
                        "left-inferior-pulmonary-vein",
                        "inferior-vena-cava",
                        "pulmonary-valve",
                    ]
                )

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
                spring_stiffness = bc_settings.valve["biventricle"].m

            elif isinstance(self.model, (FourChamber, FullHeart)):
                spring_stiffness = bc_settings.valve["fourchamber"].m

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

        # use pre-computed nodal area
        nodal_areas = self.model.mesh.point_data["nodal_areas"][attached_nodes]

        # scaling this node due to large deformation
        # nodal_areas[np.where(attached_nodes == cap.node_ids[0])[0][0]] *=len(cap.node_ids)

        # scaled spring stiffness by nodal area
        scale_factor_normal *= nodal_areas
        scale_factor_radial *= nodal_areas

        # add part, section discrete, mat spring, sd_orientiation, element discrete

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
        boundary_conditions = copy.deepcopy(self.settings.mechanics.boundary_conditions)
        boundary_conditions._remove_units()
        pericardium_settings = boundary_conditions.pericardium

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

        penalty_function = (
            -_sigmoid(
                (abs(uvc_l) - pericardium_settings["penalty_function"][0])
                * pericardium_settings["penalty_function"][1]
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
            epicardium_faces = np.vstack([epicardium_faces, surface.triangles])

        # NOTE: some duplicates may exist - fix this in preprocessor
        _, idx, counts = np.unique(epicardium_nodes, return_index=True, return_counts=True)
        if np.any(counts > 1):
            LOGGER.warning("Duplicate nodes found in pericardium")
        epicardium_nodes = epicardium_nodes[np.sort(idx)]

        # select only nodes for penalty factor > 0.0001
        pericardium_nodes = epicardium_nodes[penalty_function[epicardium_nodes] > 0.0001]
        # select surfaces containing these nodes
        pericardium_faces = epicardium_faces[
            np.any(np.isin(epicardium_faces, pericardium_nodes), axis=1)
        ]
        # some nodes on the edge must be included
        pericardium_nodes, a = np.unique(pericardium_faces, return_inverse=True)

        # build pericardium polydata
        coord = self.model.mesh.nodes[pericardium_nodes]
        connect = a.reshape(pericardium_faces.shape)
        pericardium_polydata = vtkmethods.create_vtk_surface_triangles(coord, connect, clean=False)
        # vtkmethods.write_vtkdata_to_vtkfile(pericardium_polydata,'pericardium.vtk')

        # compute normal
        cell_normal, point_normal = vtkmethods.add_normals_to_polydata(
            pericardium_polydata, return_normals=True
        )
        # normal_obj = vtkmethods.add_normals_to_polydata(pericardium_polydata)
        # vtkmethods.write_vtkdata_to_vtkfile(normal_obj,'normal.vtk')

        # use pre-computed nodal areas
        nodal_areas = self.model.mesh.point_data["nodal_areas"][pericardium_nodes]
        nodal_penalty = penalty_function[pericardium_nodes]
        # compute scale factor
        scale_factors = nodal_areas * nodal_penalty

        def __debug():
            import meshio

            meshio.write_points_cells(
                "pericardium.vtk",
                coord,
                [("triangle", connect)],
                point_data={"area": nodal_areas, "normal": point_normal, "penalty": nodal_penalty},
                cell_data={"normal": [cell_normal]},
            )

        # __debug()

        # create unique ids for keywords
        part_id = self.get_unique_part_id()
        section_id = self.get_unique_section_id()
        mat_id = self.get_unique_mat_id()

        # define part
        part_kw = keywords.Part()
        part_kw.parts = pd.DataFrame(
            {"heading": ["Pericardium"], "pid": [part_id], "secid": [section_id], "mid": [mat_id]}
        )
        # define section
        section_kw = keywords.SectionDiscrete(secid=section_id, cdl=0, tdl=0)
        # define material
        mat_kw = keywords.MatSpringElastic(mid=mat_id, k=pericardium_settings["spring_stiffness"])

        # define spring orientations
        sd_orientation_kw = create_define_sd_orientation_kw(
            vectors=point_normal, vector_id_offset=self.id_offset["vector"]
        )
        # add offset
        self.id_offset["vector"] = sd_orientation_kw.vectors["vid"].to_numpy()[-1]
        vector_ids = sd_orientation_kw.vectors["vid"].to_numpy().astype(int)

        # define spring nodes
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
                e=material_settings.cap["mu1"] * 1000,
            )

        section_kw = keywords.SectionShell(
            secid=section_id,
            elform=4,
            shrf=0.8333,
            nip=3,
            t1=material_settings.cap["thickness"],
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
                if add_mesh:
                    # Add center node
                    node_kw = keywords.Node()
                    df = pd.DataFrame(
                        data=np.insert(cap.centroid, 0, cap.centroid_id + 1).reshape(1, -1),
                        columns=node_kw.nodes.columns[0:4],
                    )
                    node_kw.nodes = df
                    # comment the line '*NODE' so nodes.k can be parsed by zerop solver correctly
                    # otherwise, these nodes will not be updated in iterations
                    s = "$" + node_kw.write()
                    self.kw_database.nodes.append(s)

                # center node constraint: average of all edge nodes
                constraint = keywords.ConstrainedInterpolation(
                    icid=len(cap_names_used) + 1,
                    dnid=cap.centroid_id + 1,
                    ddof=123,
                    ityp=1,
                    fgm=1,
                    inid=cap.nsid,
                    idof=123,
                )
                self.kw_database.cap_elements.append(constraint)

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
                file_path = pkg_resources.resource_filename(
                    "ansys.heart.writer", "templates/system_model_settings_bv.json"
                )

            elif isinstance(self.model, LeftVentricle):
                file_path = pkg_resources.resource_filename(
                    "ansys.heart.writer", "templates/system_model_settings_lv.json"
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

        self._update_main_db(add_damping=False)

        self.kw_database.main.title = self.model.info.model_type + " zero-pressure"

        self._update_node_db()
        self._update_parts_db()
        self._update_solid_elements_db(add_fibers=True)
        self._update_segmentsets_db()
        self._update_nodesets_db()
        self._update_material_db(add_active=False)

        # for boundary conditions
        self._add_cap_bc(bc_type="fix_caps")

        if self.cap_in_zerop:
            # define cap element
            self._update_cap_elements_db()

        # # Approximate end-diastolic pressures
        pressure_lv = bc_settings.end_diastolic_cavity_pressure["left_ventricle"].m
        pressure_rv = bc_settings.end_diastolic_cavity_pressure["right_ventricle"].m
        self._add_enddiastolic_pressure_bc(pressure_lv=pressure_lv, pressure_rv=pressure_rv)

        # zerop key words
        self._add_control_reference_configuration()
        #
        # # export dynain file
        # NOTE: generates a new part-set. Use part-set id 999.
        # Please note that choosing 999 as the part-set id is arbitrary,
        # and defining a new part set adding this to the main database will
        # create a part-set id of 999+1
        self.kw_database.main.append(keywords.SetPartListGenerate(sid=999, b1beg=1, b1end=999999))
        self.kw_database.main.append(
            custom_keywords.InterfaceSpringbackLsdyna(
                psid=999, nshv=999, ftype=3, rflag=1, optc="OPTCARD", ndflag=1, cflag=1, hflag=1
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

        return

    def _add_solution_controls(self):
        """Rewrite method for the zerop simulation."""
        settings = copy.deepcopy(self.settings.stress_free)
        settings._remove_units()

        self.kw_database.main.append(keywords.ControlTermination(endtim=settings.analysis.end_time))

        self.kw_database.main.append(keywords.ControlImplicitDynamics(imass=0))
        # self.kw_database.main.append(
        #     keywords.ControlImplicitDynamics(imass=1, gamma=0.6, beta=0.38)
        # )

        # add auto controls
        self.kw_database.main.append(
            keywords.ControlImplicitAuto(
                iauto=1, dtmin=settings.analysis.dtmin, dtmax=settings.analysis.dtmax
            )
        )

        # add general implicit controls
        self.kw_database.main.append(
            keywords.ControlImplicitGeneral(imflag=1, dt0=settings.analysis.dtmax)
        )

        # add implicit solution controls: Defaults are OK?
        self.kw_database.main.append(keywords.ControlImplicitSolution())

        # add implicit solver controls
        self.kw_database.main.append(custom_keywords.ControlImplicitSolver())

        # add binout for post-process
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

            # remove nodes not attached to ventricle parts
            self.model.mesh.nodes = self.model.mesh.nodes[nids]
            self._update_node_db(ids=nids + 1)

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
        self._update_nodesets_db()

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

        for part in parts:
            element_ids = part.element_ids
            em_mat_id = self.get_unique_mat_id()
            self.kw_database.material.extend(
                [
                    keywords.MatElastic(mid=em_mat_id, ro=1e-6, e=1),
                    custom_keywords.EmMat003(
                        mid=em_mat_id,
                        mtype=2,
                        sigma11=0.5,
                        sigma22=0.1,
                        sigma33=0.1,
                        beta=140,
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
        for part in self.model.parts:
            em_mat_id = part.pid
            self.kw_database.material.extend(
                [
                    keywords.MatElastic(mid=em_mat_id, ro=1e-6, e=1),
                    custom_keywords.EmMat003(
                        mid=em_mat_id,
                        mtype=2,
                        sigma11=0.5,
                        sigma22=0.1,
                        sigma33=0.1,
                        beta=140,
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
        super().__init__(model=model, settings=settings)

        self.kw_database = ElectrophysiologyDecks()
        """Collection of keywords relevant for Electrophysiology."""

    def update(self):
        """Update keyword database for Electrophysiology."""
        self._isolate_atria_and_ventricles()
        ##
        self._update_main_db()

        self._update_solution_controls()
        self._update_export_controls()
        self._update_node_db()

        self._update_parts_db()
        self._update_solid_elements_db(add_fibers=True)
        self._update_material_db()
        self._update_ep_material_db()
        self._update_cellmodels()
        self._update_segmentsets_db()
        self._update_nodesets_db()
        self._update_use_Purkinje()
        # update ep settings
        self._update_ep_settings()

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

                    tets: np.ndarray = self.model.mesh.tetrahedrons
                    tets[self.model.left_atrium.element_ids, :] = tets_atrium

                    self.model.mesh.tetrahedrons = tets

                    self.model.mesh.nodes = np.append(
                        self.model.mesh.nodes, self.model.mesh.nodes[interface_nids, :], axis=0
                    )

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
        """Add simple mechanics material for each defined part."""
        for part in self.model.parts:
            ep_mid = part.pid
            self.kw_database.material.extend(
                [
                    keywords.MatElastic(mid=ep_mid, ro=1e-6, e=1),
                ]
            )

    def _update_ep_material_db(self):
        """Add EP material for each defined part."""
        for part in self.model.parts:
            partname = part.name.lower()
            if ("atrium" in partname) or ("ventricle" in partname) or ("septum" in partname):
                ep_mid = part.pid
                self.kw_database.material.extend(
                    [
                        custom_keywords.EmMat003(
                            mid=ep_mid,
                            mtype=2,
                            sigma11=0.5,
                            sigma22=0.1,
                            sigma33=0.1,
                            beta=140,
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
            else:
                ep_mid = part.pid
                self.kw_database.material.extend(
                    [
                        keywords.EmMat001(mid=ep_mid, mtype=4, sigma=1),
                    ]
                )

    def _update_cellmodels(self):
        """Add cell model for each defined part."""
        for part in self.model.parts:
            partname = part.name.lower()
            if ("atrium" in partname) or ("ventricle" in partname) or ("septum" in partname):
                ep_mid = part.pid
                cell_kw = keywords.EmEpCellmodelTentusscher(
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
                )
                cell_kw.gas_constant = 8314.472
                cell_kw.faraday_constant = 96485.3415
                self.kw_database.cell_models.extend([cell_kw])

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
        part_ids = [None] * 7
        part_ids[0 : len(self.model.part_ids)] = self.model.part_ids
        nsid_all_parts = self.get_unique_nodeset_id()
        kw = keywords.SetNodeGeneral(
            title="All nodes",
            option="PART",
            sid=nsid_all_parts,
            e1=part_ids[0],
            e2=part_ids[1],
            e3=part_ids[2],
            e4=part_ids[3],
            e5=part_ids[4],
            e6=part_ids[5],
            e7=part_ids[6],
        )
        self.kw_database.node_sets.append(kw)

        # use defaults
        self.kw_database.ep_settings.append(custom_keywords.EmControlEp(numsplit=1))

        # max iter should be int
        self.kw_database.ep_settings.append(
            keywords.EmSolverFem(reltol=1e-6, maxite=int(1e4), precon=2)
        )

        self.kw_database.ep_settings.append(keywords.EmOutput(mats=1, matf=1, sols=1, solf=1))

        if isinstance(self.model, BiVentricle):
            node_apex_left = self.model.left_ventricle.apex_points[0].node_id
            node_apex_right = self.model.right_ventricle.apex_points[0].node_id

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
            # if isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
            #     self.model.left_atrium.apex_points
        if isinstance(self.model, (FourChamber, FullHeart)):
            node_apex_left = self.model.left_ventricle.apex_points[0].node_id
            node_apex_right = self.model.right_ventricle.apex_points[0].node_id
            stim_nodes = np.array([node_apex_left, node_apex_right])

            if self.model.right_atrium.get_point("SA_node") != None:
                stim_nodes = self.model.right_atrium.get_point("SA_node").node_id
            # TODO add more nodes to initiate wave propagation !!!!
            node_set_id_stimulationnodes = self.get_unique_nodeset_id()
            # create node-sets for apex
            node_set_kw = create_node_set_keyword(
                node_ids=stim_nodes + 1,
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
            # if isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
            #     self.model.left_atrium.apex_points

        elif isinstance(self.model, LeftVentricle):
            node_apex_left = self.model.left_ventricle.apex_points[0].node_id

            node_set_id_apex_left = self.get_unique_nodeset_id()
            # create node-sets for apex left
            node_set_kw = create_node_set_keyword(
                node_ids=[node_apex_left + 1],
                node_set_id=node_set_id_apex_left,
                title="apex node left",
            )
            self.kw_database.node_sets.append(node_set_kw)

            # TODO add more nodes to initiate wave propagation !!!!
            node_set_id_stimulationnodes = self.get_unique_nodeset_id()
            # create node-sets for apex
            node_set_kw = create_node_set_keyword(
                node_ids=[node_apex_left + 1],
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

    def _update_solution_controls(
        self,
        end_time: float = 800,
    ):
        """Add solution controls and other solver settings as keywords."""
        # add termination keywords
        self.kw_database.main.append(keywords.ControlTermination(endtim=end_time, dtmin=0.0))

        self.kw_database.main.append(keywords.ControlTimeStep(dtinit=1.0, dt2ms=1.0))
        return

    def _update_main_db(self):
        return

    def _get_list_of_includes(self):
        """Get a list of files to include in main.k. omit any empty decks."""
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

    def _update_use_Purkinje(self):
        """Update keywords for Purkinje usage."""
        if self.model.mesh.beam_network:
            self.kw_database.parts.append(keywords.SectionBeam(secid=3, elform=3, a=645))
            if self.__class__.__name__ == "ElectroMechanicsDynaWriter":
                self.kw_database.ep_settings.append(keywords.EmControlCoupling(smcoupl=0))
            else:
                self.kw_database.ep_settings.append(keywords.EmControlCoupling(smcoupl=1))
            beams_kw = keywords.ElementBeam()
            for network in self.model.mesh.beam_network:
                # It is previously defined from purkinje generation step
                # but needs to reassign part ID here
                # to make sure no conflict with 4C/full heart case.
                network.pid = self.get_unique_part_id()

                origin_coordinates = self.model.mesh.nodes[network.node_ids[0], :]
                if network.name == None:
                    node_apex_left = self.model.left_ventricle.apex_points[0].xyz
                    node_apex_right = self.model.right_ventricle.apex_points[0].xyz
                    distance = np.linalg.norm(
                        origin_coordinates - np.array([node_apex_left, node_apex_right]),
                        axis=1,
                    )
                    if np.min(distance[0]) < 1e-3:
                        network.name = "Left" + "-" + "purkinje"
                        network.nsid = self.model.left_ventricle.endocardium.id
                    elif np.min(distance[1]) < 1e-3:
                        network.name = "Right" + "-" + "purkinje"
                        network.nsid = self.model.right_ventricle.endocardium.id
                    else:
                        LOGGER.error("Point too far from apex")

                self.kw_database.main.append(
                    custom_keywords.EmEpPurkinjeNetwork2(
                        purkid=2,
                        buildnet=0,
                        ssid=network.nsid,
                        mid=network.pid,
                        pointstx=origin_coordinates[0],
                        pointsty=origin_coordinates[1],
                        pointstz=origin_coordinates[2],
                        edgelen=2,
                        ngen=50,
                        nbrinit=8,
                        nsplit=2,
                        # inodeid=node_id_start_right,
                        # iedgeid=edge_id_start_right,
                    )
                )
                part_df = pd.DataFrame(
                    {
                        "heading": [network.name],
                        "pid": [network.pid],
                        "secid": [3],
                        "mid": [network.pid],
                    }
                )
                part_kw = keywords.Part()
                part_kw.parts = part_df
                self.kw_database.parts.append(part_kw)
                self.kw_database.material.append(keywords.MatNull(mid=network.pid, ro=1e-11))
                self.kw_database.material.append(
                    keywords.EmMat001(mid=network.pid, mtype=2, sigma=10)
                )

                beams_kw = add_beams_to_kw(
                    beams=network.edges + 1,
                    beam_kw=beams_kw,
                    pid=network.pid,
                    offset=len(self.model.mesh.tetrahedrons) + len(beams_kw.elements),
                )
                cell_kw = keywords.EmEpCellmodelTentusscher(
                    mid=network.pid,
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
                self.kw_database.cell_models.extend([cell_kw])
            self.kw_database.beam_networks.append(beams_kw)

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


class ElectroMechanicsDynaWriter(MechanicsDynaWriter, ElectrophysiologyDynaWriter):
    """Class for preparing the input for LS-DYNA electromechanical simulation."""

    def __init__(
        self,
        model: HeartModel,
        settings: SimulationSettings = None,
    ) -> None:
        super().__init__(model=model, settings=settings)

        raise NotImplementedError("This writer has not been implemented yet.")
        exit()
        self.kw_database = ElectroMechanicsDecks()
        """Collection of keyword decks relevant for mechanics."""

        self.system_model_name = system_model_name
        """Name of system model to use."""

        # Depending on the system model specified give list of parameters

        return

    def update(self):
        """Update the keyword database."""
        self._update_node_db()
        self._update_parts_db()
        self._update_main_db()
        self._update_solid_elements_db(add_fibers=True)
        self._update_segmentsets_db()
        self._update_nodesets_db()
        self._update_use_Purkinje()
        self._update_material_db(add_active=True)
        self._update_ep_material_db()
        self._update_cellmodels()
        self._update_ep_settings()
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

    def _update_material_db(self, add_active: bool = True):
        """Update the database of material keywords."""
        material_settings = copy.deepcopy(self.settings.mechanics.material)
        # removes all units from settings, hence <attribute>.m not required anymore to access value.
        material_settings._remove_units()

        for part in self.model.parts:
            if "ventricle" in part.name.lower() or "septum" in part.name.lower():
                if not add_active:
                    active_dct = None
                else:
                    active_dct = {
                        "actype": material_settings.myocardium["active"]["actype"],
                        "taumax": material_settings.myocardium["active"]["tmax"],
                        "ca2ionm": material_settings.myocardium["active"]["ca2ionm"],
                    }

                myocardium_kw = MaterialHGOMyocardium(
                    mid=part.mid,
                    iso_user=material_settings.myocardium["isotropic"],
                    anisotropy_user=material_settings.myocardium["anisotropic"],
                    active_user=active_dct,
                )

                self.kw_database.material.append(myocardium_kw)

            elif "atrium" in part.name:
                # add atrium material
                # atrium_kw = MaterialAtrium(mid=part.mid)
                atrium_kw = MaterialHGOMyocardium(
                    mid=part.mid, iso_user=dict(material_settings.atrium)
                )

                self.kw_database.material.append(atrium_kw)

            else:
                LOGGER.warning("Assuming same material as atrium for: {0}".format(part.name))

                # general_tissue_kw = MaterialAtrium(mid=part.mid)
                general_tissue_kw = MaterialHGOMyocardium(
                    mid=part.mid, iso_user=dict(material_settings.atrium)
                )
                self.kw_database.material.append(general_tissue_kw)

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
        elif self.type == "ra_fiber":
            parts_to_keep = ["Right atrium"]
            #  A manual point for RA fiber
            for key, value in kwargs.items():
                if key == "raa":
                    self.right_appendage_apex = value

        # remove unnecessary mesh
        if self.type == "uvc":
            elems_to_keep = []
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
            LOGGER.error("Cannot find top node set for right atrium.")
            exit()
        for i in range(3):
            share_nodes = np.any(
                np.logical_and(x.point_data["tricuspid-valve"] == 1, x.point_data["RegionId"] == i)
            )
            # This region has no shared node with tricuspid valve
            if not share_nodes:
                mask = x.point_data["RegionId"] == i
                break

        top_ids = x["point_ids_tmp"][mask]

        atrium.cell_data.remove("cell_ids_tmp")
        atrium.point_data.remove("point_ids_tmp")

        # assign
        kw = create_node_set_keyword(top_ids + 1, node_set_id=10, title="top")
        self.kw_database.node_sets.append(kw)
        atrium["top"] = np.zeros(atrium.n_points)
        atrium["top"][top_ids] = 1

        # Find tricuspid_wall and tricuspid_septum
        id_sorter = np.argsort(atrium["point_ids"])
        # need a copied object to do clip, atrium will be corrupted otherwise
        septum, free_wall = copy.deepcopy(atrium).clip(
            origin=cut_center, normal=cut_normal, crinkle=True, return_clipped=True
        )
        # ids in full mesh
        tv_s_ids = septum["point_ids"][np.where(septum["tricuspid-valve"] == 1)]

        tv_s_ids_sub = id_sorter[np.searchsorted(atrium["point_ids"], tv_s_ids, sorter=id_sorter)]
        atrium["tv_s"] = np.zeros(atrium.n_points)
        atrium["tv_s"][tv_s_ids_sub] = 1

        kw = create_node_set_keyword(tv_s_ids_sub + 1, node_set_id=12, title="tv_s")
        self.kw_database.node_sets.append(kw)

        tv_w_ids = free_wall["point_ids"][np.where(free_wall["tricuspid-valve"] == 1)]
        tv_w_ids_sub = id_sorter[np.searchsorted(atrium["point_ids"], tv_w_ids, sorter=id_sorter)]
        # remove re constraint nodes
        tv_w_ids_sub = np.setdiff1d(tv_w_ids_sub, tv_s_ids_sub)

        atrium["tv_w"] = np.zeros(atrium.n_points)
        atrium["tv_w"][tv_w_ids_sub] = 1

        kw = create_node_set_keyword(tv_w_ids_sub + 1, node_set_id=13, title="tv_w")
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

            return set_id

        id_sorter = np.argsort(atrium["point_ids"])
        ids_edges = []
        for i, cap in enumerate(self.model.parts[0].caps):
            # node IDs in LA volume mesh
            ids_sub = id_sorter[
                np.searchsorted(atrium["point_ids"], cap.node_ids, sorter=id_sorter)
            ]
            set_id = get_nodeset_id_by_cap_name(cap)

            kw = create_node_set_keyword(ids_sub + 1, node_set_id=set_id, title=cap.name)
            self.kw_database.node_sets.append(kw)

            ids_edges.extend(ids_sub)
            atrium[cap.name] = np.zeros(atrium.n_points, dtype=int)
            atrium[cap.name][ids_sub] = i + 1

        # endo nodes ID
        ids_endo = id_sorter[
            np.searchsorted(
                atrium["point_ids"], self.model.parts[0].surfaces[0].node_ids, sorter=id_sorter
            )
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

            self.kw_database.main.append("*CASE_BEGIN_5")
            self._define_Laplace_Dirichlet_bc(set_ids=[12, 13], bc_values=[1.0, -1.0])
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
        id_sorter = np.argsort(self.target["point_ids"])

        endo_set = []
        epi_set = []
        for part in self.model.parts:
            for surf in part.surfaces:
                if "endocardium" in surf.name:
                    endo_set.extend(surf.node_ids)
                # elif "epicardium" in surf.name:
                #     epi_set.extend(surf.node_ids)

        # map IDs to sub mesh
        endo_set_new = id_sorter[
            np.searchsorted(self.target["point_ids"], endo_set, sorter=id_sorter)
        ]
        endo_set_new = np.unique(endo_set_new)

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
        apex_set = self.model._compute_uvc_apex_set()
        id_sorter = np.argsort(self.target["point_ids"])
        ids_submesh = id_sorter[
            np.searchsorted(self.target["point_ids"], apex_set, sorter=id_sorter)
        ]
        sid = self.get_unique_nodeset_id()
        kw = create_node_set_keyword(ids_submesh + 1, node_set_id=sid, title="apex")
        self.kw_database.node_sets.append(kw)
        return sid

    def _create_base_nodeset(self):
        # base
        base_set = np.array([])
        for part in self.model.parts:
            for cap in part.caps:
                if ("mitral" in cap.name) or ("tricuspid" in cap.name):
                    base_set = np.append(base_set, cap.node_ids)

        id_sorter = np.argsort(self.target["point_ids"])
        ids_submesh = id_sorter[
            np.searchsorted(self.target["point_ids"], base_set, sorter=id_sorter)
        ]
        sid = self.get_unique_nodeset_id()
        kw = create_node_set_keyword(ids_submesh + 1, node_set_id=sid, title="base")
        self.kw_database.node_sets.append(kw)
        return sid

    def _create_surface_nodeset(self, surftype: str, cavity_type: str):
        id_sorter = np.argsort(self.target["point_ids"])

        nodeset = np.array([])
        for part in self.model.parts:
            if cavity_type in part.name:
                for surf in part.surfaces:
                    if surftype in surf.name:
                        nodeset = np.append(nodeset, surf.node_ids)
        nodeset = np.unique(nodeset.astype(int))

        # map IDs to sub mesh
        ids_submesh = id_sorter[
            np.searchsorted(self.target["point_ids"], nodeset, sorter=id_sorter)
        ]

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
            part.pid = self.get_unique_part_id()
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
