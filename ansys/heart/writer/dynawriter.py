"""Contains class for writing dyna keywords based on the HeartModel
"""
import numpy as np
import pandas as pd
import os
import time
import json
from pathlib import Path
from tqdm import tqdm  # for progress bar

from ansys.heart.preprocessor._deprecated_heart_model import HeartModel
from ansys.heart.preprocessor._deprecated_cavity_module import ClosingCap

from ansys.heart.custom_logging import LOGGER
from ansys.heart.preprocessor.mesh.vtkmethods import (
    get_tetra_info_from_unstructgrid,
    vtk_surface_filter,
    compute_surface_nodal_area,
)

from ansys.heart.writer.keyword_module import (
    add_nodes_to_kw,
    create_discrete_elements_kw,
    create_element_solid_ortho_keyword,
    create_element_shell_keyword,
    create_segment_set_keyword,
    create_node_set_keyword,
    create_discrete_elements_kw,
    create_define_curve_kw,
    create_define_sd_orientation_kw,
    fast_element_writer,
    get_list_of_used_ids,
)

# import commonly used material models
from ansys.heart.writer.material_keywords import (
    MaterialCap,
    MaterialHGOMyocardium,
    MaterialAtrium,
    active_curve,
)

from ansys.heart.writer.heart_decks import (
    BaseDecks,
    MechanicsDecks,
    FiberGenerationDecks,
    PurkinjeGenerationDecks,
)

from vtk.numpy_interface import dataset_adapter as dsa  # noqa

from ansys.dyna.keywords import keywords

# import missing keywords
from ansys.heart.writer import custom_dynalib_keywords as custom_keywords


class BaseDynaWriter:
    """BaseDynaWriter class which contains features essential
    for all relevant LS-DYNA heart models
    """

    def __init__(self, model: HeartModel) -> None:
        """Initializes writer by loading a HeartModel

        Parameters
        ----------
        model : HeartModel
            HeartModel object which contains the necessary
            information for the writer, such as nodes and elements.
        """

        self.model = model
        """Contains model information necessary for creating the LS-DYNA .k files"""

        self.kw_database = BaseDecks()

        # These are general attributes useful for keeping track of ids:
        self.max_node_id: int = 0
        """Max node id"""
        self.part_ids = []
        """List of used part ids"""
        self.section_ids = []
        """List of used section ids"""
        self.mat_ids = []
        """List of used mat ids"""

        self.volume_mesh = {
            "nodes": np.empty(0),
            "tetra": np.empty(0),
            "cell_data": {},
            "point_data": {},
        }
        """Volume mesh information"""

        # keeps track of some element id offsets
        self.id_offset = {
            "part": 0,
            "section": 0,
            "material": 0,
            "vector": 0,
            "element": {"solid": 0, "discrete": 0, "shell": 0},
        }
        """Id offset for several relevant keywords"""

        self.include_files = []
        """List of .k files to include in main. This is derived from the Decks classes"""

        # read mesh information into dictionary
        self._get_mesh_info()

        if "Improved" in self.model.info.model_type:
            LOGGER.warning(
                "Changing model type from : {0} to {1}".format(
                    self.model.info.model_type, self.model.info.model_type.replace("Improved", "")
                )
            )
            self.model.info.model_type = self.model.info.model_type.replace("Improved", "")

        return

    def _get_mesh_info(self):
        """Gets nodes, element definition, cell data and point data from the
        volume mesh
        """
        nodes, tetra, cell_data, point_data = get_tetra_info_from_unstructgrid(
            self.model._mesh._vtk_volume
        )
        self.volume_mesh["nodes"] = nodes
        self.volume_mesh["tetra"] = tetra
        self.volume_mesh["cell_data"] = cell_data
        self.volume_mesh["point_data"] = point_data

        return

    def _get_list_of_includes(self):
        """Gets a list of files to include in main.k. Ommit any empty decks"""
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
        """Adds *INCLUDE keywords"""
        for include_file in self.include_files:
            filename_to_include = include_file + ".k"
            self.kw_database.main.append(keywords.Include(filename=filename_to_include))

        return

    def export_databases(self, export_directory: str):
        """Exports each of non-empty databases to a specified directory"""

        if not export_directory:
            export_directory = self.model.info.working_directory

        for deckname, deck in vars(self.kw_database).items():
            # skip empty databases:
            if deck.keywords == []:
                continue

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

            else:
                deck.export_file(filepath)
        return


class MechanicsDynaWriter(BaseDynaWriter):
    """Derived from BaseDynaWriter and derives all keywords relevant
    for simulations involving mechanics"""

    def __init__(self, model: HeartModel, system_model_name: str = "ClosedLoop") -> None:
        super().__init__(model)

        self.kw_database = MechanicsDecks()
        """Collection of keyword decks relevant for mechanics"""

        self.system_model_name = system_model_name
        """Name of system model to use"""

        # Depending on the system model specified give list of parameters

        return

    @property
    def system_model_name(self):
        """System model name. Valid options include:
        ["ConstantPreloadWindkesselAfterload",
        "ClosedLoop]"""
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
        """Formats the keywords and stores these in the
        respective keyword databases
        """

        self._update_main_db()
        self._update_node_db()
        self._update_parts_db()
        self._update_solid_elements_db(add_fibers=True)
        self._update_segmentsets_db()
        self._update_nodesets_db()
        self._update_material_db(add_active=True)

        # for boundary conditions
        # self._update_boundary_conditions_db()
        self._add_cap_bc(bc_type="springs_caps")
        self._add_pericardium_bc()
        # self._add_pericardium_bc_usr()

        # for control volume
        self._update_cap_elements_db()
        self._update_controlvolume_db()
        self._update_system_model()

        # Approximate end-diastolic pressures
        pressure_lv = 2  # kPa
        pressure_rv = 0.5333  # kPa

        self._add_enddiastolic_pressure_bc(pressure_lv=pressure_lv, pressure_rv=pressure_rv)

        self._get_list_of_includes()
        self._add_includes()

        return

    def export(self, export_directory: str):
        """Writes the model to files"""
        tstart = time.time()
        LOGGER.debug("Writing all LS-DYNA .k files...")

        if not export_directory:
            export_directory = self.model.info.working_directory

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
        LOGGER.debug("Time spend writing files: {:.2f} s".format(tend - tstart))

        return

    def _update_main_db(self, add_damping: bool = True):
        """Updates the main .k file
        Note
        -----
        Consider using a settings (json?) file as input
        """
        LOGGER.debug("Updating main keywords...")

        self.kw_database.main.title = self.model.info.model_type

        self._add_solution_controls()
        self._add_export_controls()

        if add_damping:
            self._add_damping()

        return

    def _update_node_db(self):
        """Adds nodes to the Node database"""
        LOGGER.debug("Updating node keywords...")
        node_kw = keywords.Node()
        node_kw = add_nodes_to_kw(self.volume_mesh["nodes"], node_kw)

        self.kw_database.nodes.append(node_kw)

        return

    def _update_parts_db(self):
        """Creates database of PART keywords.
        Each element set within cavity defined as separate part. Each
        cavity associated with one material
        """

        LOGGER.debug("Updating part keywords...")
        # add parts with a dataframe
        part_ids = []
        part_id = 0
        mat_id = 0
        for cavity in self.model._mesh._cavities:
            mat_id = mat_id + 1
            for element_set in cavity.element_sets:
                part_id = part_id + 1
                part_name = " ".join([cavity.name, element_set["name"]])
                part = pd.DataFrame(
                    {"heading": [part_name], "pid": [part_id], "secid": [1], "mid": [mat_id]}
                )
                part_kw = keywords.Part()
                part_kw.parts = part

                self.kw_database.parts.append(part_kw)

                # store part id for future use
                element_set["id"] = part_id
                element_set["mid"] = mat_id
                part_ids.append(part_id)

        # set up section solid for cavity myocardium
        section_kw = keywords.SectionSolid(secid=1, elform=13)

        self.kw_database.parts.append(section_kw)

        return

    def _update_solid_elements_db(self, add_fibers: bool = True):
        """Creates Solid ortho elements for all cavities.
        Each cavity contains one myocardium which consists corresponds
        to one part.
        """
        LOGGER.debug("Updating solid element keywords...")
        cell_data_fields = self.volume_mesh["cell_data"].keys()
        if "fiber" not in cell_data_fields or "sheet" not in cell_data_fields:
            raise KeyError("Mechanics writer requires fiber and sheet fields")
            # logger.warning("Not writing fiber and sheet directions")
            # add_fibers = False

        # create elements for each separate cavity
        solid_element_count = 0  # keeps track of number of solid elements already defined

        for cavity in self.model._mesh._cavities:
            LOGGER.debug("Writing elements for myocardium" " attached to cavity: " + cavity.name)

            # get list of elements to write to the database. Create new keyword
            # for each part. Parts defined by element sets
            for element_set in cavity.element_sets:
                tetra_idx = element_set["set"]
                part_id = element_set["id"]

                tetra_to_write = self.volume_mesh["tetra"][tetra_idx, :] + 1

                num_elements = len(tetra_idx)

                element_ids = np.arange(1, num_elements + 1, 1) + solid_element_count
                part_ids = np.ones(num_elements, dtype=int) * part_id

                # format the element keywords
                if not add_fibers:
                    kw_elements = keywords.ElementSolid()
                    elements = pd.DataFrame(
                        {
                            "eid": element_ids,
                            "pid": part_ids,
                            "n1": tetra_to_write[:, 0],
                            "n2": tetra_to_write[:, 1],
                            "n3": tetra_to_write[:, 2],
                            "n4": tetra_to_write[:, 3],
                            "n5": tetra_to_write[:, 3],
                            "n6": tetra_to_write[:, 3],
                            "n7": tetra_to_write[:, 3],
                            "n8": tetra_to_write[:, 3],
                        }
                    )
                    kw_elements.elements = elements

                elif add_fibers:
                    fiber = self.volume_mesh["cell_data"]["fiber"][tetra_idx]
                    sheet = self.volume_mesh["cell_data"]["sheet"][tetra_idx]

                    # normalize fiber and sheet directions:
                    norm = np.linalg.norm(fiber, axis=1)
                    fiber = fiber / norm[:, None]
                    norm = np.linalg.norm(sheet, axis=1)
                    sheet = sheet / norm[:, None]
                    kw_elements = create_element_solid_ortho_keyword(
                        elements=tetra_to_write,
                        a_vec=fiber,
                        d_vec=sheet,
                        partid=part_id,
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
        """Adds solution controls, output controls and other solver settings
        as keywords
        """
        # add termination keywords
        self.kw_database.main.append(keywords.ControlTermination(endtim=end_time, dtmin=dtmin))

        # add implict controls
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
                "Simulation type not recoqnized: Please choose " "either quasi-static or static"
            )

        self.kw_database.main.append(
            keywords.ControlImplicitDynamics(imass=imass, gamma=gamma, beta=beta)
        )

        # add auto controls
        self.kw_database.main.append(
            keywords.ControlImplicitAuto(iauto=1, dtmin=dtmin, dtmax=dtmax)
        )

        # add general implicit controls
        self.kw_database.main.append(
            keywords.ControlImplicitGeneral(imflag=1, dt0=dtmin)
        )  # imflag=1 means implicit

        # add implicit solution controls: Defaults are OK?
        self.kw_database.main.append(keywords.ControlImplicitSolution())

        # add implicit solver controls
        self.kw_database.main.append(keywords.ControlImplicitSolver())
        return

    def _add_export_controls(self, dt_output_d3plot: float = 0.05, dt_output_icvout: float = 0.001):
        """Adds solution controls to the main simulation

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
        self.kw_database.main.append(keywords.DatabaseBinaryD3Plot(dt=dt_output_d3plot))

        self.kw_database.main.append(keywords.DatabaseExtentBinary(neiph=27, strflg=1, maxint=0))

        return

    def _add_damping(self):
        """Adds damping to the main file"""
        lcid_damp = 3

        kw_damp = keywords.DampingGlobal(lcid=lcid_damp)

        kw_damp_curve = create_define_curve_kw(
            x=[0, 100],
            y=[500, 500],
            curve_name="damping",
            curve_id=lcid_damp,
            lcint=0,
        )
        self.kw_database.main.append(kw_damp)
        self.kw_database.main.append(kw_damp_curve)

        return

    def _update_segmentsets_db(self):
        """Updates the segment set database"""
        # write each segment set to the segment set database

        # NOTE 1: need to more robustly check segids that are already used?
        segment_set_id = 100
        for cavity in self.model._mesh._cavities:
            for segset in cavity.segment_sets:
                segset_name = " ".join([cavity.name, segset["name"]])
                kw = create_segment_set_keyword(
                    segments=segset["set"] + 1,
                    segid=segment_set_id,
                    title=segset_name,
                )
                # append this kw to the segment set database
                self.kw_database.segment_sets.append(kw)

                # add assigned segment set id to the model
                segset["id"] = segment_set_id
                segment_set_id = segment_set_id + 1

        return

    def _update_nodesets_db(self):
        """Updates the node set database"""
        # formats endo, epi- and septum nodeset keywords
        # do for all cavities and for all caps that are defined
        # append each of the keywords to the nodesets database
        nodeset_id = 1
        node_sets_visited = []
        for cavity in self.model._mesh._cavities:
            for nodeset in cavity.node_sets:
                set_name = " ".join([cavity.name, nodeset["name"]])
                kw = create_node_set_keyword(
                    nodeset["set"] + 1, node_set_id=nodeset_id, title=set_name
                )
                self.kw_database.node_sets.append(kw)
                # add assigned node set id to the model
                nodeset["id"] = nodeset_id
                nodeset_id = nodeset_id + 1

        nodeset_id_offset = 100
        # NOTE: could use same structure for nodeset as <class Cavity>
        set_id_visited = []
        for cavity in self.model._mesh._cavities:
            for cap in cavity.closing_caps:
                nodeset_id = cap.id + nodeset_id_offset
                if nodeset_id in set_id_visited:
                    # skip this node set
                    continue
                cap.nodeset_id = nodeset_id

                kw = create_node_set_keyword(
                    cap.node_ids_cap_edge + 1, cap.nodeset_id, title=cap.name
                )
                self.kw_database.node_sets.append(kw)
                set_id_visited.append(cap.nodeset_id)

    def _update_material_db(self, add_active: bool = True):
        """Updates the database of material keywords"""
        for cavity in self.model._mesh._cavities:
            # element set id of first element set
            # assumed correspond to material id (one material per cavity)
            # NOTE: in case of writing the zero-pressure input files the
            # active module should be off
            mat_id = cavity.element_sets[0]["mid"]

            # curve id for active module
            act_curve_id = 15

            if "ventricle" in cavity.name:
                # add ventricular materials
                myocardium_material_kw = MaterialHGOMyocardium(
                    mid=mat_id, add_anisotropy=True, add_active=add_active
                )
                myocardium_material_kw.acid = act_curve_id
                self.kw_database.material.append(myocardium_material_kw)

            if "atrium" in cavity.name:
                # add arterial material
                atrium_material = MaterialAtrium(mid=mat_id)
                self.kw_database.material.append(atrium_material)
                print("")

        if add_active:
            # write and add active curve to material database
            time_array, active_stress_array = active_curve("Strocchi2020")
            active_curve_kw = create_define_curve_kw(
                x=time_array,
                y=active_stress_array,
                curve_name="calcium_concentration",
                curve_id=act_curve_id,
                lcint=15000,
            )

            active_curve_kw.sfo = 4.35  # y scaling
            active_curve_kw.offa = 1.00  # x offset

            self.kw_database.material.append(active_curve_kw)

        return

    def _update_boundary_conditions_db(self):
        """Updates the boundary conditions keyword database"""

        # self._add_cap_bc(bc_type="fix_all_caps")
        pass
        return

    def _add_cap_bc(self, bc_type: str):
        """Adds boundary condition to the cap.

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
        if self.model.info.model_type in ["LeftVentricle", "BiVentricle"]:
            # use all caps:
            # for cavity in self.model._mesh._cavities:
            #     for cap in cavity.closing_caps:
            #         caps_to_use.append(cap.name)
            caps_to_use = [
                "Mitral valve plane",
                "Tricuspid valve plane",
            ]

        elif self.model.info.model_type in ["FourChamber"]:
            caps_to_use = [
                "Superior vena cava inlet",
                "Right inferior pulmonary vein inlet",
                "Right superior pulmonary vein inlet",
            ]

        if bc_type == "fix_caps":
            for cavity in self.model._mesh._cavities:
                for cap in cavity.closing_caps:
                    if cap.name in caps_to_use:
                        kw_fix = keywords.BoundarySpcSet()
                        kw_fix.nsid = cap.nodeset_id
                        kw_fix.dofx = 1
                        kw_fix.dofy = 1
                        kw_fix.dofz = 1

                        self.kw_database.boundary_conditions.append(kw_fix)

        # if bc type is springs -> add springs
        # NOTE add to boundary condition db or seperate spring db?
        elif bc_type == "springs_caps":

            # NOTE: Make dynamic and expose to user?

            # NOTE: Need to be made dynamic
            part_id = 200
            section_id = 200
            mat_id = 200

            # TODO: exposed to user/parameters?
            if self.model.info.model_type == "BiVentricle":
                spring_stiffness = 5  # kPa/mm
            elif self.model.info.model_type == "FourChamber":
                spring_stiffness = 20  # kPa/mm

            scale_factor_normal = 0.5
            scale_factor_radial = 1.0

            part_kw = keywords.Part()
            part = pd.DataFrame(
                {
                    "pid": [part_id],
                    "secid": [section_id],
                    "mid": [mat_id],
                    "heading": ["SupportSpring"],
                }
            )
            part_kw.parts = part

            section_kw = keywords.SectionDiscrete(secid=section_id, cdl=0, tdl=0)

            mat_kw = keywords.MatSpringElastic(mid=mat_id, k=spring_stiffness)

            self.kw_database.boundary_conditions.append(part_kw)
            self.kw_database.boundary_conditions.append(section_kw)
            self.kw_database.boundary_conditions.append(mat_kw)

            # add springs for each cap
            for cavity in self.model._mesh._cavities:
                for cap in cavity.closing_caps:
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
        cap: ClosingCap,
        part_id: int,
        scale_factor_normal: float,
        scale_factor_radial: float,
    ):
        """Adds springs to the cap nodes and appends these
        to the boundary condition database
        """
        # NOTE: may want to extent the node ids to include adjacent nodes
        num_nodes_edge = len(cap.node_ids_cap_edge)
        # -------------------------------------------------------------------

        # compute nodal areas:
        vtk_surface = vtk_surface_filter(self.model._mesh._vtk_volume, True)
        nodal_areas = compute_surface_nodal_area(vtk_surface)
        surface_obj = dsa.WrapDataObject(vtk_surface)
        surface_global_node_ids = surface_obj.PointData["GlobalPointIds"]

        # select only those nodal areas which match the cap node ids
        idx_select = np.isin(surface_global_node_ids, cap.node_ids_cap_edge)
        nodal_areas = nodal_areas[idx_select]

        # scaled spring stiffness by nodal area
        scale_factor_normal *= nodal_areas
        scale_factor_radial *= nodal_areas

        # -------------------------------------------------------------------
        LOGGER.debug("Adding spring b.c. for cap: %s" % cap.name)

        # add part, section discrete, mat spring, sd_orientiation
        # element discrete

        # compute the radial components
        sd_orientations_radial = cap.nodes_cap_edge - cap.centroid

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
            [cap.node_ids_cap_edge + 1, np.zeros(num_nodes_edge)], dtype=int
        ).T
        vector_ids_normal = np.ones(num_nodes_edge, dtype=int) * vector_id_normal

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
        """Adds the pericardium. Uses the universal ventricular longitudinal coordinate
        and a sigmoid penalty function. Strocchi et al 2020 doi: 10.1016/j.jbiomech.2020.109645.
        """

        def _sigmoid(z):
            """sigmoid function to scale spring coefficient"""
            return 1 / (1 + np.exp(-z))

        # compute penalty function
        uvc_l = self.volume_mesh["point_data"]["uvc_longitudinal"]

        if np.any(uvc_l < 0):
            LOGGER.warning(
                "Negative normalized longitudinal coordinate detected."
                "Changing {0} negative uvc_l values to 1".format(np.sum((uvc_l < 0))),
            )

        uvc_l[uvc_l < 0] = 1
        penalty = -_sigmoid((abs(uvc_l) - 0.1) * 25) + 1

        # collect all pericardium nodes:
        epicardium_nodes = np.empty(0, dtype=int)
        LOGGER.debug("Collecting epicardium nodesets:")
        for cavity in self.model._mesh._cavities:
            if cavity.name == "Right ventricle" or cavity.name == "Left ventricle":
                for nodeset in cavity.node_sets:
                    if nodeset["name"] == "epicardium":
                        LOGGER.debug("\t{0} {1}".format(cavity.name, nodeset["name"]))
                        epicardium_nodes = np.append(epicardium_nodes, nodeset["set"])

        # select only nodes that are on the epicardium and penalty factor > 0.1
        pericardium_nodes = epicardium_nodes[penalty[epicardium_nodes] > 0.001]
        # coord = self.volume_mesh["nodes"][pericardium_nodes]
        # np.savetxt(r"pericardium.txt",
        #            np.concatenate((coord,penalty[pericardium_nodes].reshape(-1,1)),axis=1))

        # TODO: exposed to user/parameters?
        spring_stiffness = 50  # kPA/mm
        # compute nodal areas:
        vtk_surface = vtk_surface_filter(self.model._mesh._vtk_volume, True)
        nodal_areas = compute_surface_nodal_area(vtk_surface)

        surface_obj = dsa.WrapDataObject(vtk_surface)
        surface_global_node_ids = surface_obj.PointData["GlobalPointIds"]

        # select only those nodal areas which match the pericardium node ids
        idx_select = np.isin(surface_global_node_ids, pericardium_nodes)
        nodal_areas = nodal_areas[idx_select]

        # compute scale factor
        scale_factors = nodal_areas * penalty[pericardium_nodes]

        # keywords
        # NOTE: Need to be made dynamic
        part_id = 201
        section_id = 201
        mat_id = 201

        part_kw = keywords.Part()
        part_kw.parts = pd.DataFrame(
            {"heading": ["Pericardium"], "pid": [part_id], "secid": [section_id], "mid": [mat_id]}
        )
        section_kw = keywords.SectionDiscrete(secid=section_id, cdl=0, tdl=0)
        mat_kw = keywords.MatSpringElastic(mid=mat_id, k=spring_stiffness)

        # 1: "omni-directional": equal springs in x,y, and z
        # 2: "apex-mitral-drection": one spring in apex-mitral valve direction
        spring_types = ["omni-directional", "apex-mitral-direction"]
        spring_type = "apex-mitral-direction"
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
            for cavity in self.model._mesh._cavities:
                if cavity.name == "Left ventricle":
                    apex_node_id = cavity.apex_id["epicardium"]
                    apex1 = self.volume_mesh["nodes"][apex_node_id, :]
                    for cap in cavity.closing_caps:
                        if cap.name == "Mitral valve plane":
                            center1 = cap.centroid
                        elif cap.name == "Aortic valve plane":
                            center2 = cap.centroid
            # Change orientation as apex- center of 2 valves plane
            center = (center2 + center1) / 2
            # define spring orientation from apex to mitral valve
            orientation = center - apex1
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

    def _add_pericardium_bc_usr(self):
        """Adds the pericardium.
        Same with _add_pericardium_bc but using *user_load, need customized LSDYNA exe!!
        """

        def _sigmoid(z):
            """sigmoid function to scale spring coefficient"""
            return 1 / (1 + np.exp(-z))

        # compute penalty function
        uvc_l = self.volume_mesh["point_data"]["uvc_longitudinal"]

        if np.any(uvc_l < 0):
            LOGGER.warning(
                "Negative normalized longitudinal coordinate detected. Changing {0}"
                "negative uvc_l values to 1".format(np.sum((uvc_l < 0))),
            )

        uvc_l[uvc_l < 0] = 1
        penalty = -_sigmoid((abs(uvc_l) - 0.15) * 25) + 1  # for all volume nodes

        # collect all pericardium nodes:
        epicardium_segment = np.empty((0, 3), dtype=int)

        LOGGER.debug("Collecting epicardium nodesets:")
        for cavity in self.model._mesh._cavities:
            if cavity.name == "Right ventricle" or cavity.name == "Left ventricle":
                for sgm_set in cavity.segment_sets:
                    if sgm_set["name"] == "epicardium":
                        LOGGER.debug("\t{0} {1}".format(cavity.name, sgm_set["name"]))
                        epicardium_segment = np.vstack((epicardium_segment, sgm_set["set"]))

        penalty = np.mean(
            penalty[epicardium_segment], axis=1
        )  # averaged for all pericardium segments

        # # debug code
        # # export pericardium segment center and also the penalty factor, can be opened in Paraview
        # coord = self.volume_mesh["nodes"][epicardium_segment]
        # center = np.mean(coord, axis=1)
        # result = np.concatenate((center,penalty.reshape(-1,1)),axis=1)
        # np.savetxt('pericardium.txt',result[result[:,3]>0.01])
        # # exit()
        # # end debug code

        # create load curve to control when pericardium is active
        load_curve_kw = keywords.DefineCurve(lcid=4)
        load_curve_kw.options["TITLE"].active = True
        load_curve_kw.title = "pericardium activation curve"
        load_curve_kw.curves = pd.DataFrame(
            {"a1": np.array([0, 1, 100]), "o1": np.array([1, 1, 1])}
        )
        self.kw_database.pericardium.append(load_curve_kw)

        cnt = 0
        load_sgm_kws = []
        segment_ids = []
        LOGGER.debug("Creating segment sets for epicardium b.c.:")

        penalty_threshold = 0.01
        for isg, sgmt in enumerate(tqdm(epicardium_segment, ascii=True)):
            if penalty[isg] > penalty_threshold:
                cnt += 1

                # coord = self.volume_mesh["nodes"][sgmt]
                # center = np.mean(coord, axis=0)
                # normal = np.cross(coord[1] - coord[0], coord[2] - coord[0])
                # normal /= np.linalg.norm(normal)
                # cs_kw = keywords.DefineCoordinateSystem(
                #     cid=cnt,
                #     xo=center[0],
                #     yo=center[1],
                #     zo=center[2],
                #     xl=center[0] + normal[0],
                #     yl=center[1] + normal[1],
                #     zl=center[2] + normal[2],
                #     xp=coord[0, 0],
                #     yp=coord[0, 1],
                #     zp=coord[0, 2],
                # )

                segment_id = 1000 + cnt
                segment_ids.append(segment_id)

                load_sgm_kw = create_segment_set_keyword(
                    segments=sgmt.reshape(1, -1) + 1, segid=segment_id
                )  # todo: auto counter

                load_sgm_kws.append(load_sgm_kw)

        self.kw_database.pericardium.extend(load_sgm_kws)

        # create user loadset keyword
        # segment_ids = 1000 + np.arange(0, np.sum( penalty > 0.01 ), 1)
        user_loadset_kw = custom_keywords.UserLoadingSet()

        # NOTE: can assign mixed scalar/array values to dataframe - scalars are assigned to each row
        user_loadset_kw.load_sets = pd.DataFrame(
            {
                "sid": segment_ids,
                "ltype": "PRESSS",
                "lcid": 4,
                "sf1": penalty[penalty > penalty_threshold],
                "iduls": 100,
            }
        )
        user_load_kw = custom_keywords.UserLoading(parm1=10.0)

        # add to pericardium deck
        self.kw_database.pericardium.extend([user_loadset_kw, user_load_kw])

        return

    def _update_cap_elements_db(self):
        """Updates the database of shell elements. Loops over all
        the defined caps/valves
        """
        # create part for each closing cap
        used_partids = get_list_of_used_ids(self.kw_database.parts, "PART")
        used_secids = get_list_of_used_ids(self.kw_database.parts, "SECTION")
        used_segids = get_list_of_used_ids(self.kw_database.segment_sets, "SET_SEGMENT")

        part_id = np.max(used_partids) + 1
        section_id = np.max(used_secids) + 1

        # NOTE should be dynamic
        mat_null_id = 100

        # material_kw = MaterialCap(mid=mat_null_id)

        material_kw = MaterialAtrium(mid=mat_null_id, rho=1e-6, poisson_ratio=0.499, c10=1000)

        section_kw = keywords.SectionShell(
            secid=section_id,
            elform=4,
            shrf=0.8333,
            nip=3,
            t1=5,
            t2=5,
            t3=5,
            t4=5,
        )

        self.kw_database.cap_elements.append(material_kw)
        self.kw_database.cap_elements.append(section_kw)

        # create new part for each cap
        for cavity in self.model._mesh._cavities:
            for cap in cavity.closing_caps:
                cap.part_id = part_id
                part_kw = keywords.Part()
                part_info = pd.DataFrame(
                    {
                        "heading": [cap.name],
                        "pid": [cap.part_id],
                        "secid": [section_id],
                        "mid": [mat_null_id],
                    }
                )
                part_kw.parts = part_info

                self.kw_database.cap_elements.append(part_kw)

                part_id = part_id + 1

        # create closing triangles for each cap
        # assumes there are no shells written yet since offset = 0
        shell_id_offset = 0
        for cavity in self.model._mesh._cavities:
            for cap in cavity.closing_caps:
                shell_kw = create_element_shell_keyword(
                    shells=cap.closing_triangles + 1,
                    part_id=cap.part_id,
                    id_offset=shell_id_offset,
                )

                self.kw_database.cap_elements.append(shell_kw)

                shell_id_offset = shell_id_offset + cap.closing_triangles.shape[0]

        # create corresponding segment sets. Store in new file?
        segset_id = np.max(used_segids) + 1

        for cavity in self.model._mesh._cavities:
            for cap in cavity.closing_caps:
                segset_kw = create_segment_set_keyword(
                    segments=cap.closing_triangles + 1,
                    segid=segset_id,
                    title=cap.name,
                )

                self.kw_database.cap_elements.append(segset_kw)

                cap.segset_id = segset_id
                segset_id = segset_id + 1

        return

    def _update_controlvolume_db(self):
        """Prepares the keywords for the control volume feature"""
        # NOTE: Assumes cavity id is reseverd for combined
        # segment set

        # set up control volume keywords and interaction of
        # cavity with ambient. Only do for ventricles
        for cavity in self.model._mesh._cavities:
            if "atrium" in cavity.name:
                continue

            cv_kw = keywords.DefineControlVolume()
            cv_kw.id = cavity.id
            cv_kw.sid = cavity.id

            self.kw_database.control_volume.append(cv_kw)

        for cavity in self.model._mesh._cavities:
            if "atrium" in cavity.name:
                continue

            cvi_kw = keywords.DefineControlVolumeInteraction()
            cvi_kw.id = cavity.id
            cvi_kw.cvid1 = cavity.id
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
        """Updates json system model settings"""
        model_type = self.model.info.model_type

        # closed loop uses a custom executable
        if self.system_model_name == "ClosedLoop":
            LOGGER.warning(
                "Note that this model type requires a custom executable that supports the Closed Loop circulation model!"
            )
            if model_type in ["FourChamber", "BiVentricle"]:
                file_path = os.path.join(
                    Path(__file__).parent.absolute(),
                    "templates",
                    "system_model_settings_bv.json",
                )

            elif model_type in ["LeftVentricle"]:
                file_path = os.path.join(
                    Path(__file__).parent.absolute(),
                    "templates",
                    "system_model_settings_lv.json",
                )

            fid = open(file_path)
            sys_settings = json.load(fid)

            # update the volumes
            for cavity in self.model._mesh._cavities:
                if "Left ventricle" in cavity.name:
                    sys_settings["SystemModelInitialValues"]["UnstressedVolumes"][
                        "lv"
                    ] = cavity.volume
                elif "Right ventricle" in cavity.name:
                    sys_settings["SystemModelInitialValues"]["UnstressedVolumes"][
                        "rv"
                    ] = cavity.volume

            self.system_model_json = sys_settings

        # otherwise add the define function
        elif self.system_model_name == "ConstantPreloadWindkesselAfterload":
            from ansys.heart.writer.system_models import define_function_windkessel

            for cavity in self.model._mesh._cavities:
                if "Left ventricle" in cavity.name:
                    constants: dict = {
                        "Rv": 5.0e-6,
                        "Ra": 1.0e-5,
                        "Rp": 1.2e-4,
                        "Ca": 2.5e4,
                        "Pven": 2,
                    }
                    initial = {"part_init": 8}
                    define_function_wk = define_function_windkessel(
                        function_id=10,
                        function_name="constant_preload_windkessel_afterload_left",
                        implicit=True,
                        constants=constants,
                        initialvalues=initial,
                    )
                    self.kw_database.control_volume.append(define_function_wk)

                elif "Right ventricle" in cavity.name:
                    constants: dict = {
                        "Rv": 5.0e-6 * 0.5,
                        "Ra": 1.0e-5 * 0.35,
                        "Rp": 1.2e-4 * 0.125,
                        "Ca": 2.5e4 * 4.5,
                        "Pven": 0.53333,
                    }
                    initial = {"part_init": 2}
                    define_function_wk = define_function_windkessel(
                        function_id=11,
                        function_name="constant_preload_windkessel_afterload_right",
                        implicit=True,
                        constants=constants,
                        initialvalues=initial,
                    )
                    self.kw_database.control_volume.append(define_function_wk)

        return

    def _export_cavity_segmentsets(self, export_directory: str):
        """Exports the actual cavity segment sets to separate files"""

        for cavity in self.model._mesh._cavities:
            filename = "cavity_" + "_".join(cavity.name.lower().split()) + ".segment"
            filepath = os.path.join(export_directory, filename)

            segments = np.empty((0, 3), dtype=int)

            # collect segment sets
            for segset in cavity.segment_sets:
                if segset["name"] in ["endocardium", "endocardium-septum"]:
                    segments = np.vstack([segments, segset["set"]])

            # append cap segments:
            for cap in cavity.closing_caps:
                segments = np.vstack([segments, cap.closing_triangles])

            # combine segment sets
            np.savetxt(filepath, segments, delimiter=",", fmt="%d")

        return

    def _add_enddiastolic_pressure_bc(self, pressure_lv: float = 1, pressure_rv: float = 1):
        """Adds end diastolic pressure boundary condition on the left and right endocardium"""

        # collect segment set ids to combine
        for cavity in self.model._mesh._cavities:
            seg_ids_to_combine = []
            # find id of endocardium
            for segset in cavity.segment_sets:
                if "endocardium" in segset["name"]:
                    seg_ids_to_combine.append(segset["id"])

            for cap in cavity.closing_caps:
                seg_ids_to_combine.append(cap.segset_id)

            # add segment set add keyword
            segadd_kw = keywords.SetSegmentAdd(sid=cavity.id)
            segadd_kw.sets._data = seg_ids_to_combine
            segadd_kw.options["TITLE"].active = True
            segadd_kw.title = cavity.name

            self.kw_database.main.append(segadd_kw)

        # create unit load curve
        load_curve_id = 2
        load_curve_kw = create_define_curve_kw(
            [0, 1, 1.001], [0, 1, 0], "unit load curve", load_curve_id, 100
        )

        # append unit curve to main.k
        self.kw_database.main.append(load_curve_kw)

        # create *LOAD_SEGMENT_SETS for each ventricular cavity
        for cavity in self.model._mesh._cavities:

            if "atrium" in cavity.name:
                continue

            if cavity.name == "Left ventricle":
                scale_factor = pressure_lv
                seg_id = 1
            elif cavity.name == "Right ventricle":
                scale_factor = pressure_rv
                seg_id = 2
            load_segset_kw = keywords.LoadSegmentSet(
                ssid=seg_id, lcid=load_curve_id, sf=scale_factor
            )
            self.kw_database.main.append(load_segset_kw)

            # logger.debug(
            #     "Adding end-diastolic pressure of {0} to {1}".format(scale_factor, cavity.name)
            # )
            #
            # seg_ids_to_use = []
            # # find id of endocardium
            # for segset in cavity.segment_sets:
            #     if "endocardium" in segset["name"]:
            #         seg_ids_to_use.append(segset["id"])
            #
            # # create load segment set for each endocardium segment
            # for seg_id in [1,2]:
            #     load_segset_kw = keywords.LoadSegmentSet(
            #         ssid=seg_id, lcid=load_curve_id, sf=scale_factor
            #     )
            #
            #     # append to main.k
            #     self.kw_database.main.append(load_segset_kw)


class ZeroPressureMechanicsDynaWriter(MechanicsDynaWriter):
    """Derived from MechanicsDynaWriter and consequently derives all keywords relevant
    for simulations involving mechanics. This class does not use write the
    control volume keywords but adds the keyword for computing the stress
    free configuration based on left/right cavity pressures instead"""

    def __init__(self, model: HeartModel) -> None:
        super().__init__(model)

        self.kw_database = MechanicsDecks()
        """Collection of keyword decks relevant for mechanics"""

        return

    def update(self):
        """Formats the keywords and stores these in the
        respective keyword databases. Overwrites the update method
        of MechanicsDynaWriter such that it yields a valid input deck
        for a zero-pressure simulation
        """

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
        # pressure_lv = 2  # kPa
        # pressure_rv = 0.5333  # kPa
        # self._add_enddiastolic_pressure_bc(pressure_lv=pressure_lv, pressure_rv=pressure_rv)

        # zerop key words
        self._add_control_reference_configuration()

        self._get_list_of_includes()
        self._add_includes()

        return

    def export(self, export_directory: str):
        """Writes the model to files"""
        tstart = time.time()
        LOGGER.debug("Writing all LS-DYNA .k files...")

        if not export_directory:
            export_directory = self.model.info.working_directory

        # export .k files
        self.export_databases(export_directory)

        tend = time.time()
        LOGGER.debug("Time spend writing files: {:.2f} s".format(tend - tstart))

        return

    def _add_export_controls(self, dt_output_d3plot: float = 0.5):
        """Rewrite method for zerop export

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
        """
        Rewrite method for the zerop simulation
        Returns
        -------

        """
        self.kw_database.main.append(keywords.ControlTermination(endtim=1.0))
        # self.kw_database.main.append(keywords.ControlImplicitDynamics(imass=0))

        self.kw_database.main.append(
            keywords.ControlImplicitDynamics(imass=1, gamma=0.6, beta=0.38)
        )

        # add auto controls
        self.kw_database.main.append(keywords.ControlImplicitAuto(iauto=1, dtmin=0.01, dtmax=0.1))

        # add general implicit controls
        self.kw_database.main.append(keywords.ControlImplicitGeneral(imflag=1, dt0=0.1))

        # add implicit solution controls: Defaults are OK?
        self.kw_database.main.append(keywords.ControlImplicitSolution())

        # add implicit solver controls
        self.kw_database.main.append(keywords.ControlImplicitSolver())

        # add binout for post-process
        self.kw_database.main.append(keywords.DatabaseNodout(dt=0.5, binary=1))
        x = keywords.SetNodeGeneral(option="part", sid=999, e1=1, e2=2, e3=3, e4=4)
        self.kw_database.main.append(x)
        self.kw_database.main.append("*DATABASE_HISTORY_NODE_SET\n999")
        return

    def _add_control_reference_configuration(self):
        """Adds control reference configuration keyword to main"""
        LOGGER.debug("Adding *CONTROL_REFERENCE_CONFIGURATION to main.k")
        kw = keywords.ControlReferenceConfiguration(maxiter=3, target="nodes.k", method=2, tol=5)

        self.kw_database.main.append(kw)

        return

    def _add_enddiastolic_pressure_bc(self, pressure_lv: float = 1, pressure_rv: float = 1):
        """Adds end diastolic pressure boundary condition on the left and right endocardium"""

        # create unit load curve
        load_curve_id = 2
        load_curve_kw = create_define_curve_kw(
            [0, 1], [0, 1], "unit load curve", load_curve_id, 100
        )

        # append unit curve to main.k
        self.kw_database.main.append(load_curve_kw)

        # create *LOAD_SEGMENT_SETS for each ventricular cavity
        for cavity in self.model._mesh._cavities:

            if "atrium" in cavity.name:
                continue

            if cavity.name == "Left ventricle":
                scale_factor = pressure_lv
            elif cavity.name == "Right ventricle":
                scale_factor = pressure_rv

            LOGGER.debug(
                "Adding end-diastolic pressure of {0} to {1}".format(scale_factor, cavity.name)
            )

            seg_ids_to_use = []
            # find id of endocardium
            for segset in cavity.segment_sets:
                if "endocardium" in segset["name"]:
                    seg_ids_to_use.append(segset["id"])

            # create load segment set for each endocardium segment
            for seg_id in seg_ids_to_use:
                load_segset_kw = keywords.LoadSegmentSet(
                    ssid=seg_id, lcid=load_curve_id, sf=scale_factor
                )

                # append to main.k
                self.kw_database.main.append(load_segset_kw)


class FiberGenerationDynaWriter(MechanicsDynaWriter):
    def __init__(self, model: HeartModel) -> None:
        super().__init__(model)

        self.kw_database = FiberGenerationDecks()
        """Collection of keywords relevant for fiber generation
        """

    def update(self):
        """Updates keyword database for Fiber generation: overwrites the inherited function"""

        ##
        self._update_main_db()  # needs updating

        self._update_node_db()  # can stay the same (could move to base class)
        if self.model.info.model_type == "FourChamber":
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
        self._update_create_fibers()

        self._get_list_of_includes()
        self._add_includes()

        return

    def export(self, export_directory: str):
        """Writes the model to files"""
        tstart = time.time()
        LOGGER.debug("Writing all LS-DYNA .k files...")

        if not export_directory:
            export_directory = self.model.info.working_directory

        # export .k files
        self.export_databases(export_directory)

        tend = time.time()
        LOGGER.debug("Time spend writing files: {:.2f} s".format(tend - tstart))

        return

    def _keep_ventricles(self):
        """Removes any cavity except the ventricular cavities"""
        # just keep ventricles in case of four chamber model
        LOGGER.warning("Just keeping ventricular-parts for fiber generation")
        cavities_to_keep = []
        for ii, cavity in enumerate(self.model._mesh._cavities):
            if "ventricle" in cavity.name:
                cavities_to_keep.append(cavity)

        self.model._mesh._cavities = cavities_to_keep
        return

    def _update_material_db(self):
        """Adds simple linear elastic and orthotropic EM material for each defined part"""
        for cavity in self.model._mesh._cavities:
            for element_set in cavity.element_sets:
                em_mat_id = element_set["id"]
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
                            lambda_=0.5,
                            a1=0,
                            a2=0,
                            a3=0,
                            d1=0,
                            d2=-1,
                            d3=0,
                        ),
                        custom_keywords.EmEpCellmodelTomek(mid=em_mat_id),
                    ]
                )

    def _update_ep_settings(self):
        """Adds the settings for the electrophysiology solver"""

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
        """Updates the keywords for fiber generation"""

        # collect relevant node and segment sets.
        # node set: apex, base
        # node set: endocardium, epicardium
        # NOTE: could be better if basal nodes are extracted in the preprocessor
        # since that would allow you to robustly extract these nodessets using the
        # input data
        # The below is relevant for all models.
        nodes_base = np.empty(0, dtype=int)
        node_set_ids_endo = []
        node_sets_ids_epi = []
        node_set_ids_epi_and_rseptum = []

        node_set_id_lv_endo = 0

        for cavity in self.model._mesh._cavities:
            if "atrium" in cavity.name:
                continue

            for cap in cavity.closing_caps:
                nodes_base = np.append(nodes_base, cap.node_ids_cap_edge)

            for node_set in cavity.node_sets:
                if "endocardium" in node_set["name"]:
                    node_set_ids_endo.append(node_set["id"])
                    if "septum" in node_set["name"]:
                        node_set_ids_epi_and_rseptum.append(node_set["id"])
                    if cavity.name == "Left ventricle":
                        node_set_id_lv_endo = node_set["id"]
                elif "epicardium" in node_set["name"]:
                    node_sets_ids_epi.append(node_set["id"])
                    node_set_ids_epi_and_rseptum.append(node_set["id"])

            if cavity.name == "Left ventricle":
                node_apex = np.array([cavity.apex_id["epicardium"]])

        # validate node set by removing any nodes that do not occur in either ventricle
        # NOTE: can be much more consice
        tet_ids_ventricles = np.empty((0), dtype=int)
        for cavity in self.model._mesh._cavities:
            for element_set in cavity.element_sets:
                if "ventricle" in cavity.name:
                    tet_ids_ventricles = np.append(tet_ids_ventricles, element_set["set"])
        tetra_ventricles = self.volume_mesh["tetra"][tet_ids_ventricles, :]

        # remove nodes that occur just in atrial part
        mask = np.isin(nodes_base, tetra_ventricles, invert=True)
        LOGGER.debug("Removing {0} nodes from base nodes".format(np.sum(mask)))
        nodes_base = nodes_base[np.invert(mask)]

        # create set parts for lv and rv myocardium
        myocardium_part_ids = []
        septum_part_ids = []
        for cavity in self.model._mesh._cavities:
            if "atrium" in cavity.name:
                continue
            for segset in cavity.element_sets:
                if segset["name"] == "myocardium":
                    myocardium_part_ids.append(segset["id"])

                if segset["name"] == "septum":
                    septum_part_ids.append(segset["id"])

        # switch between the various models to generate valid input decks
        if self.model.info.model_type in ["LeftVentricle"]:
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
            node_set_id_all_endocardium = 1000
            set_add_kw = keywords.SetNodeAdd(sid=node_set_id_all_endocardium)
            set_add_kw.options["TITLE"].active = True
            set_add_kw.title = "all_endocardium_segments"
            set_add_kw.nodes._data = node_set_ids_endo

            self.kw_database.create_fiber.append(set_add_kw)

            # combine node sets epicardium:
            node_set_id_all_epicardium = 1001
            set_add_kw = keywords.SetNodeAdd(sid=node_set_id_all_epicardium)
            set_add_kw.options["TITLE"].active = True
            set_add_kw.title = "all_epicardium_segments"
            set_add_kw.nodes._data = node_sets_ids_epi

            self.kw_database.create_fiber.append(set_add_kw)

            node_set_id_base = 200
            node_set_id_apex = 201
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

        elif self.model.info.model_type in ["BiVentricle", "FourChamber"]:
            LOGGER.warning("Model type %s under development " % self.model.info.model_type)

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
            node_set_id_all_endocardium = 1000
            set_add_kw = keywords.SetNodeAdd(sid=node_set_id_all_endocardium)

            set_add_kw.options["TITLE"].active = True
            set_add_kw.title = "all_endocardium_segments"
            set_add_kw.nodes._data = node_set_ids_endo

            self.kw_database.create_fiber.append(set_add_kw)

            # combine node sets epicardium:
            node_set_id_all_epicardium = 1001
            set_add_kw = keywords.SetNodeAdd(sid=node_set_id_all_epicardium)

            set_add_kw.options["TITLE"].active = True
            set_add_kw.title = "all_epicardium_segments"
            set_add_kw.nodes._data = node_sets_ids_epi

            self.kw_database.create_fiber.append(set_add_kw)

            # combine node sets epicardium and septum:
            node_set_all_but_left_endocardium = 1002
            set_add_kw = keywords.SetNodeAdd(sid=node_set_all_but_left_endocardium)

            set_add_kw.options["TITLE"].active = True
            set_add_kw.title = "all_but_left_endocardium"
            set_add_kw.nodes._data = node_set_ids_epi_and_rseptum

            self.kw_database.create_fiber.append(set_add_kw)

            node_set_id_base = 200
            node_set_id_apex = 201
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
    def __init__(self, model: HeartModel) -> None:
        super().__init__(model)

        self.kw_database = PurkinjeGenerationDecks()
        """Collection of keywords relevant for Purkinje generation
        """

    def update(self):
        """Updates keyword database for Purkinje generation: overwrites the inherited function"""

        ##
        self._update_main_db()  # needs updating

        self._update_node_db()  # can stay the same (could move to base class)
        if self.model.info.model_type == "FourChamber":
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
        """Writes the model to files"""
        tstart = time.time()
        LOGGER.debug("Writing all LS-DYNA .k files...")

        if not export_directory:
            export_directory = self.model.info.working_directory

        # export .k files
        self.export_databases(export_directory)

        tend = time.time()
        LOGGER.debug("Time spend writing files: {:.2f} s".format(tend - tstart))

        return

    def _keep_ventricles(self):
        """Removes any cavity except the ventricular cavities"""
        # just keep ventricles in case of four chamber model
        LOGGER.warning("Just keeping ventricular-parts for Purkinje generation")
        cavities_to_keep = []
        for ii, cavity in enumerate(self.model._mesh._cavities):
            if "ventricle" in cavity.name:
                cavities_to_keep.append(cavity)

        self.model._mesh._cavities = cavities_to_keep
        return

    def _update_material_db(self):
        """Adds simple linear elastic material for each defined part"""
        for cavity in self.model._mesh._cavities:
            for element_set in cavity.element_sets:
                em_mat_id = element_set["id"]
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
                            lambda_=0.5,
                            a1=0,
                            a2=0,
                            a3=0,
                            d1=0,
                            d2=-1,
                            d3=0,
                        ),
                    ]
                )

    def _update_ep_settings(self):
        """Adds the settings for the electrophysiology solver"""

        self.kw_database.ep_settings.append(
            keywords.EmControl(
                emsol=11, numls=4, macrodt=1, dimtype=None, nperio=None, ncylbem=None
            )
        )

        self.kw_database.ep_settings.append(keywords.EmOutput(mats=1, matf=1, sols=1, solf=1))

        return

    def _update_create_Purkinje(self):
        """Updates the keywords for Purkinje generation"""

        # collect relevant node and segment sets.
        # node set: apex, base
        # node set: endocardium, epicardium
        # NOTE: could be better if basal nodes are extracted in the preprocessor
        # since that would allow you to robustly extract these nodessets using the
        # input data
        # The below is relevant for all models.
        nodes_base = np.empty(0, dtype=int)
        node_apex_left = np.empty(0, dtype=int)
        node_apex_right = np.empty(0, dtype=int)
        edge_id_start_left = np.empty(0, dtype=int)
        edge_id_start_right = np.empty(0, dtype=int)
        for cavity in self.model._mesh._cavities:
            if cavity.name == "Left ventricle":
                node_apex_left = cavity.apex_id["endocardium"]
                for segment_set in cavity.segment_sets:
                    if "endocardium" in segment_set["name"]:
                        segment_set_ids_endo_left = segment_set["id"]
            elif cavity.name == "Right ventricle":
                node_apex_right = cavity.apex_id["endocardium"]
                for segment_set in cavity.segment_sets:
                    if "endocardium" in segment_set["name"] and "septum" not in segment_set["name"]:
                        segment_set_ids_endo_right = segment_set["id"]

        # validate node set by removing any nodes that do not occur in either ventricle
        # NOTE: can be much more consice
        tet_ids_ventricles = np.empty((0), dtype=int)
        for cavity in self.model._mesh._cavities:
            for element_set in cavity.element_sets:
                if "ventricle" in cavity.name:
                    tet_ids_ventricles = np.append(tet_ids_ventricles, element_set["set"])
        tetra_ventricles = self.volume_mesh["tetra"][tet_ids_ventricles, :]

        # remove nodes that occur just in atrial part
        mask = np.isin(nodes_base, tetra_ventricles, invert=True)
        LOGGER.debug("Removing {0} nodes from base nodes".format(np.sum(mask)))
        nodes_base = nodes_base[np.invert(mask)]

        node_set_id_apex_left = 201
        # create node-sets for apex
        node_set_apex_kw = create_node_set_keyword(
            node_ids=[node_apex_left + 1],
            node_set_id=node_set_id_apex_left,
            title="apex node left",
        )

        self.kw_database.node_sets.extend([node_set_apex_kw])

        apex_left_X = self.volume_mesh["nodes"][node_apex_left, 0]
        apex_left_Y = self.volume_mesh["nodes"][node_apex_left, 1]
        apex_left_Z = self.volume_mesh["nodes"][node_apex_left, 2]
        node_id_start_left = (
            self.volume_mesh["nodes"].shape[0] + 1
        )  # TODO seek for max id rather than number of rows

        edge_id_start_left = self.volume_mesh["tetra"].shape[0] + 1

        # Purkinje generation parameters
        self.kw_database.mainLEFT.append(
            custom_keywords.EmEpPurkinjeNetwork(
                purkid=1,
                buildnet=1,
                ssid=segment_set_ids_endo_left,
                mid=25,
                pointstx=apex_left_X,
                pointsty=apex_left_Y,
                pointstz=apex_left_Z,
                edgelen=2,
                ngen=50,
                nbrinit=8,
                nsplit=2,
                inodeid=node_id_start_left,
                iedgeid=edge_id_start_left,  # TODO check if beam elements exist in mesh
            )
        )

        # Add right purkinje only in biventricular or 4chamber models
        if self.model.info.model_type in ["BiVentricle", "FourChamber"]:
            LOGGER.warning("Model type %s in development " % self.model.info.model_type)

            node_set_id_apex_right = 202
            # create node-sets for apex
            node_set_apex_kw = create_node_set_keyword(
                node_ids=[node_apex_right + 1],
                node_set_id=node_set_id_apex_right,
                title="apex node right",
            )

            self.kw_database.node_sets.extend([node_set_apex_kw])

            apex_right_X = self.volume_mesh["nodes"][node_apex_right, 0]
            apex_right_Y = self.volume_mesh["nodes"][node_apex_right, 1]
            apex_right_Z = self.volume_mesh["nodes"][node_apex_right, 2]
            node_id_start_right = (
                2 * self.volume_mesh["nodes"].shape[0]
            )  # TODO find a solution in dyna to better handle id definition

            edge_id_start_right = 2 * self.volume_mesh["tetra"].shape[0]

            # Purkinje generation parameters
            self.kw_database.mainRIGHT.append(
                custom_keywords.EmEpPurkinjeNetwork(
                    purkid=2,
                    buildnet=1,
                    ssid=segment_set_ids_endo_right,
                    mid=26,
                    pointstx=apex_right_X,
                    pointsty=apex_right_Y,
                    pointstz=apex_right_Z,
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
        """Gets a list of files to include in main.k. Ommit any empty decks"""
        for deckname, deck in vars(self.kw_database).items():
            if deckname == "mainLEFT" or deckname == "mainRIGHT":
                continue
            # skip if no keywords are present in the deck
            if len(deck.keywords) == 0:
                LOGGER.debug("No keywords in deck: {0}".format(deckname))
                continue
            self.include_files.append(deckname)
        return

    def _add_includes(self):
        """Adds *INCLUDE keywords"""
        for include_file in self.include_files:
            filename_to_include = include_file + ".k"
            self.kw_database.mainLEFT.append(keywords.Include(filename=filename_to_include))
            self.kw_database.mainRIGHT.append(keywords.Include(filename=filename_to_include))


if __name__ == "__main__":
    print("protected")
