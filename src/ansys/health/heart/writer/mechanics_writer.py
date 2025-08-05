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
"""Module containing classes for writing LS-DYNA keyword files mechanics simulations."""

import copy
from enum import Enum
from typing import Callable, Literal, Optional

import numpy as np
import pandas as pd
import pyvista as pv

from ansys.dyna.core.keywords import keywords
from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.models import BiVentricle, FourChamber, FullHeart, HeartModel, LeftVentricle
from ansys.health.heart.objects import Cap, CapType, SurfaceMesh
import ansys.health.heart.parts as anatomy
from ansys.health.heart.settings.material.material import (
    Mat295,
    MechanicalMaterialModel,
    NeoHookean,
)
import ansys.health.heart.settings.settings as sett
from ansys.health.heart.settings.settings import SimulationSettings
from ansys.health.heart.utils.vtk_utils import compute_surface_nodal_area_pyvista
from ansys.health.heart.writer import custom_keywords as custom_keywords
from ansys.health.heart.writer._control_volume import (
    ControlVolume,
    _create_closed_loop,
    _create_open_loop,
)
from ansys.health.heart.writer.base_writer import BaseDynaWriter
from ansys.health.heart.writer.heart_decks import MechanicsDecks
from ansys.health.heart.writer.material_keywords import MaterialHGOMyocardium, MaterialNeoHook
from ansys.health.heart.writer.writer_utils import (
    create_define_curve_kw,
    create_define_sd_orientation_kw,
    create_discrete_elements_kw,
    create_element_shell_keyword,
)


class _BoundaryConditionType(Enum):
    """Boundary condition type."""

    FIX = "fix"
    ROBIN = "Robin"


class MechanicsDynaWriter(BaseDynaWriter):
    """Class for preparing the input for a mechanics LS-DYNA simulation."""

    def __init__(
        self,
        model: HeartModel,
        settings: Optional[SimulationSettings] = None,
    ) -> None:
        super().__init__(model=model, settings=settings)

        self.kw_database = MechanicsDecks()
        """Collection of keyword decks relevant for mechanics."""

        self.set_flow_area: bool = True
        """Flag indicating if the flow area is set for control volume."""

        if sett.Mechanics not in self._get_subsettings():
            raise ValueError("Expecting mechanics settings.")
        return

    def update(self, dynain_name: Optional[str] = None, robin_bcs: list[Callable] = None) -> None:
        """Update the keyword database.

        Parameters
        ----------
        dynain_name : str, default: None
            Dynain file from stress-free configuration computation.
        robin_bcs : list[Callable], default: None
            List of lambda functions to apply Robin-type coundary conditions.

        Notes
        -----
        You do not need to write mesh files if a Dynain file is given.
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
            self.include_to_main(dynain_name)
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

        if system_settings.name == "open-loop":
            lcid = self.get_unique_curve_id()
            system_map = _create_open_loop(lcid, self.model, system_settings)
        elif system_settings.name == "closed-loop":
            LOGGER.warning("Closed loop uses a recompiled version of LS-DYNA!")
            system_map = _create_closed_loop(self.model)
        else:
            msg = r"System name must be `open-loop` or `closed-loop`"
            LOGGER.error(msg)
            raise TypeError(msg)

        self._update_controlvolume_db(system_map)

        include_files = self._get_decknames_of_include()
        self.include_to_main(include_files)

        return

    def _update_main_db(self) -> None:
        """Update the main K file."""
        LOGGER.debug("Updating main keywords...")

        self.kw_database.main.append("$$- Unit system: g-mm-ms-N-MPa-mJ -$$")
        self.kw_database.main.title = self.model.__class__.__name__

        if isinstance(self, ZeroPressureMechanicsDynaWriter):
            settings = self.settings.stress_free
            self._add_solution_controls()
            self._add_export_controls(settings.analysis.dt_d3plot.m)

        elif isinstance(self, MechanicsDynaWriter):
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
    ) -> None:
        """Add solution controls, output controls, and solver settings."""
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
                "Simulation type is not recognized: Choose either 'quasi-static' or 'static'."
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

    def _add_export_controls(
        self, dt_output_d3plot: float = 0.05, dt_output_icvout: float = 0.001
    ) -> None:
        """Add solution controls to the main simulation.

        Parameters
        ----------
        dt_output_d3plot : float, default: 0.5
            Time-step spacing to write full D3PLOT results at.
        dt_output_icvout : float, default: 0.001
            Time-step spacing to write control volume results at.
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

    def _add_damping(self) -> None:
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

    def _update_material_db(self, add_active: bool = True, em_couple: bool = False) -> None:
        #
        for part in self.model.parts:
            if isinstance(part.meca_material, MechanicalMaterialModel.DummyMaterial):
                # assign material for part if it's empty
                LOGGER.info(f"Material of {part.name} is assigned automatically.")
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

    def _add_cap_bc(self, bc_type: _BoundaryConditionType) -> None:
        """Add boundary condition to the cap.

        Parameters
        ----------
        bc_type : BoundaryType
           Boundary condition type.

        """
        # create list of cap names where to add the spring b.c
        constraint_caps = self._get_contraint_caps()

        self.model.all_caps

        if bc_type == _BoundaryConditionType.FIX:
            for cap in self.model.all_caps:
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
            for cap in self.model.all_caps:
                if cap.type in constraint_caps:
                    self.kw_database.boundary_conditions.append(f"$$ spring at {cap.name}$$")
                    self._add_springs_cap_edge(
                        cap,
                        part_id,
                        scale_factor_normal,
                        scale_factor_radial,
                    )

        return

    def _get_contraint_caps(self) -> list[CapType]:
        """Get a list of constraint caps, depending on models."""
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
    ) -> None:
        """Add springs to the cap nodes.

        Notes
        -----
        This method appends these springs to the boundary condition database.
        """
        LOGGER.debug(f"Adding spring boundary condition for cap: {cap.name} of type {cap.type}")

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

        # create discrete elements
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

    def _add_pericardium_bc(self, scale: float = 1.0) -> None:
        """Add the pericardium."""
        boundary_conditions = copy.deepcopy(self.settings.mechanics.boundary_conditions)
        robin_settings = boundary_conditions.robin

        # collect all pericardium nodes:
        ventricles_epi = self._get_epi_surface(apply=anatomy._PartType.VENTRICLE)

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
            atrial_epi = self._get_epi_surface(anatomy._PartType.ATRIUM)

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
        self,
        apply: Literal[
            anatomy._PartType.VENTRICLE, anatomy._PartType.ATRIUM
        ] = anatomy._PartType.VENTRICLE,
    ) -> SurfaceMesh:
        """Get the epicardial surfaces of either the ventricle or atria."""
        LOGGER.debug(f"Collecting epicardium nodesets of {apply}:")

        targets = [part for part in self.model.parts if apply == part._part_type]

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

    def _get_longitudinal_penalty(self, pericardium_settings: dict) -> np.ndarray:
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
                "No apico-basal is found in point data. Pericardium spring won't be created."
            )
            uvc_l = np.ones(self.model.mesh.GetNumberOfPoints())
        if np.any(uvc_l < 0):
            LOGGER.warning(
                "Negative normalized longitudinal coordinate is detected."
                "Changing {0} negative uvc_l values to 1.".format(np.sum((uvc_l < 0))),
            )
        uvc_l[uvc_l < 0] = 1

        penalty_function = -_sigmoid((abs(uvc_l) - penalty_c0) * penalty_c1) + 1
        return penalty_function

    def write_robin_bc(
        self,
        robin_type: Literal["spring", "damper"],
        constant: float,
        surface: pv.PolyData,
        normal: Optional[np.ndarray] = None,
    ) -> list:
        """Create Robin boundary condition on a given surface.

        Parameters
        ----------
        robin_type : Literal["spring", "damper"]
            Create spring or damper.
        constant : float
            Stiffness (MPa/mm) or viscosity (MPa/mm*ms).
        surface : pv.PolyData
            Surface to apply boundary condition to. It must contain point data
            ``_global-point-ids``. It is scaled by the nodal area and point data
            scale factor if it exists.
        normal : np.ndarray, default: None
            Normal values. If no normal values are given, nodal normals are used.

        Returns
        -------
        list
            List of the DYNA input deck.
        """
        if surface.n_points == 0:
            LOGGER.error("Surface is empty. No Robin boundary condition is added.")
            return []

        if "_global-point-ids" not in surface.point_data:
            raise ValueError("Surface must contain point data '_global-point-ids'.")

        # global node ids where to apply the BC
        # NOTE: if we pass in a SurfaceMesh object we could use the
        # .global_node_ids attribute instead.
        nodes = surface["_global-point-ids"]

        # scale factor is nodal area
        # Add area flag in case PyVista defaults change.
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

        # create unique IDs for keywords
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

    def _update_cap_elements_db(self, add_mesh: bool = True) -> None:
        """Update the database of shell elements.

        Notes
        -----
        This method loops over all the defined caps and valves.
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

        # create new part for each cap
        cap_names_used = []
        for cap in self.model.all_caps:
            if cap.name in cap_names_used:
                # avoid to write mitral valve and triscupid valve twice
                LOGGER.debug("Already created material for {}. Skipping.".format(cap.name))
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
                    LOGGER.error("Cap nodeset ID is not yet assigned.")
                    raise ValueError("Cap nodeset ID is not yet assigned.")

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
            for cap in self.model.all_caps:
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

    def _update_controlvolume_db(self, system_map: list[ControlVolume]) -> None:
        """Prepare the keywords for the control volume feature.

        Parameters
        ----------
        system_map : list[ControlVolume]
            List of control volumes.
        """

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
                if len(sets) % 8 == 0:  # PyDYNA keywords bug when length is 8,16,...
                    sets.append(0)
                self.kw_database.control_volume.append(keywords.SetSegmentAdd(sid=sid, sets=sets))

                # TODO: use PyDYNA keywords: keywords.DefineControlVolumeFlowArea()
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
                self.kw_database.control_volume.append(interaction._define_function_keyword())

        return


class ZeroPressureMechanicsDynaWriter(MechanicsDynaWriter):
    """
    Class for preparing the input for a stress-free LS-DYNA simulation.

    Notes
    -----
    This class is derived from the ``MechanicsDynaWriter`` class and consequently
    derives all keywords relevant for simulations involving mechanics. This class
    does not write the control volume keywords but rather adds the keyword for computing
    the stress-free configuration based on left/right cavity pressures instead.

    """

    def __init__(
        self,
        model: HeartModel,
        settings: Optional[SimulationSettings] = None,
    ) -> None:
        super().__init__(model=model, settings=settings)

        self.kw_database = MechanicsDecks()
        """Collection of keyword decks relevant for mechanics."""

        return

    def update(self, robin_bcs: list[Callable] = None) -> None:
        """Update the keyword database.

        Parameters
        ----------
        robin_bcs : list[Callable], default: None
            List of lambda functions to apply Robin-type boundary conditions.
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

        for cap in self.model.all_caps:
            if cap.pid is not None:  # MV,TV for atrial parts get None
                save_part_ids.append(cap.pid)

        partset_id = self.get_unique_partset_id()
        kw = keywords.SetPartList(sid=partset_id)
        # kw.parts._data = save_part_ids
        # NOTE: when len(save_part_ids) = 8/16, PyDYNA keywords bugs
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

        include_files = self._get_decknames_of_include()
        self.include_to_main(include_files)

        return

    def _add_export_controls(self, dt_output_d3plot: float = 0.5) -> None:
        """Rewrite the method for zerop export.

        Parameters
        ----------
        dt_output_d3plot : float, default: 0.5
            Time-space spacing to write full D3PLOT results at.
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

    def _add_solution_controls(self) -> None:
        """Rewrite the method for the zerop simulation."""
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

    def _add_control_reference_configuration(self) -> None:
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
    #     LS-DYNA stress reference configuration leads to a bug with this load.
    #     It seems due to define function and must be investigated.
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

    def _add_enddiastolic_pressure_bc(self) -> None:
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
        cavities = self.model.cavities
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
