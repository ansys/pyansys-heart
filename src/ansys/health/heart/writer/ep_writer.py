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
"""Module containing classes for writing LS-DYNA keyword files for electrophysiology simulations."""

from typing import Union

import numpy as np
import pandas as pd
import pyvista as pv
import scipy.spatial as spatial

from ansys.dyna.core.keywords import keywords
from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.models import BiVentricle, FourChamber, FullHeart, HeartModel, LeftVentricle
from ansys.health.heart.models_utils import LandMarks
from ansys.health.heart.pre.conduction_path import ConductionPathType
from ansys.health.heart.settings.material.ep_material import CellModel, EPMaterial
import ansys.health.heart.settings.settings as sett
from ansys.health.heart.settings.settings import SimulationSettings, Stimulation
from ansys.health.heart.writer import custom_keywords as custom_keywords
from ansys.health.heart.writer.base_writer import BaseDynaWriter
from ansys.health.heart.writer.heart_decks import ElectrophysiologyDecks, PurkinjeGenerationDecks
from ansys.health.heart.writer.writer_utils import (
    add_beams_to_kw,
    add_nodes_to_kw,
    create_node_set_keyword,
    create_segment_set_keyword,
)


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

        if sett.Purkinje not in self._get_subsettings():
            raise ValueError("Expecting Purkinje settings.")

    def update(self) -> None:
        """Update keyword database.

        This method overwrites the inherited function.
        """
        ##
        self._update_main_db()  # needs updating

        self._update_node_db()  # can stay the same (could move to base class)
        if isinstance(self.model, (FourChamber, FullHeart)):
            LOGGER.warning(
                "Atrium are present in the model. "
                "These are removed for ventricle Purkinje generation."
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

        include_files = self._get_decknames_of_include()
        self.include_to_main(include_files)

        return

    def _update_material_db(self) -> None:
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

    def _update_ep_settings(self) -> None:
        """Add the settings for the electrophysiology solver."""
        self.kw_database.ep_settings.append(
            keywords.EmControl(
                emsol=11, numls=4, macrodt=1, dimtype=None, nperio=None, ncylbem=None
            )
        )

        self.kw_database.ep_settings.append(keywords.EmOutput(mats=1, matf=1, sols=1, solf=1))

        return

    def _update_create_Purkinje(self) -> None:  # noqa N802
        """Update the keywords for Purkinje generation."""
        # collect relevant node and segment sets.
        # nodeset: apex, base
        # nodeset: endocardium, epicardium
        # NOTE: could be better if basal nodes are extracted in the preprocessor
        # since that would allow you to robustly extract these nodessets using the
        # input data
        # What follows is relevant for all models.

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

    def _update_main_db(self) -> None:
        return


class ElectrophysiologyDynaWriter(BaseDynaWriter):
    """Class for preparing the input for an electrophysiology LS-DYNA simulation."""

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
        """Collection of keywords relevant for electrophysiology."""

        if sett.Electrophysiology not in self._get_subsettings():
            raise ValueError("Expecting electrophysiology settings.")

    def update(self) -> None:
        """Update keyword database for electrophysiology."""
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

        # TODO: check if no existing nodeset ids conflict with surface ids
        # For now, new nodesets should be created after calling
        # self._update_nodesets_db()
        self._update_nodesets_db()
        self._update_parts_cellmodels()

        if self.model.conduction_mesh.number_of_cells != 0:
            # with smcoupl=1, mechanical coupling is disabled
            # with thcoupl=1, thermal coupling is disabled
            self.kw_database.ep_settings.append(keywords.EmControlCoupling(thcoupl=1, smcoupl=1))
            beam_pid = self._update_use_Purkinje()
        else:
            beam_pid = None

        # update ep settings
        self._update_ep_settings(beam_pid)
        self._update_stimulation()

        if self.model._add_blood_pool:
            self._update_blood_settings()

        if hasattr(self.model, "electrodes") and len(self.model.electrodes) != 0:
            self._update_ECG_coordinates()

        include_files = self._get_decknames_of_include()
        self.include_to_main(include_files)

        return

    def _update_dummy_material_db(self) -> None:
        """Add simple mechanics material for each defined part."""
        for part in self.model.parts:
            ep_mid = part.pid
            self.kw_database.material.append(
                keywords.MatElastic(mid=ep_mid, ro=1e-6, e=1),
            )

    def _update_ep_material_db(self) -> None:
        """Add electrophysiology material for each defined part."""
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
                LOGGER.info(f"Material of {part.name} is assigned automatically.")
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

    def _update_parts_cellmodels(self) -> None:
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

    def _create_myocardial_nodeset_layers(self) -> tuple[int, int, int]:
        """Create myocardial node set layers."""
        percent_endo = self.settings.electrophysiology.material.myocardium["percent_endo"]
        percent_mid = self.settings.electrophysiology.material.myocardium["percent_mid"]
        values = self.model.mesh.point_data["transmural"]
        # Values from experimental data. See:
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

    def _add_cell_model_keyword(self, matid: int, cellmodel: CellModel) -> None:
        """Add cell model keyword to the database."""
        if isinstance(cellmodel, CellModel.Tentusscher):
            self._add_Tentusscher_keyword(matid=matid, params=cellmodel.to_dictionary())
        else:
            raise NotImplementedError

    def _add_Tentusscher_keyword(self, matid: int, params: dict) -> None:  # noqa N802
        cell_kw = keywords.EmEpCellmodelTentusscher(**{**params})
        cell_kw.mid = matid
        # Note: bug in EmEpCellmodelTentusscher
        # the following 2 parameters cannot be assigned by above method
        cell_kw.gas_constant = 8314.472
        cell_kw.faraday_constant = 96485.3415

        self.kw_database.cell_models.append(cell_kw)

    def _update_ep_settings(self, beam_pid: list[int]) -> None:
        """Add the settings for the electrophysiology solver."""
        save_part_ids = []
        for part in self.model.parts:
            save_part_ids.append(part.pid)
        if beam_pid is not None:
            save_part_ids.extend(beam_pid)
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
                "EP timestep > Mechanics timestep. Setting EP timestep to Mechanics timestep."
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

    def _update_stimulation(self) -> None:
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

    def _add_stimulation_keyword(
        self, stim: Stimulation
    ) -> tuple[keywords.SetNodeList, Union[custom_keywords.EmEpTentusscherStimulus, str]]:
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
            # TODO: get the right part set ID
            # setpart_kwds = self.kw_database.ep_settings.get_kwds_by_type()
            # ID of the eikonal solver (different eikonal solves
            # can be performed in different parts of the model)
            eikonal_id = 1
            psid = 1
            eikonal_stim_content += f"{eikonal_id:>10d}{psid:>10d}{nsid:>10d}"
            if solvertype == "ReactionEikonal":
                eikonal_stim_content += "\n$ footType     footT     footA  footTauf   footVth\n"
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

        1/2 apex points for the left/bi-ventricle model.

        Sinoatrial node for four-chamber or full-heart model.

        Returns
        -------
        list[int]
            List of 0-based node IDs to stimulate.
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

            if ConductionPathType.SAN_AVN in [beam.name for beam in self.model.conduction_paths]:
                # Active SA node (belong to both solid and beam)
                stim_nodes = list(self.model.mesh.find_closest_point(LandMarks.SA_NODE.xyz, n=5))

                # add 1 more beam node to initiate wave propagation
                p = self.model.conduction_mesh.find_closest_point(LandMarks.SA_NODE.xyz, n=2)
                # take the second point, the first point is SA node itself
                pointid = self.model.conduction_mesh["_shifted_id"][p[1]]
                stim_nodes.append(pointid)

        # stimulate entire elements for Eikonal
        if self.settings.electrophysiology.analysis.solvertype in [
            "Eikonal",
            "ReactionEikonal",
        ]:
            stim_cells = np.where(np.isin(self.model.mesh.tetrahedrons, stim_nodes))[0]
            stim_nodes = np.unique(self.model.mesh.tetrahedrons[stim_cells].ravel())

        return stim_nodes

    def _update_blood_settings(self) -> None:
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

    def _update_ECG_coordinates(self) -> None:  # noqa N802
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

    def _update_solution_controls(self) -> None:
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

    def _update_main_db(self) -> None:
        pass

    def _update_use_Purkinje(self, associate_to_segment: bool = True) -> list[int]:  # noqa N802
        """Update keywords for Purkinje use."""
        if not isinstance(self.model, (FullHeart, FourChamber, BiVentricle, LeftVentricle)):
            LOGGER.error("Model type is not recognized.")
            return

        # write section
        sid = self.get_unique_section_id()
        self.kw_database.beam_networks.append(keywords.SectionBeam(secid=sid, elform=3, a=645))

        # id can be offset due to spring-type elements in mechanical
        beam_elem_id_offset = self.id_offset["element"]["discrete"]

        # write beam nodes
        new_nodes = self.model.conduction_mesh.points[
            (np.where(np.logical_not(self.model.conduction_mesh["_is-connected"]))[0])
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

        # loop for each beam
        beam_pid = []
        registered_surfaces = [surf for part in self.model.parts for surf in part.surfaces]
        for beam in self.model.conduction_paths:
            if isinstance(beam.ep_material, EPMaterial.DummyMaterial):
                epmat = self._get_default_beam_ep_material()
            else:
                epmat = beam.ep_material

            pid = self.get_unique_part_id()
            beam_pid.append(pid)

            name = beam.name.value

            result = next((x for x in registered_surfaces if x == beam.relying_surface), None)
            if result is not None:
                _node_set_id = result._seg_set_id
            else:
                _node_set_id = self._add_segment_from_surface(
                    beam.relying_surface, name + "_relying_segment"
                )

            # overwrite nsid if beam should not follow the motion of segment
            if not associate_to_segment:
                _node_set_id = -1

            # write Purkinje keyword
            self.kw_database.beam_networks.append(f"$$ {name} $$")
            origin_coordinates = beam.mesh.points[0]
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

            # write part
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

            # write material
            self.kw_database.beam_networks.append(keywords.MatNull(mid=pid, ro=1e-11))
            kw = self._get_ep_material_kw(pid, epmat)
            self.kw_database.beam_networks.append(kw)

            # write cell model
            self._add_cell_model_keyword(matid=pid, cellmodel=epmat.cell_model)

            # build element connectivity
            line_ids = self.model.conduction_mesh["_line-id"] == beam.id
            # get connectivity in _conduction_system
            edges_org = self.model.conduction_mesh.cells.reshape(-1, 3)[line_ids, 1:]
            # get shifted ID
            shift_id = self.model.conduction_mesh.point_data["_shifted_id"]
            edges = shift_id[edges_org]

            # write element
            beams_kw = add_beams_to_kw(
                beams=edges + 1,
                beam_kw=keywords.ElementBeam(),
                pid=pid,
                offset=beam_elem_id_offset,
            )
            self.kw_database.beam_networks.append(beams_kw)
            # offset beam element ID
            beam_elem_id_offset += len(edges)

        self.id_offset["element"]["discrete"] = beam_elem_id_offset

        return beam_pid

    def _get_default_beam_ep_material(self) -> EPMaterial.ActiveBeam:
        """Get default EP material for conduction beams."""
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

        return default_epmat

    def _add_segment_from_surface(self, surface: pv.PolyData, name: str) -> int:
        """Add a segment into keywords and return its segment ID.

        Parameters
        ----------
        surface : pv.PolyData
            surface
        name : str
            surface name

        Returns
        -------
        int
            segment ID
        """
        # surface_id = int(np.max(self.model.mesh.surface_ids) + 1)
        # self.model.mesh.add_surface(surface, surface_id, name=name)
        # self.model.mesh = self.model.mesh.clean()
        # surface = self.model.mesh.get_surface_by_name(name)
        seg_id = self.get_unique_segmentset_id()
        faces = surface.faces.reshape(-1, 4)[:, 1:]
        points = surface.points[faces]
        tree = spatial.cKDTree(self.model.mesh.points)
        a, b = tree.query(points)

        kw = create_segment_set_keyword(
            segments=b + 1,
            segid=seg_id,
            title=name,
        )
        # append this kw to the segment set database
        self.kw_database.segment_sets.append(kw)

        return seg_id

    def _update_export_controls(self) -> None:
        """Add solution controls to the main simulation."""
        self.kw_database.main.append(
            keywords.DatabaseBinaryD3Plot(dt=self.settings.electrophysiology.analysis.dt_d3plot.m)
        )

        return

    def _get_ep_material_kw(
        self, ep_mid: int, ep_material: EPMaterial
    ) -> Union[custom_keywords.EmMat001, custom_keywords.EmMat003]:
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
                # LS-DYNA bug prevents using isotropic mat (EMMAT001) for active isotropic case
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
    """Class for preparing the input for an electrophysiology LS-DYNA simulation with beams only."""

    def __init__(self, model: HeartModel, settings: SimulationSettings = None) -> None:
        super().__init__(model=model, settings=settings)
        self.kw_database = ElectrophysiologyDecks()
        """Collection of keywords relevant for electrophysiology."""

    def update(self) -> None:
        """Update keyword database for electrophysiology."""
        # self._isolate_atria_and_ventricles()

        ##
        self._update_main_db()
        self._update_solution_controls()
        self._update_export_controls()

        self._update_node_db()

        if self.model.conduction_mesh.number_of_cells != 0:
            # with smcoupl=1, coupling is disabled
            self.kw_database.ep_settings.append(keywords.EmControlCoupling(thcoupl=1, smcoupl=1))
            beam_pid = self._update_use_Purkinje(associate_to_segment=False)
        else:
            beam_pid = None

        # update ep settings
        self._update_ep_settings(beam_pid)
        self._update_stimulation()

        include_files = self._get_decknames_of_include()
        self.include_to_main(include_files)

        return
