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
"""Module containing classes for writing LS-DYNA keyword files for laplace simulations."""

import copy
from typing import Literal

import numpy as np
import pandas as pd
import pyvista as pv
import scipy.spatial as spatial

from ansys.dyna.core.keywords import keywords
from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.models import BiVentricle, FourChamber, FullHeart, HeartModel, LeftVentricle
from ansys.health.heart.objects import CapType, SurfaceMesh
from ansys.health.heart.writer.base_writer import BaseDynaWriter
from ansys.health.heart.writer.writer_utils import (
    create_element_solid_keyword,
    create_node_keyword,
    create_node_set_keyword,
)


class LaplaceWriter(BaseDynaWriter):
    """Class for preparing the input for a Laplace LS-DYNA simulation."""

    # constant nodeset ID for atrial valves/caps
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
    ) -> None:
        """Write thermal input to set up a Laplace Dirichlet problem.

        Parameters
        ----------
        model : HeartModel
            Heart model.
        type : Literal["uvc", "la_fiber", "ra_fiber", "D-RBM"]
            Simulation type.
        """
        super().__init__(model=model)
        self.type = type
        """Problem type."""
        self.landmarks = kwargs
        """Landmarks are ``laa``, ``raa``,  and ``top``."""
        self.target: pv.UnstructuredGrid = None
        """Target mesh related to the problem."""

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
                elems_to_keep.extend(model.left_ventricle.get_element_ids(model.mesh))
            else:
                elems_to_keep.extend(model.left_ventricle.get_element_ids(model.mesh))
                elems_to_keep.extend(model.right_ventricle.get_element_ids(model.mesh))
                elems_to_keep.extend(model.septum.get_element_ids(model.mesh))

            # model.mesh.clear_data()
            model.mesh["cell_ids"] = np.arange(0, model.mesh.n_cells, dtype=int)
            model.mesh["point_ids"] = np.arange(0, model.mesh.n_points, dtype=int)

            self.target = model.mesh.extract_cells(elems_to_keep)

        elif self.type == "la_fiber" or self.type == "ra_fiber":
            self._keep_parts(parts_to_keep)
            # model.mesh.clear_data()
            model.mesh["cell_ids"] = np.arange(0, model.mesh.n_cells, dtype=int)
            model.mesh["point_ids"] = np.arange(0, model.mesh.n_points, dtype=int)

            self.target = model.mesh.extract_cells(model.parts[0].get_element_ids(model.mesh))

    def _update_ra_top_nodeset(self, atrium: pv.UnstructuredGrid) -> None:
        """
        Define right atrium top nodeset with nodeset ID 10.

        Parameters
        ----------
        atrium : pv.UnstructuredGrid
            Right atrium PyVista object.
        """
        if "top" in self.landmarks.keys():
            top_ids = self._find_top_nodeset_by_geodesic(atrium)
        else:
            top_ids = self._find_top_nodeset_by_cut(atrium)

        # assign top nodeset
        kw = create_node_set_keyword(top_ids + 1, node_set_id=10, title="top")
        self.kw_database.node_sets.append(kw)

    def _find_top_nodeset_by_cut(self, atrium: pv.UnstructuredGrid) -> np.ndarray:
        """
        Define right atrium top nodeset.

        Cut through the center of TV, IVC, and SVC, expecting to result in
        three unconnected regions and the farthest is top.
        This method might fail with varying geometries. If so, the user
        must define the top landmarks.
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
            raise ValueError("Define top start and end points and then re-run.")

        # get tricuspid-valve name
        tv_name = CapType.TRICUSPID_VALVE_ATRIUM.value

        # compare closest point with TV nodes, top region should be far with TV nodeset
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

    def _find_top_nodeset_by_geodesic(self, atrium: pv.UnstructuredGrid) -> np.ndarray:
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

    def _define_ra_cut(self) -> tuple[np.ndarray, np.ndarray]:
        """Define a cutplane using the three caps of the right atrium."""
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

    def _update_ra_tricuspid_nodeset(self, atrium: pv.UnstructuredGrid) -> None:
        """Define the nodeset for the tricuspid wall and septum."""
        # get tricuspid-valve name
        tv_name = CapType.TRICUSPID_VALVE_ATRIUM.value

        # cut_normal is determined so that the first part is the septum and the second is free
        cut_center, cut_normal = self._define_ra_cut()

        # need a copied object to do clip, atrium is corrupted otherwise
        septum, free_wall = copy.deepcopy(atrium).clip(
            origin=cut_center, normal=cut_normal, crinkle=True, return_clipped=True
        )
        # IDs in full mesh
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

    def _update_atrial_caps_nodeset(self, atrium: pv.UnstructuredGrid) -> None:
        """Define nodesets for the caps."""
        for cap in self.model.parts[0].caps:
            # get node IDs for atrium mesh
            cap._mesh = self.model.mesh.get_surface(cap._mesh.id)
            ids_sub = np.where(np.isin(atrium["point_ids"], cap.global_node_ids_edge))[0]
            # create nodeset
            set_id = self._CAP_NODESET_MAP[cap.type]

            if set_id:  # Can be None for LEFT_ATRIUM_APPENDAGE
                kw = create_node_set_keyword(ids_sub + 1, node_set_id=set_id, title=cap.name)
                self.kw_database.node_sets.append(kw)

                # Add info to PyVista object, which is necessary for right atrial fibers.
                atrium[cap.type.value] = np.zeros(atrium.n_points, dtype=int)
                atrium[cap.type.value][ids_sub] = 1

        return

    def _update_la_bc(self) -> None:
        atrium = self.target

        def get_laa_nodes(atrium: pv.UnstructuredGrid, laa: np.ndarray) -> np.ndarray:
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
        endo_nodes = self._get_update_global_ids(self.model.left_atrium.endocardium.name)
        epi_nodes = self._get_update_global_ids(self.model.left_atrium.epicardium.name)
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

    def _update_ra_bc(self) -> None:
        atrium = self.target
        # caps
        self._update_atrial_caps_nodeset(atrium)

        # endo/epi
        endo_nodes = self._get_update_global_ids(self.model.right_atrium.endocardium.name)
        epi_nodes = self._get_update_global_ids(self.model.right_atrium.epicardium.name)
        epi_nodes = np.setdiff1d(epi_nodes, endo_nodes)

        self._add_nodeset(endo_nodes, "endocardium", nodeset_id=100)
        self._add_nodeset(epi_nodes, "epicardium", nodeset_id=200)

        # Find appendage apex
        tree = spatial.cKDTree(atrium.points)
        raa_ids = np.array(tree.query_ball_point(self.landmarks["raa"], self._LANDMARK_RADIUS))
        if len(raa_ids) == 0:
            LOGGER.error("No node is identified as right atrium appendage apex.")
            raise ValueError("No node is identified as right atrium appendage apex.")

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

    def update(self) -> None:
        """Update the keyword database."""
        # nodes
        node_kw = create_node_keyword(self.target.points)
        self.kw_database.nodes.append(node_kw)

        # part and materials
        self._update_parts_materials_db()

        # elems
        kw_elements = create_element_solid_keyword(
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

        include_files = self._get_decknames_of_include()
        self.include_to_main(include_files)

    def _get_update_global_ids(self, name: str) -> np.ndarray:
        """Get the update global IDs of a surface from its name."""
        # Note: This is temporary fix to make sure node IDs are correctly traced.
        surface1 = self.model.mesh.get_surface_by_name(name)
        return surface1.global_node_ids_triangles

    def _update_uvc_bc(self) -> None:
        # transmural uvc
        endo_nodes = self._get_update_global_ids(self.model.left_ventricle.endocardium.name)
        epi_nodes = self._get_update_global_ids(self.model.left_ventricle.epicardium.name)

        if not isinstance(self.model, LeftVentricle):
            rv_endo = self._get_update_global_ids(self.model.right_ventricle.endocardium.name)
            septum_endo = self._get_update_global_ids(self._get_rv_septum_endo_surface().name)
            rv_epi = self._get_update_global_ids(self.model.right_ventricle.epicardium.name)

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

    def _get_uvc_rotation_bc(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Select the nodeset on the long axis plane."""
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

    def _update_parts_materials_db(self) -> None:
        """Loop over parts defined in the model and create keywords."""
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

    def _update_main_db(self) -> None:
        self.kw_database.main.append(keywords.ControlSolution(soln=1))
        self.kw_database.main.append(keywords.ControlThermalSolver(atype=0, ptype=0, solver=11))
        self.kw_database.main.append(keywords.DatabaseBinaryD3Plot(dt=1.0))
        self.kw_database.main.append(keywords.DatabaseGlstat(dt=1.0))
        self.kw_database.main.append(keywords.DatabaseMatsum(dt=1.0))
        self.kw_database.main.append(keywords.DatabaseTprint(dt=1.0))
        self.kw_database.main.append(keywords.DatabaseExtentBinary(therm=2))  # save heat flux
        self.kw_database.main.append(keywords.ControlTermination(endtim=1, dtmin=1.0))

    def _add_nodeset(self, nodes: np.ndarray, title: str, nodeset_id: int = None) -> int:
        """Convert to local node ID and add to the nodeset.

        Parameters
        ----------
        nodes : np.ndarray
            Nodes global IDsx
        title : str
            Nodeset title.
        nodeset_id : int, default: None
            Attribute a nodeset ID if one is not given.

        Returns
        -------
        int
            Nodeset ID.
        """
        # get node IDs of submesh
        nodes = np.where(np.isin(self.target["point_ids"], nodes))[0]
        if nodeset_id is None:
            nodeset_id = self.get_unique_nodeset_id()
        # LS-DYNA ID starts with 1
        kw = create_node_set_keyword(nodes + 1, node_set_id=nodeset_id, title=title)
        self.kw_database.node_sets.append(kw)
        return nodeset_id

    def _update_drbm_bc(self) -> None:
        """Update D-RBM boundary conditions."""

        def clean_node_set(nodes: np.ndarray, exclude_nodes: np.ndarray = None) -> np.ndarray:
            """Ensure no duplicate or excluded nodes to avoid a thermal boundary condition error."""
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
        lv_endo_nodes = self._get_update_global_ids(self.model.left_ventricle.endocardium.name)
        lv_endo_nodes = clean_node_set(lv_endo_nodes, rings_nodes)
        # LV epi
        epi_nodes = self._get_update_global_ids(self.model.left_ventricle.epicardium.name)
        epi_nodes = clean_node_set(epi_nodes, np.hstack((lv_endo_nodes, rings_nodes)))
        # LV apex
        la_node = self.model.get_apex_node_set(part="left")

        if not isinstance(self.model, LeftVentricle):
            # Right ventricle endocardium
            septum_endo = self._get_rv_septum_endo_surface()
            rv_endo_nodes = np.hstack(
                (
                    self._get_update_global_ids(self.model.right_ventricle.endocardium.name),
                    self._get_update_global_ids(septum_endo.name),
                )
            )
            rv_endo_nodes = clean_node_set(rv_endo_nodes, rings_nodes)

            # append RV epi
            epi_nodes = np.hstack(
                (
                    epi_nodes,
                    self._get_update_global_ids(self.model.right_ventricle.epicardium.name),
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

    def _get_rv_septum_endo_surface(self) -> SurfaceMesh:
        """Get the right ventricle septal surface."""
        try:
            return self.model.right_ventricle.septum
        except Exception as e:
            raise ValueError(f"Septum endocardium surface is not found in right ventricle. {e}")

    def _update_ventricular_caps_nodes(
        self,
    ) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray], bool]:
        combined_av_mv = False  # combined mitral and aortic valve
        mv_nodes = av_nodes = tv_nodes = pv_nodes = None

        for cap in self.model.all_caps:
            cap._mesh = self.model.mesh.get_surface(cap._mesh.id)
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

    def add_case(
        self, case_id: int, case_name: str, set_ids: list[int], bc_values: list[float]
    ) -> None:
        """Add a case to the keyword database.

        Parameters
        ----------
        case_id : int
           Case ID.
        case_name : str
            Case name, which is the d3plot filename.
        set_ids : list[int]
            List of nodeset IDs for boundary conditions.
        bc_values : list[float]
            List of boundary condition values.
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
