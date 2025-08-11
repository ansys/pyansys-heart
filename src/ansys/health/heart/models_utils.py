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

"""Stateless methods for the heart model."""

from dataclasses import dataclass
import os
from typing import Literal

import numpy as np
import pyvista as pv

from ansys.health.heart import LOG as LOGGER
import ansys.health.heart.models as models
from ansys.health.heart.objects import CapType, Point
from ansys.health.heart.pre.conduction_path import ConductionPath, ConductionPathType


@dataclass
class LandMarks:
    """Heart anatomical points."""

    SA_NODE = Point("SA_node", xyz=None, node_id=None)
    AV_NODE = Point("AV_node", xyz=None, node_id=None)
    HIS_BIF_NODE = Point("His_bifurcation", xyz=None, node_id=None)
    HIS_LEFT_END_NODE = Point("His_left_end", xyz=None, node_id=None)
    HIS_RIGHT_END_NODE = Point("His_right_end", xyz=None, node_id=None)
    BACHMAN_END_NODE = Point("Bachman_end", xyz=None, node_id=None)
    LEFT_FASCILE_END = Point("Left_fasicle_end", xyz=None, node_id=None)

    # NOTE bring APEX here
    LEFT_APEX = Point("Left_apex", xyz=None, node_id=None)
    RIGHT_APEX = Point("Right_apex", xyz=None, node_id=None)


class HeartModelUtils:
    """Stateless methods for the heart model."""

    @staticmethod
    def define_sino_atrial_node(
        model: models.FullHeart | models.FourChamber, target_coord: np.ndarray | list = None
    ) -> LandMarks | None:
        """Define Sino-atrial node.

        Parameters
        ----------
        model : models.FullHeart | models.FourChamber
            Heart model.
        target_coord : np.ndarray | list, default: None
            If ``None``, the target coordinate is computed as the midpoint between
            the centroids of the superior and inferior vena cavae. If a coordinate is provided,
            the closest point on the right atrium endocardium surface to that coordinate is used.

        Returns
        -------
        LandMarks | None
            SA node.
        """
        try:
            right_atrium_endo = model.mesh.get_surface(model.right_atrium.endocardium.id)
        except AttributeError:
            LOGGER.error("Cannot find right atrium to create SinoAtrial node")
            return

        if target_coord is None:
            sup_vcava_centroid = next(
                cap.centroid
                for cap in model.right_atrium.caps
                if cap.type == CapType.SUPERIOR_VENA_CAVA
            )
            inf_vcava_centroid = next(
                cap.centroid
                for cap in model.right_atrium.caps
                if cap.type == CapType.INFERIOR_VENA_CAVA
            )

            # define SinoAtrial node:
            target_coord = sup_vcava_centroid - (inf_vcava_centroid - sup_vcava_centroid) / 2

        target_id = pv.PolyData(
            model.mesh.points[right_atrium_endo.global_node_ids_triangles, :]
        ).find_closest_point(target_coord)

        sino_atrial_node_id = right_atrium_endo.global_node_ids_triangles[target_id]

        LandMarks.SA_NODE.xyz = model.mesh.points[sino_atrial_node_id, :]
        LandMarks.SA_NODE.node_id = sino_atrial_node_id

        return LandMarks.SA_NODE

    @staticmethod
    def define_atrio_ventricular_node(
        model: models.FullHeart | models.FourChamber, target_coord: np.ndarray | list = None
    ) -> LandMarks | None:
        """Define Atrio-ventricular node.

        Parameters
        ----------
        model : models.FullHeart | models.FourChamber
            Heart model.
        target_coord : np.ndarray | list, default: None
            If ``None``, the target coordinate is computed as the closest point on the right atrium
            endocardium surface to the right ventricle septum. If a coordinate is provided, the
            closest point on the right atrium endocardium surface to that coordinate is used.

        Returns
        -------
        LandMarks | None
            AV node.
        """
        try:
            right_atrium_endo = model.mesh.get_surface(model.right_atrium.endocardium.id)
        except AttributeError:
            LOGGER.error("Cannot find right atrium to create SinoAtrial node")
            return

        if target_coord is None:
            right_septum = model.mesh.get_surface(model.right_ventricle.septum.id)
            # define AtrioVentricular as the closest point to septum
            target_id = pv.PolyData(
                model.mesh.points[right_atrium_endo.global_node_ids_triangles, :]
            ).find_closest_point(right_septum.center)

        else:
            target_id = pv.PolyData(
                model.mesh.points[right_atrium_endo.global_node_ids_triangles, :]
            ).find_closest_point(target_coord)

        # assign a point
        av_id = right_atrium_endo.global_node_ids_triangles[target_id]
        LandMarks.AV_NODE.xyz = model.mesh.points[av_id, :]
        LandMarks.AV_NODE.node_id = av_id

        return LandMarks.AV_NODE

    @staticmethod
    def define_his_bundle_bifurcation_node(
        model: models.FourChamber | models.FullHeart, target_coord: np.ndarray | list = None
    ) -> LandMarks | None:
        """Define His bundle bifurcation node.

        Parameters
        ----------
        model : models.FourChamber | models.FullHeart
            Heart model.
        target_coord : np.ndarray | list, default: None
            If ``None``, the target coordinate is computed as the closest point in the septum to
            the AV node. If a coordinate is provided, the closest point in the septum to that
            coordinate is used.

        Returns
        -------
        LandMarks | None
            HIS bifurcation node.
        """
        if target_coord is None:
            av_coord = LandMarks.AV_NODE.xyz
            if av_coord is None:
                LOGGER.error("AV node need to be defined before.")
                return
            target_coord = av_coord

        septum_point_ids = np.unique(
            np.ravel(model.mesh.tetrahedrons[model.septum.get_element_ids(model.mesh)])
        )

        # remove nodes on surface, to make sure His bundle nodes are inside of septum
        septum_point_ids = np.setdiff1d(
            septum_point_ids,
            model.mesh.get_surface(model.left_ventricle.endocardium.id).global_node_ids_triangles,
        )
        septum_point_ids = np.setdiff1d(
            septum_point_ids,
            model.mesh.get_surface(model.right_ventricle.septum.id).global_node_ids_triangles,
        )

        septum_pointcloud = pv.PolyData(model.mesh.points[septum_point_ids, :])

        # Define start point: closest to artria
        pointcloud_id = septum_pointcloud.find_closest_point(target_coord)

        bifurcation_id = septum_point_ids[pointcloud_id]
        LandMarks.HIS_BIF_NODE.xyz = model.mesh.points[bifurcation_id, :]
        LandMarks.HIS_BIF_NODE.node_id = bifurcation_id

        return LandMarks.HIS_BIF_NODE

    @staticmethod
    def define_his_bundle_end_node(
        model: models.FullHeart | models.FourChamber,
        target_coord: np.ndarray | list = None,
        side: Literal["left", "right"] = "left",
        n_close: int = 20,
    ) -> LandMarks | None:
        """Define His bundle end node.

        Parameters
        ----------
        model : models.FullHeart | models.FourChamber
            Heart model.
        target_coord : np.ndarray | list, default: None
            If ``None``, the target coordinate is computed as the n-th closest point
            on the endocardium to the His bundle bifurcation node.
            Not implemented yet if a coordinate is provided.
        side : Literal[&quot;left&quot;, &quot;right&quot;], default: "left"
            Side of the heart to define the end node for.
        n_close : int, default: 20
            n-th closest point to the bifurcation node, to avoid too close to the bifurcation node.

        Returns
        -------
        LandMarks | None
            End node of His left or right bundle.
        """
        if side == "left":
            endo = model.mesh.get_surface(model.left_ventricle.endocardium.id)
        elif side == "right":
            endo = model.mesh.get_surface(model.right_ventricle.septum.id)

        if target_coord is not None:
            LOGGER.error("Do not support user defined point.")
            return
        else:
            # find n-th closest point to bifurcation
            bifurcation_coord = LandMarks.HIS_BIF_NODE.xyz
            if bifurcation_coord is None:
                LOGGER.error("AV node need to be defined before.")
                return
            temp_id = pv.PolyData(
                model.mesh.points[endo.global_node_ids_triangles, :]
            ).find_closest_point(bifurcation_coord, n=n_close)[n_close - 1]

            his_end_id = endo.global_node_ids_triangles[temp_id]

        if side == "left":
            LandMarks.HIS_LEFT_END_NODE.node_id = his_end_id
            LandMarks.HIS_LEFT_END_NODE.xyz = model.mesh.points[his_end_id, :]
            return LandMarks.HIS_LEFT_END_NODE

        elif side == "right":
            LandMarks.HIS_RIGHT_END_NODE.node_id = his_end_id
            LandMarks.HIS_RIGHT_END_NODE.xyz = model.mesh.points[his_end_id, :]

            return LandMarks.HIS_RIGHT_END_NODE

    @staticmethod
    def define_bachman_bundle_end_node(
        model: models.FullHeart | models.FourChamber, target_coord=None
    ) -> LandMarks | None:
        """Define Bachmann bundle end node."""
        NotImplementedError

    @staticmethod
    def define_fascile_bundle_end_node(
        model: models.FullHeart | models.FourChamber, target_coord=None
    ) -> LandMarks | None:
        """Define fascile bundle end node."""
        NotImplementedError

    @staticmethod
    def define_full_conduction_system(
        model: models.FullHeart | models.FourChamber, purkinje_folder: str
    ) -> list[ConductionPath]:
        """Define the full conduction system.

        Parameters
        ----------
        model : models.FullHeart | models.FourChamber
            Heart model.
        purkinje_folder : str
            Folder with LS-DYNA's Purkinje generation.

        Returns
        -------
        list[ConductionPath]
            List of conduction paths.
        """
        left_purkinje = ConductionPath.create_from_k_file(
            ConductionPathType.LEFT_PURKINJE,
            k_file=os.path.join(purkinje_folder, "purkinjeNetwork_001.k"),
            id=1,
            base_mesh=model.left_ventricle.endocardium,
            model=model,
        )

        right_purkinje = ConductionPath.create_from_k_file(
            ConductionPathType.RIGHT_PURKINJE,
            k_file=os.path.join(purkinje_folder, "purkinjeNetwork_002.k"),
            id=2,
            base_mesh=model.right_ventricle.endocardium,
            model=model,
        )

        sa = HeartModelUtils.define_sino_atrial_node(model)
        av = HeartModelUtils.define_atrio_ventricular_node(model)

        sa_av = ConductionPath.create_from_keypoints(
            name=ConductionPathType.SAN_AVN,
            keypoints=[sa.xyz, av.xyz],
            id=3,
            base_mesh=model.right_atrium.endocardium,
            line_length=None,
            connection="first",
        )

        his_bif = HeartModelUtils.define_his_bundle_bifurcation_node(model)
        his_left_point = HeartModelUtils.define_his_bundle_end_node(model, side="left")
        his_right_point = HeartModelUtils.define_his_bundle_end_node(model, side="right")

        his_top = ConductionPath.create_from_keypoints(
            name=ConductionPathType.HIS_TOP,
            keypoints=[av.xyz, his_bif.xyz],
            id=4,
            base_mesh=model.mesh.extract_cells_by_type(10),
        )
        his_top.up_path = sa_av

        his_left = ConductionPath.create_from_keypoints(
            name=ConductionPathType.HIS_LEFT,
            keypoints=[his_bif.xyz, his_left_point.xyz],
            id=5,
            base_mesh=model.mesh.extract_cells_by_type(10),
        )
        his_left.up_path = his_top

        his_right = ConductionPath.create_from_keypoints(
            name=ConductionPathType.HIS_RIGHT,
            keypoints=[his_bif.xyz, his_right_point.xyz],
            id=6,
            base_mesh=model.mesh.extract_cells_by_type(10),
        )
        his_right.up_path = his_top

        left_bundle = ConductionPath.create_from_keypoints(
            name=ConductionPathType.LEFT_BUNDLE_BRANCH,
            keypoints=[his_left_point.xyz, model.left_ventricle.apex_points[0].xyz],
            id=7,
            base_mesh=model.left_ventricle.endocardium,
            line_length=None,
            center=True,
        )

        # Create Purkinje junctions on the lower part of the left bundle branch,
        # so depolarization begins in the left side of the septum.
        pmj_list = list(
            range(
                int((0.4 * left_bundle.mesh.n_points)),
                int((0.9 * left_bundle.mesh.n_points)),
                4,  # every 4 nodes
            )
        )
        left_bundle.add_pmj_path(pmj_list, merge_with="cell")
        left_bundle.up_path = his_left
        left_bundle.down_path = left_purkinje

        surface_ids = [model.right_ventricle.endocardium.id, model.right_ventricle.septum.id]
        endo_surface = model.mesh.get_surface(surface_ids)

        right_bundle = ConductionPath.create_from_keypoints(
            name=ConductionPathType.RIGHT_BUNDLE_BRANCH,
            keypoints=[his_right_point.xyz, model.right_ventricle.apex_points[0].xyz],
            id=8,
            base_mesh=endo_surface,
            line_length=None,
        )
        right_bundle.up_path = his_right
        right_bundle.down_path = right_purkinje
        return [
            left_purkinje,
            right_purkinje,
            sa_av,
            his_top,
            his_left,
            his_right,
            left_bundle,
            right_bundle,
        ]
