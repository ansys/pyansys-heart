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

"""
Full-heart conduction system example
------------------------------------
This example demonstrates how to build a conduction system for a full-heart model,
including the Purkinje network and major conduction paths.

The Purkinje network is generated using a fractal tree algorithm, which is a
third-party tool (`fractal_tree`). It creates a
realistic, branching Purkinje system on the endocardial surfaces of the ventricles.

In addition to the Purkinje fibers, this example defines other key conduction pathways,
such as the sinoatrial (SA) node to atrioventricular (AV) node paths, the Bachmann bundle,
the His bundle and its bifurcation, as well as the left and right bundle branches.
Each conduction path is constructed using anatomical landmarks and keypoints, and
integrated into the heart model for simulation.
"""

###############################################################################
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from pathlib import Path

import numpy as np
import pyvista as pv

try:
    from fractal_tree import FractalTreeParameters, Mesh, generate_fractal_tree

    NO_FRACTAL_TREE = False

except ImportError:
    print(
        "The 'fractal_tree' package is not a requirement of this package, "
        "you can install it using 'pip install fractal_tree'\n"
        "Otherwise, the precomputed results will be used."
    )
    NO_FRACTAL_TREE = True

from ansys.health.heart.examples import (
    get_fractal_tree_purkinje,
    get_preprocessed_fullheart,
)
import ansys.health.heart.models as models
from ansys.health.heart.models_utils import HeartModelUtils
from ansys.health.heart.pre.conduction_path import ConductionPath, ConductionPathType

# Set the working directory and path to the model.
workdir = Path.home() / "pyansys-heart" / "downloads" / "Rodero2021" / "01" / "FullHeart"
path_to_model, path_to_partinfo, _ = get_preprocessed_fullheart(resolution="2.0mm")

###############################################################################
# Load the full-heart model
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Load the full-heart model.
model: models.FullHeart = models.FullHeart(working_directory=workdir)
model.load_model_from_mesh(path_to_model, path_to_partinfo)


def run_fractal_tree(
    surface: pv.PolyData, dir: np.ndarray, apex_node: np.ndarray
) -> tuple[pv.PolyData, pv.PolyData]:
    """
    Generate a fractal tree (Purkinje network) on a given surface.

    Parameters
    ----------
    surface : pv.PolyData
        The endocardial surface to grow the tree on.
    dir : np.ndarray
        Initial direction vector for tree growth.
    apex_node : np.ndarray
        Apex node coordinates.

    Returns
    -------
    pv.PolyData
        The generated Purkinje network as a PolyData object.
    """
    verts = surface.points
    connectivity = surface.faces.reshape(-1, 4)[:, 1:]
    index = np.linalg.norm(np.subtract(verts, apex_node), axis=1).argmin()

    param = FractalTreeParameters(
        N_it=10,
        length=8,
        initial_direction=dir,
        l_segment=1,
        fascicles_length=(10, 10),
        repulsitivity=0.2,
        save=False,
    )
    # Generate the fractal tree
    tree = generate_fractal_tree(
        Mesh(verts=verts, connectivity=connectivity, init_node=verts[index, :]),
        param,
    )

    # Convert the tree to a PolyData object
    purkinje = pv.PolyData(tree.nodes, lines=np.insert(tree.lines, 0, 2, axis=1))
    # Add the 'end' point data to indicate leaf nodes.
    p1 = np.unique(np.array(tree.lines)[:, 0])
    p2 = np.unique(np.array(tree.lines)[:, 1])
    leaf_nodes = np.setdiff1d(p2, p1)
    purkinje.point_data["end"] = np.zeros(purkinje.number_of_points)
    purkinje.point_data["end"][leaf_nodes] = 1

    return purkinje


def generate_purkinje():
    """
    Generate left and right Purkinje networks using a fractal tree algorithm.

    Returns
    -------
    tuple
        (left_purkinje, right_purkinje) as pyvista.PolyData objects with
        'end' point data.
    """
    np.random.seed(1234)  # For reproducibility of the fractal tree generation

    # Left
    lv_endo = model.left_ventricle.endocardium
    lv_apex = model.left_ventricle.apex_points[0].xyz
    base_node = model.left_ventricle.caps[0].centroid

    dir = -base_node + lv_apex
    dir = dir / np.linalg.norm(dir)
    l_purkinje = run_fractal_tree(lv_endo, dir, lv_apex)

    # Right
    rv_endo = model.right_ventricle.endocardium
    rv_apex = model.right_ventricle.apex_points[0].xyz
    r_purkinje = run_fractal_tree(rv_endo, dir, rv_apex)

    return l_purkinje, r_purkinje


def build_conudciton_system(model: models.HeartModel):
    """Build the conduction system for the heart model, including 3
    sa-av paths, Bachmann bundle, and the His bundle with its bifurcation,
    left and right bundle branches, and the Purkinje fibers

    Parameters
    ----------
    model : models.HeartModel
        Heart model instance to which the conduction system will be added.
    """
    if NO_FRACTAL_TREE:
        # Use precomputed Purkinje
        left, right = get_fractal_tree_purkinje()
    else:
        # Generate the Purkinje network for both ventricles.
        left, right = generate_purkinje()

    # Create conduction paths for the Purkinje network
    left_pirkinje = ConductionPath(
        ConductionPathType.LEFT_PURKINJE,
        left,
        is_connected=np.zeros(left.n_points),
        id=1,
        relying_surface=model.left_ventricle.endocardium,
    )
    # create purkinje junctions (PMJ) on leaf nodes
    left_pirkinje.add_pmj_path(pmj_list=np.where(left["end"] == 1)[0])

    # Create conduction paths for the right Purkinje network
    right_pirkinje = ConductionPath(
        ConductionPathType.RIGHT_PURKINJE,
        right,
        is_connected=np.zeros(right.n_points),
        id=2,
        relying_surface=model.right_ventricle.endocardium,
    )
    right_pirkinje.add_pmj_path(pmj_list=np.where(right["end"] == 1)[0])

    # Define the conduction paths for the heart's conduction system
    sa = HeartModelUtils.define_sino_atrial_node(model, target_coord=[6, 66, 88])
    av = HeartModelUtils.define_atrio_ventricular_node(model)

    # Create the SA-AV node conduction path
    sa_av = ConductionPath.create_from_keypoints(
        name=ConductionPathType.SAN_AVN,
        keypoints=[sa.xyz, av.xyz],
        id=3,
        base_mesh=model.right_atrium.endocardium,
        connection="first",
        line_length=None,
        center=True,
    )
    sa_av.add_pmj_path(list(range(1, sa_av.mesh.n_points - 1, 4)))

    # Create the His bundle and its bifurcation
    his_bif = HeartModelUtils.define_his_bundle_bifurcation_node(model)
    his_left_point = HeartModelUtils.define_his_bundle_end_node(model, side="left")
    his_right_point = HeartModelUtils.define_his_bundle_end_node(model, side="right")
    his_top = ConductionPath.create_from_keypoints(
        name=ConductionPathType.HIS_TOP,
        keypoints=[av.xyz, his_bif.xyz],
        id=4,
        base_mesh=model.mesh,
    )
    his_top.up_path = sa_av

    # TODO : we can set a delay by his_top.ep_material

    his_left = ConductionPath.create_from_keypoints(
        name=ConductionPathType.HIS_LEFT,
        keypoints=[his_bif.xyz, his_left_point.xyz],
        id=5,
        base_mesh=model.mesh,
    )
    his_left.up_path = his_top

    his_right = ConductionPath.create_from_keypoints(
        name=ConductionPathType.HIS_RIGHT,
        keypoints=[his_bif.xyz, his_right_point.xyz],
        id=6,
        base_mesh=model.mesh,
    )
    his_right.up_path = his_top

    # Create the left and right bundle branches
    left_bundle = ConductionPath.create_from_keypoints(
        name=ConductionPathType.LEFT_BUNDLE_BRANCH,
        keypoints=[his_left_point.xyz, model.left_ventricle.apex_points[0].xyz],
        id=7,
        base_mesh=model.left_ventricle.endocardium,
        line_length=None,
        center=True,
    )
    # Add PMJ path for left bundle branch
    pmj_list = list(
        range(
            int((0.4 * left_bundle.mesh.n_points)),
            int((0.9 * left_bundle.mesh.n_points)),
            4,  # every 4 nodes
        )
    )
    left_bundle.add_pmj_path(pmj_list, merge_with="node")

    left_bundle.up_path = his_left
    left_bundle.down_path = left_pirkinje

    surface_ids = [
        model.right_ventricle.endocardium.id,
        model.right_ventricle.septum.id,
    ]
    endo_surface = model.mesh.get_surface(surface_ids)
    right_bundle = ConductionPath.create_from_keypoints(
        name=ConductionPathType.RIGHT_BUNDLE_BRANCH,
        keypoints=[his_right_point.xyz, model.right_ventricle.apex_points[0].xyz],
        id=8,
        base_mesh=endo_surface,
        line_length=None,
    )
    right_bundle.up_path = his_right
    right_bundle.down_path = right_pirkinje

    # Create the Bachmann bundle
    surface_ids = [
        model.left_atrium.epicardium.id,
        model.right_atrium.epicardium.id,
    ]
    surface = model.mesh.get_surface(surface_ids)
    bachman_bundle = ConductionPath.create_from_keypoints(
        name=ConductionPathType.BACHMANN_BUNDLE,
        keypoints=[sa.xyz, [46, 102, 97]],
        id=9,
        base_mesh=surface,
        line_length=None,
        center=True,
    )
    bachman_bundle.add_pmj_path(list(range(1, bachman_bundle.mesh.n_points - 1, 4)))
    bachman_bundle.up_path = sa_av

    # Create the mid and post SA-AV node conduction paths
    mid_sa_av = ConductionPath.create_from_keypoints(
        name=ConductionPathType.MID_SAN_AVN,
        keypoints=[
            sa.xyz,
            [10, 79, 64],
            [18, 95, 41],
            [32, 93, 31],
            [43, 88, 26],
            av.xyz,
        ],
        id=10,
        base_mesh=model.right_atrium.endocardium,
        line_length=None,
        center=True,
    )
    mid_sa_av.add_pmj_path(list(range(5, mid_sa_av.mesh.n_points - 5, 4)))
    mid_sa_av.up_path = sa_av
    mid_sa_av.down_path = sa_av

    post_sa_av = ConductionPath.create_from_keypoints(
        name=ConductionPathType.POST_SAN_AVN,
        keypoints=[sa.xyz, [2, 65, 53], [6, 73, 34], [25, 75, 26], av.xyz],
        id=11,
        base_mesh=model.right_atrium.endocardium,
        line_length=None,
        center=True,
    )
    post_sa_av.add_pmj_path(list(range(5, post_sa_av.mesh.n_points - 5, 4)))
    post_sa_av.up_path = sa_av
    post_sa_av.down_path = sa_av

    # Assign conduction paths to the model
    model.assign_conduction_paths(
        [
            left_pirkinje,
            right_pirkinje,
            sa_av,
            his_top,
            his_left,
            his_right,
            left_bundle,
            right_bundle,
            bachman_bundle,
            mid_sa_av,
            post_sa_av,
        ]
    )

    return


# Generate the conduction system
build_conudciton_system(model)

###############################################################################
# Plot the conduction system
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model.plot_purkinje()
