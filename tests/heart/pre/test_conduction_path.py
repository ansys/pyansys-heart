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

"""Test for reading purkinje network as a beam mesh. Uses a mock Mesh object."""

import os

import numpy as np
import pyvista as pv

from ansys.health.heart.models_utils import HeartModelUtils
from ansys.health.heart.pre.conduction_path import ConductionPath, ConductionPathType
from ansys.health.heart.settings.material.ep_material import EPMaterial
from tests.heart.conftest import get_assets_folder, get_fourchamber, get_fullheart


def test_compute_sa_node():
    fourchamber = get_fourchamber()

    sa_node = HeartModelUtils.define_sino_atrial_node(fourchamber)

    assert np.allclose(sa_node.xyz, np.array([-48.80218005, 107.90170883, 423.33688959]))
    assert sa_node.node_id == 105021


def test_compute_av_node():
    fourchamber = get_fourchamber()
    av_node = HeartModelUtils.define_atrio_ventricular_node(fourchamber)
    assert np.allclose(av_node.xyz, np.array([-10.16353107, 108.95410155, 371.9505145]))
    assert av_node.node_id == 100501


def test_compute_his_bif_node():
    fourchamber = get_fourchamber()
    bif_node = HeartModelUtils.define_his_bundle_bifurcation_node(
        fourchamber, target_coord=np.array([-10.16353107, 108.95410155, 371.9505145])
    )
    assert np.allclose(bif_node.xyz, np.array([1.22510233, 110.31896126, 364.402475]))
    assert bif_node.node_id == 25326


def test_compute_his_end_node():
    fourchamber = get_fourchamber()
    # need pre steps
    HeartModelUtils.define_atrio_ventricular_node(fourchamber)
    HeartModelUtils.define_his_bundle_bifurcation_node(fourchamber)
    # test for left
    left = HeartModelUtils.define_his_bundle_end_node(fourchamber, side="left")

    assert np.allclose(left.xyz, np.array([4.15421613, 113.63743565, 369.27104019]))
    assert left.node_id == 49464
    # test for right
    right = HeartModelUtils.define_his_bundle_end_node(fourchamber, side="right")

    assert np.allclose(right.xyz, np.array([2.93215687, 106.09459183, 365.20590901]))
    assert right.node_id == 43585


def test_create_conductionbeams_on_surface():
    """Test conductionbeams can be initialized correctly on a surface."""
    model = get_fourchamber()
    sa = HeartModelUtils.define_sino_atrial_node(model)
    av = HeartModelUtils.define_atrio_ventricular_node(model)

    sa_av = ConductionPath.create_from_keypoints(
        name=ConductionPathType.SAN_AVN,
        keypoints=[sa.xyz, av.xyz],
        id=2,
        base_mesh=model.right_atrium.endocardium,
        connection="none",
    )

    assert sa_av.mesh.n_lines == 48
    assert np.isclose(sa_av.length, 64.36592438345)


def test_create_conductionbeams_on_surface_with_refinement():
    """Test conductionbeams can be initialized correctly on a surface."""
    model = get_fourchamber()
    sa = HeartModelUtils.define_sino_atrial_node(model)
    av = HeartModelUtils.define_atrio_ventricular_node(model)

    sa_av = ConductionPath.create_from_keypoints(
        name=ConductionPathType.SAN_AVN,
        keypoints=[sa.xyz, av.xyz],
        id=2,
        base_mesh=model.right_atrium.endocardium,
        connection="all",
        line_length=0.5,
    )
    sa_av.relying_surface.save("surf.vtp")
    sa_av.mesh.save("line.vtp")

    assert sa_av.mesh.n_lines == 140
    assert np.isclose(sa_av.length, 64.36592438345)
    # check only necessary points are connected
    assert np.sum(sa_av.is_connected) == 49


def test_create_conductionbeams_in_solid():
    model = get_fourchamber()
    av = HeartModelUtils.define_atrio_ventricular_node(model)
    bif = HeartModelUtils.define_his_bundle_bifurcation_node(model)
    his_top = ConductionPath.create_from_keypoints(
        name=ConductionPathType.HIS_TOP,
        keypoints=[av.xyz, bif.xyz],
        id=1,
        base_mesh=model.mesh,
        connection="none",
    )
    assert np.isclose(his_top.length, 14.276232139149878)
    assert his_top.relying_surface.n_cells == 9


def meshes_equal(mesh1: pv.DataSet, mesh2: pv.DataSet) -> bool:
    if isinstance(mesh1, pv.PolyData):
        mesh1 = mesh1.cast_to_unstructured_grid()
    if isinstance(mesh2, pv.PolyData):
        mesh2 = mesh2.cast_to_unstructured_grid()

    if not np.array_equal(mesh1.points, mesh2.points):
        return False
    if not np.array_equal(mesh1.cells, mesh2.cells):
        return False
    if mesh1.point_data.keys() != mesh2.point_data.keys():
        return False
    for key in mesh1.point_data.keys():
        if not np.array_equal(mesh1.point_data[key], mesh2.point_data[key]):
            return False
    if mesh1.cell_data.keys() != mesh2.cell_data.keys():
        return False
    for key in mesh1.cell_data.keys():
        if not np.array_equal(mesh1.cell_data[key], mesh2.cell_data[key]):
            return False
    return True


def test_conduction():
    model = get_fullheart()
    folder = os.path.join(
        get_assets_folder(), "reference_models", "strocchi2020", "01", "conduction"
    )
    # f1 = os.path.join(folder, "purkinjeNetwork_001.k")
    # left_purkjinje = model.add_purkinje_from_kfile(f1, _ConductionType.LEFT_PURKINJE.value)
    # ref0 = pv.read(os.path.join(folder, "left_purkinje.vtp"))
    # assert meshes_equal(ref0, left_purkjinje)

    # new method
    beam_list = HeartModelUtils.define_default_conduction_system(model, purkinje_folder=folder)
    model.assign_conduction_paths(beam_list)
    res = model.conduction_mesh

    # old method
    # from ansys.health.heart.pre.conduction_beam import _compute_heart_conductionsystem

    # f1 = os.path.join(folder, "purkinjeNetwork_001.k")
    # f2 = os.path.join(folder, "purkinjeNetwork_002.k")
    # model.add_purkinje_from_kfile(f1, _ConductionType.LEFT_PURKINJE.value)
    # model.add_purkinje_from_kfile(f2, _ConductionType.RIGHT_PURKINJE.value)
    # _compute_heart_conductionsystem(model, 1.5)

    ref = pv.read(os.path.join(folder, "conduction.vtu"))

    assert res.n_cells == ref.n_cells
    assert res.n_points == ref.n_points
    assert np.allclose(res.points, ref.points, atol=1e-3)
    # old method has HIS together but new method split it into 3 part
    # assert np.array_equal(res["_line-id"], ref["_line-id"])
    assert np.array_equal(res["_is-connected"], ref["_is-connected"])

    # test ID shift with merging to solid
    assert np.sum(model._shifted_id()) == 587209415


def test_conductionbeams_from_k():
    """Test conductionbeams can be initialized correctly from a k file."""
    model = get_fullheart()
    folder = os.path.join(
        get_assets_folder(), "reference_models", "strocchi2020", "01", "conduction"
    )
    f1 = os.path.join(folder, "purkinjeNetwork_001.k")

    l_pj = ConductionPath.create_from_k_file(
        name=ConductionPathType.LEFT_PURKINJE,
        k_file=f1,
        id=1,
        base_mesh=model.left_ventricle.endocardium,
        model=model,
    )
    assert l_pj.name == ConductionPathType.LEFT_PURKINJE
    assert l_pj.ep_material == EPMaterial.DummyMaterial()
    ref0 = pv.read(os.path.join(folder, "left_purkinje.vtp"))
    assert meshes_equal(l_pj.mesh, ref0)
