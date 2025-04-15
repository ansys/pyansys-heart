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

import os

import numpy as np
import pyvista as pv

from ansys.health.heart.models_utils import HeartModelUtils
from ansys.health.heart.objects import _ConductionType
from ansys.health.heart.pre.conduction_beam2 import ConductionBeams, ConductionBeamType
from ansys.health.heart.settings.material.ep_material import EPMaterial
from tests.heart.conftest import get_assets_folder, get_fullheart


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
    model.add_conduction_beam(beam_list)
    res = model._conduction_system

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
    assert np.array_equal(res["_line-id"], ref["_line-id"])
    assert np.array_equal(res["_is-connected"], ref["_is-connected"])


def test_conductionbeams_init():
    model = get_fullheart()
    folder = os.path.join(
        get_assets_folder(), "reference_models", "strocchi2020", "01", "conduction"
    )
    f1 = os.path.join(folder, "purkinjeNetwork_001.k")
    left_purkjinje = model.add_purkinje_from_kfile(f1, _ConductionType.LEFT_PURKINJE.value)

    l_pj = ConductionBeams(
        name=ConductionBeamType.LEFT_PURKINJE,
        mesh=left_purkjinje,
        id=1,
        is_connected=left_purkjinje["_is-connected"],
        relying_surface=model.left_ventricle.endocardium,
    )

    assert l_pj.name == ConductionBeamType.LEFT_PURKINJE
    assert l_pj.ep_material == EPMaterial.DummyMaterial()
    ref0 = pv.read(os.path.join(folder, "left_purkinje.vtp"))
    assert meshes_equal(l_pj.mesh, ref0)
