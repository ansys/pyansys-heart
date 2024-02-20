# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
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

"""Functional test to determine whether generated biventricle model has all the
expected features."""

import os
import shutil

import ansys.heart.preprocessor.models.v0_1.models as models
from ansys.heart.simulator.support import run_preprocessor
import numpy as np
import pytest

from .common import (
    _deprecated_compare_caps_nodeids,
    _deprecated_compare_caps_num_nodeids,
    _deprecated_compare_cavity_topology,
    _deprecated_compare_surface_faces,
    compare_cavity_volume,
    compare_part_element_ids,
    compare_part_names,
    compare_surface_names,
)
from .conftest import download_asset, get_assets_folder, get_workdir


# run this fixture first
@pytest.fixture(autouse=True, scope="module")
def extract_bi_ventricle():
    """Extract BiVentricular model which is similar to the reference model.

    Note
    ----
    Do this once as fixture.
    """

    assets_folder = get_assets_folder()
    path_to_case = os.path.join(assets_folder, "cases", "strocchi2020", "01", "01.case")

    if not os.path.isfile(path_to_case):
        path_to_case = download_asset("Strocchi2020", casenumber=1)

    # load model to test against
    path_to_reference_model = os.path.join(
        assets_folder,
        "reference_models",
        "strocchi2020",
        "01",
        "BiVentricle",
        "heart_model.pickle",
    )
    assert os.path.isfile(path_to_case)
    assert os.path.isfile(path_to_reference_model)

    global reference_model
    reference_model = models.HeartModel.load_model(path_to_reference_model)

    assert isinstance(reference_model, models.BiVentricle), (
        "Reference model should be of type %s" % models.BiVentricle.__class__.__name__
    )

    workdir = os.path.join(get_workdir(), reference_model.__class__.__name__)
    path_to_model = os.path.join(workdir, "heart_model.pickle")

    global model
    model = run_preprocessor(
        model_type=reference_model.__class__,
        database="Strocchi2020",
        path_original_mesh=path_to_case,
        work_directory=workdir,
        path_to_model=path_to_model,
        mesh_size=reference_model.info.mesh_size,
        add_blood_pool=True,
    )

    yield

    # cleanup
    shutil.rmtree(workdir)


@pytest.mark.skip(reason="Not re-implemented yet.")
def test_blood_pools():
    """Test if blood fluid_mesh a mesh."""
    assert isinstance(model.fluid_mesh.tetrahedrons, np.ndarray)
    assert isinstance(model.fluid_mesh.nodes, np.ndarray)
    # may be too strict
    assert np.shape(model.fluid_mesh.tetrahedrons)[0] > 450000
    assert np.shape(model.fluid_mesh.tetrahedrons)[1] == 4
    assert np.shape(model.fluid_mesh.nodes)[0] > 120000
    assert np.shape(model.fluid_mesh.nodes)[1] == 3

    return


@pytest.mark.skip(reason="Not re-implemented yet.")
def test_part_names():
    compare_part_names(model, reference_model)
    pass


@pytest.mark.skip(reason="Not re-implemented yet.")
def test_part_element_ids():
    compare_part_element_ids(model, reference_model)
    pass


@pytest.mark.skip(reason="Not re-implemented yet.")
def test_surface_names():
    compare_surface_names(model, reference_model)
    pass


@pytest.mark.skip(reason="Not re-implemented yet.")
def test_surface_faces():
    _deprecated_compare_surface_faces(model, reference_model)
    pass


@pytest.mark.skip(reason="Not re-implemented yet.")
def test_cavities_topology():
    _deprecated_compare_cavity_topology(model, reference_model)
    pass


def test_cavities_volumes():
    compare_cavity_volume(model, reference_model)
    pass


@pytest.mark.skip(reason="Not re-implemented yet.")
def test_caps_nodeids():
    _deprecated_compare_caps_nodeids(model, reference_model)
    pass


@pytest.mark.skip(reason="Not re-implemented yet.")
def test_caps_num_nodeids():
    _deprecated_compare_caps_num_nodeids(model, reference_model)
    pass
