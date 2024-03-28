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

import glob
import os
import pathlib
import shutil

os.environ["ANSYS_HEART_MODEL_VERSION"] = "v0.2"

import ansys.heart.preprocessor.models.v0_2.models as models
from ansys.heart.writer.dynawriter import (
    ElectrophysiologyDynaWriter,
    FiberGenerationDynaWriter,
    MechanicsDynaWriter,
    PurkinjeGenerationDynaWriter,
    ZeroPressureMechanicsDynaWriter,
)
import pytest

from tests.heart.common import (
    compare_cavity_volume,
    compare_generated_mesh,
    compare_part_names,
    compare_surface_names,
)
from tests.heart.conftest import get_assets_folder, get_workdir
from tests.heart.end2end.compare_k import read_file

# marks all tests with the 'requires_fluent' tag after this line
# pytest.mark.requires_fluent
pytestmark = pytest.mark.requires_fluent


# run this fixture first
@pytest.fixture(autouse=True, scope="module")
def extract_fullheart():
    """Extract FullHeart model which is similar to the reference model.

    Note
    ----
    Do this once as fixture.
    """
    assets_folder = get_assets_folder()

    path_to_input_vtp = os.path.join(
        assets_folder, "cases", "Strocchi2020", "01", "input_polydata.vtp"
    )
    part_def_file = os.path.join(
        assets_folder, "cases", "Strocchi2020", "01", "part_definitions_fullheart.json"
    )

    path_to_reference_stats = os.path.join(
        assets_folder,
        "reference_models",
        "strocchi2020",
        "01",
        "FullHeart",
        "stats_reference_model_FullHeart.json",
    )

    assert os.path.isfile(path_to_input_vtp)
    assert os.path.isfile(path_to_reference_stats)

    global ref_stats
    import json

    with open(path_to_reference_stats, "r") as openfile:
        ref_stats = json.load(openfile)

    workdir = os.path.join(get_workdir(), "v02", "FullHeart")

    with open(part_def_file, "r") as f:
        part_definitions = json.load(f)

    global model
    info = models.ModelInfo(
        path_to_input_vtp,
        "boundary-id",
        part_definitions=part_definitions,
        work_directory=workdir,
        mesh_size=2.0,
    )
    model = models.FullHeart(info)
    model.load_input()
    model.mesh_volume(use_wrapper=True)
    model._update_parts()

    yield

    # cleanup
    try:
        shutil.rmtree(workdir)
        print()
    except:
        print("Failed to cleanup.")


@pytest.mark.models_v2
def test_part_names():
    compare_part_names(model, ref_stats)
    pass


@pytest.mark.models_v2
def test_surface_names():
    compare_surface_names(model, ref_stats)
    pass


@pytest.mark.models_v2
def test_cavities_volumes():
    compare_cavity_volume(model, ref_stats)
    pass


@pytest.mark.models_v2
@pytest.mark.xfail(reason="Comparing against v0.1 model - uses different meshing method(s)")
def test_mesh():
    """Test the number of tetrahedrons and triangles in the volume mesh and surface meshes"""
    compare_generated_mesh(model, ref_stats)


@pytest.mark.models_v2
@pytest.mark.skip(reason="This test is currently too strict - rewrite or disable")
@pytest.mark.parametrize(
    "writer_class",
    [
        ElectrophysiologyDynaWriter,
        MechanicsDynaWriter,
        ZeroPressureMechanicsDynaWriter,
        FiberGenerationDynaWriter,
        PurkinjeGenerationDynaWriter,
    ],
)
def test_writers(writer_class):
    """Test whether all writers yield the same .k files as the reference model.

    Notes
    -----
    This skips over most .k files that contain mesh related info.
    """
    writer = writer_class(model)
    ref_folder = os.path.join(
        get_assets_folder(),
        "reference_models",
        "strocchi2020",
        "01",
        "FullHeart",
        "k_files",
        writer_class.__name__,
    )
    to_test_folder = os.path.join(get_workdir(), "fullheart", writer_class.__name__)
    writer.update()
    writer.export(to_test_folder)

    ref_files = glob.glob(os.path.join(ref_folder, "*.k"))
    # compare each of the reference files to the files that were generated.
    for ref_file in ref_files:
        file_to_compare = os.path.join(to_test_folder, pathlib.Path(ref_file).name)
        assert read_file(ref_file) == read_file(
            file_to_compare
        ), f"File {pathlib.Path(ref_file).name} does not match."

    return
