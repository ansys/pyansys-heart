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
import copy
import glob
import os
import pathlib
import shutil
import tempfile

from ansys.heart.preprocessor.helpers import model_summary
import ansys.heart.preprocessor.models.v0_1.models as models
from ansys.heart.simulator.support import run_preprocessor
from ansys.heart.writer.dynawriter import (
    ElectrophysiologyDynaWriter,
    FiberGenerationDynaWriter,
    MechanicsDynaWriter,
    PurkinjeGenerationDynaWriter,
    ZeroPressureMechanicsDynaWriter,
)
import pytest
import yaml

pytestmark = pytest.mark.skip(reason="v0.1 HeartModel removed")
pytestmark = pytest.mark.requires_fluent

from tests.heart.common import compare_stats_mesh, compare_stats_names, compare_stats_volumes
from tests.heart.conftest import download_asset, get_assets_folder, get_workdir
from tests.heart.end2end.compare_k import read_file


# run this fixture first
@pytest.fixture(autouse=True, scope="module")
def extract_bi_ventricle():
    """Extract BiVentricular model which is similar to the reference model.

    Notes
    -----
    Do this once as fixture.
    """

    assets_folder = get_assets_folder()
    path_to_case = download_asset("Strocchi2020", casenumber=1, clean_folder=True)

    path_to_reference_stats = os.path.join(
        get_assets_folder(),
        "reference_models",
        "strocchi2020",
        "01",
        "stats_biventricle_v0-1.yml",
    )

    assert os.path.isfile(path_to_case)
    assert os.path.isfile(path_to_reference_stats)

    global stats_ref
    with open(path_to_reference_stats, "r") as f:
        stats_ref = yaml.load(f, yaml.SafeLoader)

    # Use a temp directory for generating the model.
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tmpdirname:
        workdir = str(pathlib.Path(tmpdirname))
        path_to_model = os.path.join(workdir, "heart_model.pickle")

        # copy mesh file
        shutil.copyfile(
            os.path.join(
                assets_folder,
                "reference_models",
                "strocchi2020",
                "01",
                "BiVentricle",
                "fluent_volume_mesh.msh.h5",
            ),
            os.path.join(workdir, "fluent_volume_mesh.msh.h5"),
        )

        global model
        model = run_preprocessor(
            model_type=models.BiVentricle,
            database="Strocchi2020",
            path_original_mesh=path_to_case,
            work_directory=workdir,
            path_to_model=path_to_model,
            mesh_size=2,
            skip_meshing=True,
        )

        global stats
        stats = model_summary(model)

    return


@pytest.mark.models_v1
def test_names():
    """Test if relevant features are present in model."""
    compare_stats_names(stats, stats_ref)
    pass


@pytest.mark.models_v1
def test_cavity_volumes():
    """Test consistency of cavity volumes."""
    compare_stats_volumes(stats, stats_ref)
    pass


@pytest.mark.models_v1
@pytest.mark.xfail(reason="Different Fluent versions may yield slightly different meshing results")
def test_mesh_stats():
    compare_stats_mesh(stats, stats_ref)
    pass


# functional tests to determine whether any change was made to
# LS-DYNA input.
@pytest.mark.models_v1
@pytest.mark.xfail(reason="This test is currently too strict - rewrite or disable")
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
    writer = writer_class(copy.deepcopy(model))
    ref_folder = os.path.join(
        get_assets_folder(),
        "reference_models",
        "strocchi2020",
        "01",
        "BiVentricle",
        "k_files",
        writer_class.__name__,
    )
    to_test_folder = os.path.join(get_workdir(), "biventricle", writer_class.__name__)
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
