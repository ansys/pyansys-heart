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

"""Functional test to determine whether generated biventricle and fullheart
model has all the expected features and writes consistent .k files."""

import copy
import glob
import json
import os
import pathlib
import tempfile
from typing import Union

import pytest
import yaml

import ansys.heart.preprocessor.models as models
import ansys.heart.writer.dynawriter as writers
from tests.heart.common import compare_stats_mesh, compare_stats_names, compare_stats_volumes
from tests.heart.conftest import get_assets_folder
from tests.heart.end2end.compare_k import read_file

#! Note: should run fast tests before slow tests.


# get the input files from the assets directory.
# given the specified model type.
def _get_inputs(model_type: Union[models.BiVentricle, models.FullHeart]):
    """Get the input files from the assets directory given the specified model type."""
    model_path = os.path.join(
        get_assets_folder(),
        "reference_models",
        "strocchi2020",
        "01",
    )
    input_path = os.path.join(model_path, "inputs")

    input_polydata = os.path.join(input_path, "input_polydata.vtp")

    if model_type == models.BiVentricle:
        part_definitions_file = os.path.join(input_path, "part_definitions_biventricle.json")
        mesh_file = os.path.join(model_path, "biventricle-post-mesh.vtu")
        stats_file = os.path.join(model_path, "stats_biventricle_v0-2.yml")

    elif model_type == models.FullHeart:
        part_definitions_file = os.path.join(input_path, "part_definitions_fullheart.json")
        mesh_file = os.path.join(model_path, "fullheart-post-mesh.vtu")
        stats_file = os.path.join(model_path, "stats_fullheart_v0-2.yml")

    with open(stats_file, "r") as f:
        ref_stats = yaml.load(f, yaml.SafeLoader)

    with open(part_definitions_file, "r") as f:
        part_definitions = json.load(f)

    return (input_polydata, part_definitions, ref_stats, mesh_file)


# Set up testing parameter combinations
# first item: model type, second item flag indicating whether to try to remesh
# the model.
test_params_fast = [
    pytest.param([models.BiVentricle, False], marks=pytest.mark.extract_models),
    pytest.param([models.FullHeart, False], marks=pytest.mark.extract_models),
]
test_params_fast_ids = ["BiVentricle-NoRemesh", "FullHeart-NoRemesh"]

test_params_slow = [
    pytest.param(
        [models.BiVentricle, True],
        marks=[pytest.mark.extract_models, pytest.mark.requires_fluent],
    ),
    pytest.param(
        [models.FullHeart, True],
        marks=[pytest.mark.extract_models, pytest.mark.requires_fluent],
    ),
]
test_params_slow_ids = [
    "BiVentricle-Remesh",
    "FullHeart-Remesh",
]


# default params are the fast parameters.
@pytest.fixture(
    scope="module",
    params=test_params_fast,
    ids=test_params_fast_ids,
)
def extract_model(request):

    if list(request.param):
        model_type = request.param[0]
        mesh_volume = request.param[1]
    else:
        raise TypeError("Expecting list of 2 input parameters.")

    if not model_type in (models.BiVentricle, models.FullHeart):
        raise TypeError(f"Expecting model to be of types: {models.BiVentricle, models.FullHeart}")
        return

    inputs = _get_inputs(model_type)  # model_type)
    ref_stats = inputs[2]
    mesh_file = inputs[3]

    # global workdir
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:

        info = models.ModelInfo(
            input=inputs[0],
            scalar="boundary-id",
            part_definitions=inputs[1],
            mesh_size=2.0,
            work_directory=workdir,
        )
        model = model_type(info)
        if not isinstance(model, (models.BiVentricle, models.FullHeart)):
            exit()

        model.load_input()
        # model.mesh_volume(wrapper=True) # could use this: but requires fluent
        if mesh_volume:
            model.mesh_volume(use_wrapper=True)
        else:
            model.mesh.load_mesh(mesh_file)
        model._update_parts()

        yield model, ref_stats

    return


@pytest.mark.parametrize(
    "extract_model",
    test_params_slow + test_params_fast,
    ids=test_params_slow_ids + test_params_fast_ids,
    indirect=["extract_model"],
)
def test_names(extract_model):
    """Test if relevant features are present in model."""
    model, ref_stats = extract_model
    stats = model.summary()
    compare_stats_names(stats, ref_stats)
    pass


@pytest.mark.parametrize(
    "extract_model",
    test_params_slow + test_params_fast,
    ids=test_params_slow_ids + test_params_fast_ids,
    indirect=["extract_model"],
)
def test_cavity_volumes(extract_model):
    """Test consistency of cavity volumes."""
    model, ref_stats = extract_model
    stats = model.summary()
    compare_stats_volumes(stats, ref_stats)
    pass


@pytest.mark.parametrize(
    "extract_model",
    test_params_slow + test_params_fast,
    ids=test_params_slow_ids + test_params_fast_ids,
    indirect=["extract_model"],
)
@pytest.mark.xfail(
    reason="Different Fluent versions or os's may yield slightly different meshing results"
)
def test_mesh_stats(extract_model):
    model, ref_stats = extract_model
    stats = model.summary()
    compare_stats_mesh(stats, ref_stats)
    pass


@pytest.fixture(autouse=True, scope="module")
def unpack_k_files():
    """Unpacks the .k files of a specific model if necessary."""
    import zipfile

    zip_file = os.path.join(
        get_assets_folder(),
        "reference_models",
        "strocchi2020",
        "01",
        "k_files_biventricle_fullheart.zip",
    )

    with zipfile.ZipFile(zip_file, "r") as zip_ref:
        zip_ref.extractall(os.path.dirname(zip_file))

    yield

    # cleanup
    try:
        import shutil

        shutil.rmtree(os.path.join(os.path.dirname(zip_file), "_BiVentricle"))
        shutil.rmtree(os.path.join(os.path.dirname(zip_file), "_FullHeart"))
    except:
        pass

    return


@pytest.mark.parametrize(
    "writer_class",
    [
        writers.ElectrophysiologyDynaWriter,
        writers.ElectroMechanicsDynaWriter,
        writers.MechanicsDynaWriter,
        writers.ZeroPressureMechanicsDynaWriter,
        writers.FiberGenerationDynaWriter,
        writers.PurkinjeGenerationDynaWriter,
    ],
)
@pytest.mark.k_file_writer
@pytest.mark.xfail(
    reason="""Testing .k files is mesh sensitive and subject to changes in model configuration.
    If no changes to the model are expected than this test should pass"""
)
def test_writers(extract_model, writer_class):
    """Test whether all writers yield the same .k files as the reference model.

    Notes
    -----
    This skips over most .k files that contain mesh related info.
    """
    model, _ = extract_model
    writer = writer_class(copy.deepcopy(model))

    if isinstance(model, models.BiVentricle):
        ref_folder = os.path.join(
            get_assets_folder(),
            "reference_models",
            "strocchi2020",
            "01",
            "_BiVentricle",
            "k_files1",
            writer_class.__name__,
        )
    elif isinstance(model, models.FullHeart):
        ref_folder = os.path.join(
            get_assets_folder(),
            "reference_models",
            "strocchi2020",
            "01",
            "_FullHeart",
            "k_files1",
            writer_class.__name__,
        )

    # with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:
        to_test_folder = os.path.join(workdir, writer_class.__name__)
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


@pytest.mark.parametrize(
    "writer_class",
    [
        writers.ElectrophysiologyDynaWriter,
        writers.ElectroMechanicsDynaWriter,
        writers.MechanicsDynaWriter,
        writers.ZeroPressureMechanicsDynaWriter,
        writers.FiberGenerationDynaWriter,
        writers.PurkinjeGenerationDynaWriter,
    ],
)
@pytest.mark.k_file_writer
@pytest.mark.xfail(
    reason="""Testing .k files is mesh sensitive and subject to changes in model configuration.
    If no changes to the model are expected than this test should pass"""
)
def test_writers_after_load_model(extract_model, writer_class):
    """Test whether all writers yield the same .k files as the reference model.

    Notes
    -----
    This tests the .k files after saving and loading a model.
    """
    model, _ = extract_model

    if isinstance(model, models.BiVentricle):
        ref_folder = os.path.join(
            get_assets_folder(),
            "reference_models",
            "strocchi2020",
            "01",
            "_BiVentricle",
            "k_files1",
            writer_class.__name__,
        )
    elif isinstance(model, models.FullHeart):
        ref_folder = os.path.join(
            get_assets_folder(),
            "reference_models",
            "strocchi2020",
            "01",
            "_FullHeart",
            "k_files1",
            writer_class.__name__,
        )

    # with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:

        model_path = os.path.join(workdir, model.__class__.__name__ + ".vtu")
        partinfo = model_path.replace(".vtu", ".partinfo.json")

        model.save_model(model_path)

        model1 = type(model)(models.ModelInfo(work_directory=workdir))
        model1.load_model_from_mesh(model_path, partinfo)
        model1._extract_apex()
        model1.compute_left_ventricle_anatomy_axis()
        model1.compute_left_ventricle_aha17()

        writer = writer_class(copy.deepcopy(model1))

        to_test_folder = os.path.join(workdir, writer_class.__name__)
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
