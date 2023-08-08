"""Functional test to determine whether generated biventricle model has all the
expected features."""
import os
import shutil
import sys

import ansys.heart.preprocessor.models_new as models
from ansys.heart.simulator.support import (
    get_input_geom_and_part_defintions_from_public_database,
    preprocess_model,
)
import pytest

from .common import (
    compare_cavity_volume,
    compare_generated_mesh,
    compare_part_names,
    compare_surface_names,
)
from .conftest import download_asset, get_assets_folder, get_workdir


# run this fixture first
@pytest.fixture(autouse=True, scope="module")
def extract_fullheart():
    """Extract FullHeart model which is similar to the reference model.

    Note
    ----
    Do this once as fixture.
    """

    assets_folder = get_assets_folder()
    path_to_case = os.path.join(assets_folder, "cases", "strocchi2020", "01", "01.case")

    if not os.path.isfile(path_to_case):
        path_to_case = download_asset("Strocchi2020", casenumber=1)

    path_to_reference_stats = path_to_reference_stats = os.path.join(
        assets_folder,
        "reference_models",
        "strocchi2020",
        "01",
        "FullHeart",
        "stats_reference_model_FullHeart.json",
    )
    assert os.path.isfile(path_to_case)
    assert os.path.isfile(path_to_reference_stats)

    global ref_stats
    import json

    with open(path_to_reference_stats, "r") as openfile:
        ref_stats = json.load(openfile)

    workdir = os.path.join(get_workdir(), "FullHeart")
    path_to_model = os.path.join(workdir, "heart_model.pickle")

    global model

    input_geom, part_definitions = get_input_geom_and_part_defintions_from_public_database(
        path_to_case, model_type="", database="Strocchi2020"
    )

    info = models.ModelInfo(
        input=input_geom,
        scalar="surface-id",
        part_definitions=part_definitions,
        work_directory=workdir,
        mesh_size=1.5,
    )

    model = preprocess_model(
        info=info, model_type="FullHeart", clean_workdir=False, use_wrapper=True
    )

    yield

    # cleanup
    shutil.rmtree(workdir)


def test_part_names():
    compare_part_names(model, ref_stats)
    pass


def test_surface_names():
    compare_surface_names(model, ref_stats)
    pass


def test_cavities_volumes():
    compare_cavity_volume(model, ref_stats)
    pass


@pytest.mark.xfail(
    sys.platform == "linux", reason="Mesh generation slightly different for Linux version of Fluent"
)
def test_mesh():
    """Test the number of tetrahedrons and triangles in the volume mesh and surface meshes"""
    compare_generated_mesh(model, ref_stats)
