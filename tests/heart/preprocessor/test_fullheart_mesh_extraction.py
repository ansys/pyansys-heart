"""Functional test to determine whether generated biventricle model has all the
expected features."""
import os
import shutil
import sys

import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.support import run_preprocessor
import pytest

from tests.heart.common import (
    compare_cavity_volume,
    compare_generated_mesh,
    compare_part_names,
    compare_surface_names,
)
from tests.heart.conftest import download_asset, get_assets_folder, get_workdir

# marks all tests with the 'requires_fluent' tag after this line
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
    model = run_preprocessor(
        model_type=models.FullHeart,
        database="Strocchi2020",
        path_original_mesh=path_to_case,
        work_directory=workdir,
        path_to_model=path_to_model,
        mesh_size=2,
    )

    yield

    # cleanup
    try:
        shutil.rmtree(workdir)
    except:
        print("Failed to cleanup.")


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
