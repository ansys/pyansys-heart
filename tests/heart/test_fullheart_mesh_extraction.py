"""Functional test to determine whether generated fullheart model has all the
expected features."""
import os
import shutil

import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.support import run_preprocessor
import pytest

from .common import (
    _deprecated_compare_caps_nodeids,
    _deprecated_compare_caps_num_nodeids,
    _deprecated_compare_cavity_topology,
    compare_cavity_volume,
    compare_part_element_ids,
    compare_part_names,
    _deprecated_compare_surface_faces,
    compare_surface_names,
)
from .conftest import download_asset, get_assets_folder, get_workdir


@pytest.fixture(autouse=True, scope="module")
def extract_full_heart():
    """Extracts FullHeart model which is similar to the reference model.
    Do this once as fixture
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
        "FullHeart",
        "heart_model.pickle",
    )

    global reference_model
    reference_model = models.HeartModel.load_model(path_to_reference_model)

    assert isinstance(reference_model, models.FullHeart), (
        "Reference model should be of type %s" % models.FullHeart.__class__.__name__
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
    )

    yield

    # cleanup
    shutil.rmtree(workdir)


def test_part_names():
    compare_part_names(model, reference_model)
    pass


@pytest.mark.xfail
def test_part_element_ids():
    compare_part_element_ids(model, reference_model)
    pass


def test_surface_names():
    compare_surface_names(model, reference_model)
    pass


@pytest.mark.xfail
def test_surface_faces():
    _deprecated_compare_surface_faces(model, reference_model)
    pass


@pytest.mark.xfail
def test_cavities_topology():
    _deprecated_compare_cavity_topology(model, reference_model)
    pass


def test_cavities_volumes():
    compare_cavity_volume(model, reference_model)
    pass


@pytest.mark.xfail
def test_caps_nodeids():
    _deprecated_compare_caps_nodeids(model, reference_model)
    pass


def test_caps_num_nodeids():
    _deprecated_compare_caps_num_nodeids(model, reference_model)
    pass
