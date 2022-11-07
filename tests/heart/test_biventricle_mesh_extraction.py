import os
import shutil
import sys

import pytest

if not sys.platform.startswith("win"):
    pytest.skip("Skipping windows-only tests", allow_module_level=True)

import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.support import run_preprocessor
from common import compare_caps, compare_cavities, compare_parts, compare_surfaces
from conftest import get_assets_folder, get_workdir


# run this fixture first
@pytest.fixture(autouse=True, scope="module")
def extract_bi_ventricle():
    """Extracts BiVentricular model which is similar to the reference model.
    Do this once as fixture
    """

    assets_folder = get_assets_folder()
    path_to_case = os.path.join(assets_folder, "cases", "strocchi2020", "01", "01.case")

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
    )

    yield

    # cleanup
    shutil.rmtree(workdir)


def test_parts():
    compare_parts(model, reference_model)

    pass


def test_surfaces():
    compare_surfaces(model, reference_model)
    pass


def test_cavities():
    compare_cavities(model, reference_model)
    pass


def test_caps():
    compare_caps(model, reference_model)
    pass
