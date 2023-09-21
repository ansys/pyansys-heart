"""Compares generated k-files of mechanics against a reference."""
import pytest

pytestmark = pytest.mark.local

import os

import ansys.heart.preprocessor.models as models
import ansys.heart.writer.dynawriter as writers
import pytest

from tests.heart.conftest import compare_string_with_file, get_assets_folder


# NOTE: Using a pickled reference model won't work when we make
# changes to the underlying objects. Pickle won't be able to reconstruct
# the object.
@pytest.fixture(autouse=True, scope="module")
def initialize_model():
    """Load a model that can be converted to LS-DYNA input."""
    pickle_file = os.path.join(
        get_assets_folder(),
        "reference_models",
        "strocchi2020",
        "01",
        "BiVentricle",
        "heart_model.pickle",
    )
    model = models.HeartModel.load_model(pickle_file)
    assert isinstance(model, models.BiVentricle), "Expecting a BiVentricle model"

    global WRITER_MECHANICS
    WRITER_MECHANICS = writers.MechanicsDynaWriter(model)
    WRITER_MECHANICS.update()
    return


@pytest.mark.xfail(reason="Settings in main subject to change")
def test_main():
    """Compare contents of main to a reference file."""
    ref_folder = os.path.join(get_assets_folder(), "k_files", "mechanics")
    ref_file = os.path.join(ref_folder, "main.k")

    string_to_test = WRITER_MECHANICS.kw_database.main.write()

    compare_string_with_file(string_to_test, ref_file)


@pytest.mark.xfail(
    reason="Uses pickled model. Pickled load will fail if we make changes to the object itself."
)
def test_parts():
    """Compare contents of main to a reference file."""
    ref_folder = os.path.join(get_assets_folder(), "k_files", "mechanics")
    ref_file = os.path.join(ref_folder, "parts.k")

    string_to_test = WRITER_MECHANICS.kw_database.parts.write()

    compare_string_with_file(string_to_test, ref_file)


@pytest.mark.xfail(reason="Material parameters subject to change")
def test_material():
    """Compare contents of main to a reference file."""
    ref_folder = os.path.join(get_assets_folder(), "k_files", "mechanics")
    ref_file = os.path.join(ref_folder, "material.k")

    string_to_test = WRITER_MECHANICS.kw_database.material.write()

    compare_string_with_file(string_to_test, ref_file)
