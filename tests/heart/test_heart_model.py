"""Tests shigh-level heart model class. """
import glob as glob
import json
import os

import ansys.heart.preprocessor.models as models
from conftest import create_directory, get_workdir
import numpy as np
import pytest

if os.name != "nt":
    skip_test = True
else:
    skip_test = False


def _get_test_model_info() -> models.ModelInfo:
    """Get a test model info and populates it."""
    info = models.ModelInfo(
        database="Strocchi2020",
        work_directory=get_workdir(),
        path_to_case="path-to-case",
        path_to_simulation_mesh="path-to-simulation-mesh",
        mesh_size=2.0,
    )
    return info


def _get_test_model(model_type: type) -> models.HeartModel:
    """Get a test heart model of specific type."""
    info = _get_test_model_info()
    model = model_type(info)
    return model


@pytest.fixture(autouse=True)
def cleanup_working_directory_after_tests():
    """Remove any files that were added during the test."""
    # execute before tests
    list_of_files_before = glob.glob(os.path.join(get_workdir(), "*"))
    yield
    # execute after tests
    list_of_files_after = glob.glob(os.path.join(get_workdir(), "*"))
    files_to_remove = list(set(list_of_files_after) - set(list_of_files_before))
    for file in files_to_remove:
        os.remove(file)
    return


@pytest.mark.skipif(skip_test, reason="Requires Windows")
def test_model_info_dump():
    """Test dumping of model info to json."""
    info = _get_test_model_info()
    create_directory(get_workdir())
    path_to_model_info = os.path.join(get_workdir(), "model_info.json")

    info.dump_info(path_to_model_info)

    path_to_model_info
    with open(path_to_model_info) as json_file:
        json_data = json.load(json_file)

    assert info.database == json_data["_database"], "Database not the same"
    assert info.workdir == json_data["workdir"], "Workdir not the same"
    assert (
        info.path_to_original_mesh == json_data["path_to_original_mesh"]
    ), "Path to original mesh not the same"
    assert (
        info.path_to_simulation_mesh == json_data["path_to_simulation_mesh"]
    ), "Path to simulation mesh not the same"
    assert info.mesh_size == json_data["mesh_size"], "Mesh size not the same"

    pass


@pytest.mark.skipif(skip_test, reason="Requires Windows")
def test_dump_model_001():
    """Test dumping of model to disk: using path in ModelInfo."""
    info = _get_test_model_info()
    info.workdir = get_workdir()
    create_directory(info.workdir)
    model = models.BiVentricle(info)
    expected_path = os.path.join(model.info.workdir, "heart_model.pickle")
    if os.path.isfile(expected_path):
        os.remove(expected_path)

    model.info.path_to_model = None
    model.dump_model()
    assert os.path.isfile(expected_path)


@pytest.mark.skipif(skip_test, reason="Requires Windows")
def test_dump_model_002():
    """Test dumping of model to disk: using specific path in ModelInfo."""
    info = _get_test_model_info()
    info.workdir = get_workdir()
    create_directory(info.workdir)
    model = models.BiVentricle(info)
    expected_path = os.path.join(model.info.workdir, "heart_model1.pickle")
    if os.path.isfile(expected_path):
        os.remove(expected_path)

    model.info.path_to_model = expected_path
    model.dump_model()
    assert os.path.isfile(expected_path)


@pytest.mark.skipif(skip_test, reason="Requires Windows")
def test_dump_model_003():
    """Test dumping of model to disk: using specific path."""
    info = _get_test_model_info()
    info.workdir = get_workdir()
    create_directory(info.workdir)
    model = models.BiVentricle(info)
    model.info.path_to_model = os.path.join(model.info.workdir, "heart_model2.pickle")
    expected_path = os.path.join(model.info.workdir, "heart_model3.pickle")
    if os.path.isfile(expected_path):
        os.remove(expected_path)

    model.dump_model(expected_path)
    assert os.path.isfile(expected_path)


@pytest.mark.skipif(skip_test, reason="Requires Windows")
def test_model_load():
    """Test loading model from pickle."""
    model: models.BiVentricle = _get_test_model(models.BiVentricle)
    # populate model
    model.left_ventricle.element_ids = np.array([1, 2, 3, 4], dtype=int)
    model.right_ventricle.element_ids = np.array([11, 66, 77, 88], dtype=int)

    model.left_ventricle.endocardium.faces = np.array([[1, 2, 3], [1, 2, 4]], dtype=int)
    model.right_ventricle.endocardium.faces = np.array([[11, 22, 33], [11, 22, 44]], dtype=int)

    model.mesh.tetrahedrons = np.array([[1, 2, 3, 4], [1, 2, 3, 5]], dtype=int)
    model.mesh.nodes = np.array([[0.0, 0.0, 0.1], [1.0, 1.0, 1.1]], dtype=float)

    # dump model to disk
    path_to_heart_model = os.path.join(get_workdir(), "heart_model.pickle")
    model.dump_model(path_to_heart_model)

    assert os.path.isfile(path_to_heart_model), "File does not exist"

    # load model
    model1 = models.HeartModel.load_model(path_to_heart_model)

    assert isinstance(model1, models.BiVentricle), "Expecting model of type BiVentricle"

    # compare contents to original
    assert np.all(model1.left_ventricle.element_ids == model.left_ventricle.element_ids)
    assert np.all(model1.right_ventricle.element_ids == model.right_ventricle.element_ids)

    assert np.all(model1.left_ventricle.endocardium.faces == model.left_ventricle.endocardium.faces)
    assert np.all(
        model1.right_ventricle.endocardium.faces == model.right_ventricle.endocardium.faces
    )

    assert np.all(model1.mesh.tetrahedrons == model.mesh.tetrahedrons)
    assert np.allclose(model1.mesh.nodes, model.mesh.nodes, atol=1e-8)

    pass


@pytest.mark.skipif(skip_test, reason="Requires Windows")
@pytest.mark.parametrize(
    "model_type",
    [models.LeftVentricle, models.BiVentricle, models.FourChamber, models.FullHeart],
)
def test_model_part_names(model_type):
    """Test whether all parts exist in the model."""
    info = _get_test_model_info()
    model = model_type(info)

    if isinstance(model, models.LeftVentricle):
        assert model.part_names == ["Left ventricle"]
    elif isinstance(model, models.BiVentricle):
        assert model.part_names == ["Left ventricle", "Right ventricle", "Septum"]
    elif isinstance(model, models.FourChamber):
        assert model.part_names == [
            "Left ventricle",
            "Right ventricle",
            "Septum",
            "Left atrium",
            "Right atrium",
        ]
    elif isinstance(model, models.FullHeart):
        assert model.part_names == [
            "Left ventricle",
            "Right ventricle",
            "Septum",
            "Left atrium",
            "Right atrium",
            "Aorta",
            "Pulmonary artery",
        ]
