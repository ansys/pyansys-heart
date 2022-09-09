"""Tests high-level heart model class """
import ansys.heart.preprocessor.models as models


def _get_test_model_info() -> models.ModelInfo:
    """Gets a test model info and populates it"""

    # create model info
    info = models.ModelInfo(
        database="Strocchi2020",
        work_directory="my-work-dir",
        path_to_case="path-to-case",
        path_to_simulation_mesh="path-to-simulation-mesh",
        mesh_size=2.0,
    )

    return info


def _get_test_model(model_type: type) -> models.HeartModel:
    """Gets the heart model of inherited class and initializes, and populates"""
    info = _get_test_model()

    model = model_type(info)

    if isinstance(model, models.LeftVentricle):
        print("do something")
    elif isinstance(model, models.BiVentricle):
        print("do something")
    elif isinstance(model, models.FourChamber):
        print("do something")
    elif isinstance(model, models.FullHeart):
        print("do something")

    return model


def test_model_info_dump():
    """Tests dumping of model info to json"""

    pass


def test_model_dump():
    """Tests dumping of model to disk"""
    pass


def test_model_load():
    """Tests loading the model from file"""
    pass
