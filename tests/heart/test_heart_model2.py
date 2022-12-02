"""unit test for HeartModel class"""
import os

from ansys.heart.preprocessor.models import HeartModel, ModelInfo
from conftest import get_assets_folder, get_workdir
import numpy as np
import pytest


@pytest.fixture(autouse=True, scope="module")
def get_test_model():
    """Load a model for tests."""

    global test_model

    # init the model with dummy information
    test_model = HeartModel(
        ModelInfo(
            database="Strocchi2020",
            work_directory=get_workdir(),
            path_to_case="path-to-case",
            path_to_simulation_mesh="path-to-simulation-mesh",
            mesh_size=2.0,
        )
    )

    # load model information from asset folder
    assets_folder = get_assets_folder()
    path_to_reference_model = os.path.join(
        assets_folder,
        "reference_models",
        "strocchi2020",
        "01",
        "BiVentricleRefactored",
        "heart_model.pickle",
    )
    test_model = test_model.load_model(path_to_reference_model)

    # if necessary, write tests results in this directory
    global workdir
    workdir = os.path.join(".", "test_model")
    os.makedirs(workdir, exist_ok=True)
    test_model.mesh.write_to_vtk(os.path.join(workdir, "model.vtk"))
    # yield
    # if os.path.isdir(workdir):
    #     shutil.rmtree(workdir)


def test_compute_left_ventricle_anatomy_axis():
    test_model.compute_left_ventricle_anatomy_axis()
    print()
    print(test_model.l4cv_axis)
    print(test_model.short_axis)
    print(test_model.l2cv_axis)

    # can be visualized by heart/misc/paraview_marco/plot_anatomy_axis.pvpy

    assert np.allclose(
        test_model.l4cv_axis["normal"], np.array([0.55653251, 0.09943442, 0.82485414])
    )
