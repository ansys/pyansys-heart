"""unit test for HeartModel class"""
import os

from ansys.heart.preprocessor.models import HeartModel, ModelInfo
from conftest import get_assets_folder, get_workdir
import meshio
import numpy as np
import pytest

test_model: HeartModel


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
        test_model.l4cv_axis["normal"], np.array([-0.55653251, -0.09943442, -0.82485414])
    )


def test_compute_left_ventricle_aha17():
    test_model.compute_left_ventricle_anatomy_axis()
    test_model.compute_left_ventricle_aha17()
    assert len(np.where(test_model.aha_ids == 1)[0]) == 10343
    assert len(np.where(test_model.aha_ids == 7)[0]) == 8100
    assert len(np.where(test_model.aha_ids == 17)[0]) == 2189

    test_model.compute_left_ventricle_aha17(p_junction=np.array([0.85, 117, 347]))
    assert len(np.where(test_model.aha_ids == 1)[0]) == 11120
    assert len(np.where(test_model.aha_ids == 7)[0]) == 9885
    assert len(np.where(test_model.aha_ids == 17)[0]) == 2189

    meshio.write_points_cells(
        test_model.mesh.nodes,
        [("tetra", test_model.mesh.tetrahedrons)],
        cell_data={"aha17": [test_model.aha_ids]},
    )


def test_compute_left_ventricle_element_cs():
    test_model.compute_left_ventricle_anatomy_axis()
    test_model.compute_left_ventricle_aha17()
    e_l, e_r, e_c = test_model.compute_left_ventricle_element_cs()

    # get only the relevant mesh: lower part of left ventricle
    ele_ids = np.where(~np.isnan(test_model.aha_ids))[0]
    elems = test_model.mesh.tetrahedrons[ele_ids]
    nodes = test_model.mesh.nodes[np.unique(elems.ravel())]
    _, a = np.unique(elems, return_inverse=True)
    connect = a.reshape(elems.shape)

    # write in to vtk
    meshio.write_points_cells(
        os.path.join(workdir, "lv_lrc_coordinate.vtk"),
        nodes,
        [("tetra", connect)],
        cell_data={"e_l": [e_l], "e_r": [e_r], "e_c": [e_c]},
    )

    assert np.allclose(e_l[0], np.array([0.61060497, -0.72221245, -0.32491653]))
    assert np.allclose(e_r[0], np.array([0.07067868, -0.35894692, 0.93067805]))
    assert np.allclose(e_c[0], np.array([-0.78877506, -0.59124132, -0.16812975]))
