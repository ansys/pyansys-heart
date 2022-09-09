import os
import shutil

from ansys.heart.preprocessor._deprecated_heart_model import HeartModel
from ansys.heart.preprocessor._deprecated_model_information import ModelInformation
from ansys.heart.preprocessor.mesh.vtkmethods import vtk_read_mesh_file
from conftest import create_directory, get_assets_folder, get_workdir
import numpy as np
import pytest
from vtk.numpy_interface import dataset_adapter as dsa

# import shutil



# cavity extraction for LeftVentricle, BiVentricle, and FourChamber models
@pytest.mark.parametrize(
    "model_type",
    [
        "LeftVentricle",
        "BiVentricle",
        "FourChamber",
    ],
)
def test_map_data(model_type):

    work_dir = os.path.join(get_workdir(), model_type)

    create_directory(work_dir)

    # get asset paths
    reference_directory = os.path.join(get_assets_folder(), "reference_models", model_type)

    model_info_path = os.path.join(reference_directory, "model_info.json")
    # specify needed paths
    path_to_raw_mesh = os.path.join(get_assets_folder(), "cases", "strocchi2020", "01", "01.case")
    path_output_meshing = os.path.join(os.path.dirname(model_info_path), "output_meshing.vtk")

    # Initialize model information
    model_info = ModelInformation()
    model_info.load_from_file(model_info_path)
    model_info.path_original_mesh = path_to_raw_mesh
    model_info.working_directory = work_dir

    model = HeartModel(model_info)

    # perform some necessary steps
    model._mesh.load_raw_mesh()
    model._mesh.extract_parts()
    model._mesh.add_cavities()

    # skip remeshing - but set volume mesh directly from file
    model._mesh.set_volume_mesh_vtk(path_output_meshing)
    # map data from original dataset to this volume mesh
    model._mesh.map_data_to_remeshed_volume()

    # simulation mesh is reference - load to compare
    path_to_reference_model = os.path.join(reference_directory, "simulation_mesh.vtk")

    # TODO: compare reference model with model after mapping
    reference_model_dsa = dsa.WrapDataObject(vtk_read_mesh_file(path_to_reference_model))
    model_dsa = dsa.WrapDataObject(model._mesh._vtk_volume)

    assert np.all(model_dsa.Points == reference_model_dsa.Points), "Points are not the same"
    assert np.all(model_dsa.Cells == reference_model_dsa.Cells), "Cells are not the same"
    assert np.all(
        sorted(model_dsa.PointData.keys()) == sorted(reference_model_dsa.PointData.keys())
    ), "PointData keys are not the same"
    assert np.all(
        sorted(model_dsa.CellData.keys()) == sorted(reference_model_dsa.CellData.keys())
    ), "CellData keys are not the same"

    for key in model_dsa.PointData.keys():
        value = model_dsa.PointData[key]
        reference_value = reference_model_dsa.PointData[key]
        assert np.all(value == reference_value), "Value of PointData %s does not match" % key

    for key in model_dsa.CellData.keys():
        value = model_dsa.CellData[key]
        reference_value = reference_model_dsa.CellData[key]
        assert np.all(value == reference_value), "Value of PointData %s does not match" % key

    # cleanup
    shutil.rmtree(work_dir)

    return
