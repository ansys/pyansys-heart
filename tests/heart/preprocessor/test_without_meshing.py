import pytest
import os
import numpy as np
import shutil

from conftest import (
    get_assets_folder,
)

from ansys.heart.preprocessor.model_information import ModelInformation
from ansys.heart.preprocessor.heart_model import HeartModel

# cavity extraction for LeftVentricle, BiVentricle, and FourChamber models
@pytest.mark.parametrize(
    "model_type",
    [
        "LeftVentricle",
        # "BiVentricle",
        # "FourChamber",
    ],
)
@pytest.mark.skip(reason="Placeholder test")
def test_map_data(model_type):

    # get asset paths
    reference_directory = os.path.join(
        get_assets_folder(), "reference_models", model_type )

    model_info_path = os.path.join(reference_directory, "model_info.json")
    # specify needed paths
    path_to_raw_mesh = os.path.join( get_assets_folder(), "cases", "strocchi2020", "01", "01.case" )
    path_output_meshing = os.path.join( os.path.dirname(model_info_path), "output_meshing.vtk")

    # Initialize model information
    model_info = ModelInformation()
    model_info.load_from_file( model_info_path )
    model_info.path_original_mesh = path_to_raw_mesh

    model = HeartModel(model_info)

    # load raw mesh data
    model._mesh.load_raw_mesh()

    # skip remeshing - but set volume mesh directly from file
    model._mesh.set_volume_mesh_vtk(path_output_meshing)
    model._mesh.map_data_to_remeshed_volume()

    # TODO: do assert: What conditions determine test pass/fail?

    return
