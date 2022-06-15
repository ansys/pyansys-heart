import pytest
import os
import json
import numpy as np
import copy
import pickle

from ansys.heart.preprocessor.model_information import ModelInformation
from ansys.heart.preprocessor.heart_model import HeartModel
from ansys.heart.preprocessor.vtk_module import read_vtk_polydata_file, read_vtk_unstructuredgrid_file
from ansys.heart.preprocessor.vtk_module import get_tetra_info_from_unstructgrid, get_tri_info_from_polydata
from vtk.numpy_interface import dataset_adapter as dsa 

from conftest import get_assets_folder, get_workdir, clean_directory, create_directory, remove_keys_from_dict, read_file

from common import workflow_extract_mesh

# fixture: run once
@pytest.fixture(autouse=True, scope="module")
def extraction_bi_ventricle():
    global output_dir
    output_dir = workflow_extract_mesh("BiVentricle", get_workdir() )    

    return

# STL Consistency
@pytest.mark.parametrize("stl_file", ["caps_Left_ventricle.stl",
                                    "caps_Right_ventricle.stl",
                                    "closed_volume_left_ventricle.stl",
                                    "closed_volume_right_ventricle.stl"])
def test_stl_consistency(stl_file):
    """Test consistency of caps. Compares the stl files"""
        
    stl_path = os.path.join( output_dir, stl_file )
    stl_path_ref = os.path.join( get_assets_folder(), "reference_models", "BiVentricle", stl_file )

    stl_string = read_file( stl_path )
    stl_string_ref = read_file ( stl_path_ref )

    assert stl_string == stl_string_ref


# VTK Consistency
@pytest.mark.parametrize("vtk_file", ["simulation_mesh.vtk",
                                    "points_endocardium_left_ventricle.vtk",
                                    "points_endocardium_right_ventricle.vtk",                                    
                                    "points_epicardium_left_ventricle.vtk",
                                    "points_epicardium_right_ventricle.vtk",
                                    "points_epicardium-septum_left_ventricle.vtk"])
def test_vtk_consistency(vtk_file):
    vtk_path = os.path.join( output_dir, vtk_file )
    vtk_path_ref = os.path.join(get_assets_folder(), "reference_models", "BiVentricle", vtk_file )

    # read vtk data
    if vtk_file in ["simulation_mesh.vtk"]:
        vtk_object = read_vtk_unstructuredgrid_file( vtk_path )
        vtk_object_ref = read_vtk_unstructuredgrid_file ( vtk_path_ref )
    elif vtk_file.startswith("points"):
        vtk_object = read_vtk_polydata_file(vtk_path )
        vtk_object_ref = read_vtk_polydata_file ( vtk_path_ref )        
    
    vtk_object_dsa = dsa.WrapDataObject( vtk_object )
    vtk_object_dsa_ref = dsa.WrapDataObject( vtk_object_ref )
    
    # compare number of points
    assert vtk_object_dsa.Points.shape[0] == vtk_object_dsa_ref.Points.shape[0], "Number of points not equal"

    # compare cell connectivity
    # NOTE: Compare only vtks that have connectivity table.
    if vtk_file == "simulation_mesh.vtk":
        connectivity = get_tetra_info_from_unstructgrid( vtk_object )[1]
        connectivity_ref = get_tetra_info_from_unstructgrid( vtk_object_ref )[1]
        assert np.all( connectivity == connectivity_ref ), "Connectivity tables not the same"


    # compare point coordinates
    node_to_node_distances = np.linalg.norm( vtk_object_dsa.Points - vtk_object_dsa_ref.Points, axis = 1 )
    epsilon = 1e-3 # max distance (mm) allowed between two nodes

    assert np.max(node_to_node_distances) < epsilon, "Point coordinates not equal"
        