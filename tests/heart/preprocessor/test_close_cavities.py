import pytest
import os
import numpy as np
import shutil

from conftest import (
    get_assets_folder,
    get_workdir,
    read_file,
)

from common import workflow_extract_mesh

# fixture: run once
# @pytest.fixture(autouse=True, scope="module")
# def extraction_left_ventricle():
#     global output_dir
#     output_dir = workflow_extract_mesh("LeftVentricle", get_workdir())

#     yield
#     # cleanup
#     shutil.rmtree(output_dir)


# cavity extraction for LeftVentricle, BiVentricle, and FourChamber models
@pytest.mark.parametrize(
    "model_type",
    [
        "LeftVentricle",
        "BiVentricle",
        "FourChamber",
    ],
)
@pytest.mark.skip(reason="Not finished")
def test_close_cavities(model_type):
    
    # vtk_path_ref = os.path.join(get_assets_folder(), "reference_models", "LeftVentricle", vtk_file)

    # # read vtk data
    # if vtk_file in ["simulation_mesh.vtk"]:
    #     vtk_object = read_vtk_unstructuredgrid_file(vtk_path)
    #     vtk_object_ref = read_vtk_unstructuredgrid_file(vtk_path_ref)
    # elif vtk_file.startswith("points"):
    #     vtk_object = read_vtk_polydata_file(vtk_path)
    #     vtk_object_ref = read_vtk_polydata_file(vtk_path_ref)

    # vtk_object_dsa = dsa.WrapDataObject(vtk_object)
    # vtk_object_dsa_ref = dsa.WrapDataObject(vtk_object_ref)

    # # compare number of points
    # assert (
    #     vtk_object_dsa.Points.shape[0] == vtk_object_dsa_ref.Points.shape[0]
    # ), "Number of points not equal"

    # # compare cell connectivity
    # # NOTE: Compare only vtks that have connectivity table.
    # if vtk_file == "simulation_mesh.vtk":
    #     connectivity = get_tetra_info_from_unstructgrid(vtk_object)[1]
    #     connectivity_ref = get_tetra_info_from_unstructgrid(vtk_object_ref)[1]
    #     assert np.all(connectivity == connectivity_ref), "Connectivity tables not the same"

    # # compare point coordinates
    # node_to_node_distances = np.linalg.norm(
    #     vtk_object_dsa.Points - vtk_object_dsa_ref.Points, axis=1
    # )
    # epsilon = 1e-3  # max distance (mm) allowed between two nodes

    # assert np.max(node_to_node_distances) < epsilon, "Point coordinates not equal"
