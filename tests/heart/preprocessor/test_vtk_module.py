import pytest
import os
import vtk
import numpy as np

from conftest import ROOT_FOLDER
from conftest import get_assets_folder, get_workdir
from ansys.heart.preprocessor.vtk_module import compute_volume_stl, vtk_surface_to_stl

@pytest.mark.skip(reason="Example test")
def test_compute_volume():    
    """Tests compute volume function from stl. Use unit-sphere as reference"""
    radius = 1
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius( radius )
    sphere.GenerateNormalsOn()
    sphere.SetPhiResolution(32)
    sphere.SetThetaResolution(32)
    sphere.Update()
    workdir = get_workdir()
    file_path = os.path.join( workdir,  "sphere.stl" )
    # write to stl
    vtk_surface_to_stl( sphere.GetOutput(), file_path )
    # compute volume
    volume_stl = compute_volume_stl ( file_path )
    os.remove( file_path )

    # compute reference volume
    volume_ref = (4/3) * np.pi * radius ** 3
    epsilon = 0.01
    epsilon_abs =  epsilon * volume_ref  # within 1% of reference volume
    assert abs( volume_stl - volume_ref ) < epsilon_abs