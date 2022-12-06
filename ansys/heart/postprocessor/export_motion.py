"""Export motion."""
import os
import subprocess

from ansys.heart.preprocessor.mesh.vtkmethods import (
    read_vtk_polydata_file,
    vtk_cutter,
    write_vtkdata_to_vtkfile,
)
from ansys.heart.preprocessor.models import HeartModel

# import pkg_resources

model: HeartModel
model = HeartModel.load_model(r"D:\Heart20\last_model\heart_model.pickle")
nb_cut = 6

pvpython = r"C:\Program Files\ParaView 5.9.0-Windows-Python3.8-msvc2017-64bit\bin\pvpython.exe"
# pkg_resources.resource_filename("ansys.heart.misc.paraview_marco", "d3plot_to_vtk.pvpy")
pv_script = r"D:\pyheart-lib\ansys\heart\misc\paraview_marco\d3plot_to_vtk.pvpy"
d3plot_file = "D:/Heart20/last_model/main/d3plot"
out_file = "D:/Heart20/last_model/lv_surface/lv.vtk"
subprocess.call([pvpython, pv_script, d3plot_file, out_file, "1"])

model.compute_left_ventricle_anatomy_axis()
p_start = model.short_axis["center"]
p_end = model.left_ventricle.apex_points[1].xyz

w_dir = f"D:\Heart20\last_model\l4cv_axe"
os.makedirs(w_dir, exist_ok=True)

for iframe in range(41):
    f = f"D:\Heart20\last_model\lv_surface\lv_{iframe}.vtk"
    surface = read_vtk_polydata_file(f)
    res = vtk_cutter(
        surface, {"center": model.l4cv_axis["center"], "normal": model.l4cv_axis["normal"]}
    )
    write_vtkdata_to_vtkfile(res, os.path.join(w_dir, f"contour_{iframe}.vtk"))

w_dir = f"D:\Heart20\last_model\l2cv_axe"
os.makedirs(w_dir, exist_ok=True)

for iframe in range(41):
    f = f"D:\Heart20\last_model\lv_surface\lv_{iframe}.vtk"
    surface = read_vtk_polydata_file(f)
    res = vtk_cutter(
        surface, {"center": model.l2cv_axis["center"], "normal": model.l2cv_axis["normal"]}
    )
    write_vtkdata_to_vtkfile(res, os.path.join(w_dir, f"contour_{iframe}.vtk"))

for icut in range(nb_cut):
    w_dir = f"D:\Heart20\last_model\short_axe{icut}"
    os.makedirs(w_dir, exist_ok=True)
    p_cut = p_start + (p_end - p_start) * icut / nb_cut

    for iframe in range(41):
        f = f"D:\Heart20\last_model\lv_surface\lv_{iframe}.vtk"
        surface = read_vtk_polydata_file(f)
        res = vtk_cutter(surface, {"center": p_cut, "normal": model.short_axis["normal"]})
        write_vtkdata_to_vtkfile(res, os.path.join(w_dir, f"contour_{iframe}.vtk"))
