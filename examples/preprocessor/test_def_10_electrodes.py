import os
from pathlib import Path

from ansys.heart.misc.downloader import download_case, unpack_case
import ansys.heart.preprocessor.models as models

import pyvista
import numpy as np
import vtk

from vtkmodules.vtkCommonDataModel import vtkIterativeClosestPointTransform

os.environ["USE_OLD_HEART_MODELS"] = "1"

# __file__ = r"c:\Users\xuhu\pyheart-lib\examples\preprocessor\doc_ECG_coordinates,py"

case_file = str(
    Path(Path(__file__).resolve().parents[2], "downloads", "Strocchi2020", "01", "01.case")
)
download_folder = str(Path(Path(__file__).resolve().parents[2], "downloads"))
workdir = str(
    Path(Path(__file__).resolve().parents[2], "downloads", "Strocchi2020", "01", "Biv")
)
path_to_model = str(Path(workdir, "heart_model.pickle"))


if not os.path.isfile(case_file):
    path_to_downloaded_file = download_case(
        "Strocchi2020", 1, download_folder=download_folder, overwrite=False
    )
    unpack_case(path_to_downloaded_file)

model: models.BiVentricle = models.BiVentricle.load_model(path_to_model)


move_points = np.array([
    [94.35242091,75.99022315,213.31654731], # mitral-valve
    [81.90321388,57.90000882,205.76663367], # aortic-valve
    [121.58912558,89.76497459,223.29557159], # tricuspid-valve
    [67.14045655,102.49380179,216.20654707], # pulmonary-valve
    [76.04229182019685,66.53094359081156,297.7182142431582], # left endo
    [75.08606835375224,66.33759424571653,302.2811669120656], # left epi
    [70.87069056682236,84.83837198547876,295.6765864478138], # right endo
    [70.54655746919204,84.50457846174797,297.2737993295601], # right epi   
])

electrode_positions = np.array([
    [x, y, z] for x, y, z in [
        [91.69106809237354, 167.4055272828183, 251.0416838617331],  # V1
        [114.07772933063883, 123.13918227704727, 291.5674747053013],  # V2
        [97.01364431022192, 109.7927312489943, 317.44575378942824],  # V3
        [81.88880486815538, 71.3859176743268, 349.4243690358569],  # V4
        [98.99550734014053, 15.879947224122954, 348.26233938958114],  # V5
        [106.23537044908527, -44.085603837273695, 329.16199248487465],  # V6
        [170.75567690191764, 234.51300755277248, 77.85629801985534],  # RA
        [262.9539413249434, -2.9189733795951724, 261.5745131716608],  # LA
        [-134.9640236606803, 197.05542132895272, 257.6409644703581],  # RL
        [-70.4506948627224, 22.20437217827933, 400.2792599184596]  # LL
    ]
])

# # Create a vtkPoints object and add the points to it
# vtk_points = vtk.vtkPoints()
# for point in np.concatenate((move_points, electrode_positions)):
#     vtk_points.InsertNextPoint(point.tolist())

# # Create a vtkPolyData object
# poly_data = vtk.vtkPolyData()
# poly_data.SetPoints(vtk_points)

# # Write the poly_data to a VTK file
# writer = vtk.vtkPolyDataWriter()
# writer.SetFileName(r"C:\Users\xuhu\OneDrive - ANSYS, Inc\Desktop\Temp\Pyansys-heart\test_p.vtk")
# writer.SetInputData(poly_data)
# writer.Write()

# plotter = pyvista.Plotter()

# plotter.add_mesh(move_points,color="red")
# # p2.add_mesh(electrode_positions, color="blue", opacity=0.3)
# # p2.add_mesh(model.mesh, color="red", opacity=0.3)
# plotter.add_mesh(electrode_positions, color="blue")

# # Set the background color and show the plotter
# plotter.background_color = "white"
# plotter.show()

transformed_electrodes = model.define_ECG_coordinates(move_points, electrode_positions)

print(model.electrodes)

plotter = pyvista.Plotter()

plotter.add_mesh(model.mesh,color="red", opacity=0.3)
# p2.add_mesh(electrode_positions, color="blue", opacity=0.3)
# p2.add_mesh(model.mesh, color="red", opacity=0.3)
plotter.add_mesh(transformed_electrodes, color="blue", opacity=1)

# Set the background color and show the plotter
plotter.background_color = "white"
plotter.show()


# Create a vtkPoints object and add the points to it
# vtk_points = vtk.vtkPoints()
# for point in np.concatenate((transformed_electrodes)):
#     vtk_points.InsertNextPoint(point.tolist())

# # Create a vtkPolyData object
# poly_data = vtk.vtkPolyData()
# poly_data.SetPoints(transformed_electrodes)

# # Write the poly_data to a VTK file
# writer = vtk.vtkPolyDataWriter()
# writer.SetFileName(r"C:\Users\xuhu\OneDrive - ANSYS, Inc\Desktop\Temp\Pyansys-heart\test_p.vtk")
# writer.SetInputData(poly_data)
# writer.Write()


'''
# ECG positions for Strocchi heart:
electrode_positions = np.array([
    [-29.893285751342773, 27.112899780273438, 373.30865478515625],   # V1
    [33.68170928955078, 30.09606170654297, 380.5427551269531],       # V2
    [56.33562469482422, 29.499839782714844, 355.533935546875],       # V3
    [100.25729370117188, 43.61333465576172, 331.07635498046875],     # V4
    [140.29800415039062, 81.36004638671875, 349.69970703125],        # V5
    [167.9899139404297, 135.89862060546875, 366.18634033203125],     # V6
    [-176.06332397460938, 57.632076263427734, 509.14202880859375],   # RA
    [133.84518432617188, 101.44053649902344, 534.9176635742188],     # LA
    [-103.60343933105469, 64.02100372314453, 160.0018310546875],     # RL
    [128.9441375732422, 92.85327911376953, 173.07363891601562]       # LL
])
'''