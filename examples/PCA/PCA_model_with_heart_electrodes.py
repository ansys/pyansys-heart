import pandas as pd
from scipy.io import loadmat, savemat
import numpy as np
import pyvista
import pyvista as pv
import numpy as np
import vtk
from scipy.spatial.transform import Rotation
from PCA_model import PCA_model, registration_to_strocchi, save_points_to_vtk

# input files in shared folder
# caeser data (humanshape model) from Internet
# evalues are the eigenvalues for each point
path_to_evalues = r'D:\Data\SSM_Input\input\evalues.mat'


# file use for extract torso part from meanshape
torso_mat_file_path = r'D:\Data\SSM_Input\input\New_torso.mat'
torso_evectors_mat_file_path = r'D:\Data\SSM_Input\input\Torso_evectors.mat'

PCA_predefined_anatomical_points_path = r"D:\Data\SSM_Input\input\PCA_define_ap.csv"

pca_torso_model = PCA_model()
pca_torso_model.initialize_PCA_model(torso_mat_file_path, path_to_evalues, torso_evectors_mat_file_path)



# define pca anatomical points (8)
pca_ap = pd.read_csv(PCA_predefined_anatomical_points_path)
# ['Points_0', 'Points_1', 'Points_2'] correspond to x y z coordinates of aps
ap_points_array = pca_ap[['Points_0', 'Points_1', 'Points_2']].to_numpy()
# use these points as attribute of thje pca model
pca_torso_model.anatomical_points = ap_points_array


# show pca model
body_points = np.array(list(pca_torso_model.body.values()))
body_polydata = pyvista.PolyData(body_points)

anatomical_point_polydata = pyvista.PolyData(np.array(pca_torso_model.anatomical_points))


plotter = pyvista.Plotter()

# plotter.add_mesh(anatomical_point_polydata, color='green', point_size=10)
plotter.add_mesh(body_polydata, color='white', label='Torso', point_size=3)

electrode_colors = {
    "V1": 'red',
    "V2": 'coral',
    "V3": 'orange',
    "V4": 'gold',
    "V5": 'yellow',
    "V6": 'lightgreen',
    "RA": 'green',
    "LA": 'turquoise',
    "LL": 'blue',
    "RL": 'purple',
}

for name, index in pca_torso_model.electrodes.items():
    electrode_point = np.array(pca_torso_model.body[index])
    electrode_polydata = pv.PolyData(electrode_point)
    color = electrode_colors.get(name, 'grey')
    plotter.add_mesh(electrode_polydata, color=color, point_size=10, label=name)

# plotter.add_mesh(model.mesh, color="blue", opacity=0.3)
plotter.add_mesh(pca_torso_model.anatomical_points, color="black", point_size=10)

# plotter.add_legend()
    
plotter.add_legend()
plotter.camera_position = 'xz'
plotter.background_color = 'grey'

plotter.show()


# PCA deformation
'''
Sigma range:
Mode ID = 1 \ Sigma: -+0.00075
Mode ID = 2 \ Sigma: -+0.0008
Mode ID = 3 \ Sigma: -+0.001
'''
mode_id=1
sigma=0.00075
deformed_pca_save_file = f'D:\\Data\\SSM_Input\\output\\mode_{mode_id}_sigma_{sigma}.vtk'

pca_torso_model.deformation(mode_id=mode_id, sigma=sigma, save=False, filename=deformed_pca_save_file, path_to_evectors=torso_evectors_mat_file_path)



# extract electrodes position from defromed torso
electrode_points = []
for name, index in pca_torso_model.electrodes.items():
    electrode_point = np.array(pca_torso_model.body[index])
    electrode_points.append(electrode_point)

electrode_points_array = np.array(electrode_points)



# registration to patient's heart
import os
from pathlib import Path
from ansys.heart.preprocessor.mesh.objects import Point
import ansys.heart.preprocessor.models.v0_1.models as models
from ansys.heart.simulator.simulator import DynaSettings, EPSimulator
import numpy as np
import pandas as pd

__file__ = r'D:\REPOS\pyansys-heart\examples\simulator\doc_EP_simulator_fourchamber.py'

workdir = Path(
    Path(__file__).resolve().parents[2], "downloads", "Strocchi2020", "01", "FourChamber"
)

path_to_model = os.path.join(workdir, "heart_model.pickle")

if not os.path.isfile(path_to_model):
    raise FileExistsError(f"{path_to_model} not found")

model: models.FourChamber = models.HeartModel.load_model(path_to_model)


# download a patient model: case 01
fix_points = [cap.centroid for cap in model.left_ventricle.caps]
fix_points += [cap.centroid for cap in model.right_ventricle.caps]
# Extract apex coordinates from the left and right ventricles
fix_points += [apex.xyz for apex in model.left_ventricle.apex_points]
fix_points += [apex.xyz for apex in model.right_ventricle.apex_points]
# Convert the list of points to a NumPy array
fix_points = np.array(fix_points)

moving_points = pca_torso_model.anatomical_points

transformed_model = registration_to_strocchi(fix_points, moving_points, electrode_points_array)
# save_points_to_vtk(transformed_model,'new_electrodes1.vtk')

print(transformed_model)


# Reading the CSV content
PCA_deformation_scenario_df = pd.read_csv(r"D:\Data\SSM_Input\input\pca_scenarii.csv")

transformed_model_save = []
# Printing each row
for index, row in PCA_deformation_scenario_df.iterrows():

    pca_torso_model.initialize_PCA_model(torso_mat_file_path, path_to_evalues, torso_evectors_mat_file_path)
    mode_id = int(row['mode_id'])
    sigma = float(row['sigma'])

    # deformation
    deformed_pca_save_file = f'D:\Data\SSM_Input\output\electrodes_mode_{mode_id}_sigma_{sigma}.vtk'
    pca_torso_model.deformation(mode_id=mode_id, sigma=sigma, save=False, filename=deformed_pca_save_file, path_to_evectors=torso_evectors_mat_file_path)

    electrode_points = []
    for name, index in pca_torso_model.electrodes.items():
        electrode_point = np.array(pca_torso_model.body[index])
        electrode_points.append(electrode_point)

    electrode_points_array = np.array(electrode_points)

    # registration
    transformed_model = registration_to_strocchi(fix_points, moving_points, electrode_points_array)
    transformed_model_str = np.array2string(transformed_model, separator=',')
    transformed_model_save.append({'mode_id': mode_id, 'sigma': sigma, 'transformed_electrodes': transformed_model_str})

    print(transformed_model)
    # save transformed electrodes
    save_points_to_vtk(transformed_model,deformed_pca_save_file)
    
df_transformed_model_save = pd.DataFrame(transformed_model_save)

csv_save_path = r'D:\Data\SSM_Input\output\pca_electrodes_results3.csv'
df_transformed_model_save.to_csv(csv_save_path, index=False)


import numpy as np

electrodes_positions = [[-79.63069638, 8.3915733, 408.82855283], [15.61357471, 7.79292347, 423.30215084], [71.22635286, -0.81217516, 362.57101383], [111.43449808, 12.25370598, 332.27745863], [152.46161301, 48.62191579, 339.36410283], [198.75978967, 126.29888517, 389.04616454], [-259.08174793, 109.44650165, 533.73951965], [209.88567613, 83.15165099, 563.61574278], [-121.92311261, 19.60510742, 64.27745959], [185.18313656, 54.02360627, 92.02302318]]

electrodes_array = np.array(electrodes_positions)

shape = electrodes_array.shape

print("Shape of electrodes array:", shape)




# Reading the CSV content
PCA_deformation_scenario_df = pd.read_csv(r"D:\Data\SSM_Input\input\pca_scenarii.csv")

transformed_model_save = []
# Printing each row
for index, row in PCA_deformation_scenario_df.iterrows():

    pca_torso_model.initialize_PCA_model(torso_mat_file_path, path_to_evalues, torso_evectors_mat_file_path)
    mode_id = int(row['mode_id'])
    sigma = float(row['sigma'])

    # deformation
    deformed_pca_save_file = f'D:\Data\SSM_Input\output\electrodes_mode_{mode_id}_sigma_{sigma}.vtk'
    pca_torso_model.deformation(mode_id=mode_id, sigma=sigma, save=False, filename=deformed_pca_save_file, path_to_evectors=torso_evectors_mat_file_path)

    electrode_points = []
    for name, index in pca_torso_model.electrodes.items():
        electrode_point = np.array(pca_torso_model.body[index])
        electrode_points.append(electrode_point)

    electrode_points_array = np.array(electrode_points)

    # registration
    transformed_model = registration_to_strocchi(fix_points, moving_points, electrode_points_array)
    transformed_model_str = np.array2string(transformed_model, separator=',')
    transformed_model_save.append({'mode_id': mode_id, 'sigma': sigma, 'transformed_electrodes': transformed_model_str})

    print(transformed_model)
    # save transformed electrodes
    save_points_to_vtk(transformed_model,deformed_pca_save_file)
    
df_transformed_model_save = pd.DataFrame(transformed_model_save)

csv_save_path = r'D:\Data\SSM_Input\output\pca_electrodes_results3.csv'
df_transformed_model_save.to_csv(csv_save_path, index=False)