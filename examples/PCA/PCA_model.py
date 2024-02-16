import pandas as pd
import numpy as np
import pyvista
import pyvista as pv
import numpy as np
import vtk

from scipy.io import loadmat, savemat
from scipy.spatial.transform import Rotation
from scipy.optimize import minimize


def is_close_relative(value_a, value_b, tolerance=1e-4):
    """Check if two values are close to each other within a certain relative tolerance"""
    if value_a == value_b == 0:
        return True
    relative_error = abs(value_a - value_b) / max(abs(value_a), abs(value_b))
    return relative_error < tolerance

def filter_data_relative_precision(mean_shape, eigenvectors, csv_points, tolerance=1e-4):
    """Extract torso points and eigenvectors from mean humanshape model based on a given csv"""
    points_from_csv = csv_points[['Points_0', 'Points_1', 'Points_2']].values
    num_points = len(mean_shape) // 3
        
    # identifying points to keep
    indexes_to_keep = []    
    for i in range(num_points):
        point = mean_shape[i], mean_shape[i + num_points], mean_shape[i + 2*num_points]
        if any(all(is_close_relative(point[dim], p[dim], tolerance) for dim in range(3)) for p in points_from_csv):
            indexes_to_keep.append(i)
    
    # filtering mean_shape
    filtered_mean_shape = np.concatenate([mean_shape[indexes_to_keep], 
                                          mean_shape[num_points:2*num_points][indexes_to_keep], 
                                          mean_shape[2*num_points:][indexes_to_keep]])
    
    # filtering eigenvectors
    filtered_eigenvectors = np.concatenate([eigenvectors[:, indexes_to_keep], 
                                            eigenvectors[:, num_points:2*num_points][:, indexes_to_keep], 
                                            eigenvectors[:, 2*num_points:][:, indexes_to_keep]], axis=1)
    
    return filtered_mean_shape, filtered_eigenvectors

def save_points_to_vtk(points, filename):
    """Save the given points to a vtk file"""
    num_points = points.shape[0]

    vtk_points = vtk.vtkPoints()
    for i in range(num_points):
        vtk_points.InsertPoint(i, points[i, 0], points[i, 1], points[i, 2])

    vertices = vtk.vtkCellArray()
    for i in range(num_points):
        vertices.InsertNextCell(1)
        vertices.InsertCellPoint(i)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.SetVerts(vertices)

    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(polydata)
    writer.Write()


def rigid_transformation_scipy(fix_points, moving_points, moving_model):
    """Rigid transformation with scipy"""
    # define the initial transformation parameters
    random_quaternion = Rotation.random().as_quat()
    initial_params = np.zeros(7)
    initial_params[:3] = np.random.rand(3)  # Random translation
    initial_params[3:] = random_quaternion

    # constrain quaternion components
    constraints = ({'type': 'eq', 'fun': lambda params: 1.0 - np.sum(params[3:] ** 2)})

    # rigid transform
    def rigid_transform(params, points):
        translation = params[:3]
        quaternion = params[3:]
        quaternion /= np.linalg.norm(quaternion)
        rotation_matrix = Rotation.from_quat(quaternion).as_matrix()
        transformed_points = np.dot(points - translation, rotation_matrix.T)
        return transformed_points

    # calculate the optimized parameter
    def objective_function(params, fixed_points, moving_points):
        transformed_points = rigid_transform(params, moving_points)
        distance = np.sum(np.square(transformed_points - fixed_points))
        return distance
    
    # optimized parameters
    result = minimize(
        objective_function, 
        initial_params, 
        args=(fix_points, moving_points), 
        method='L-BFGS-B', 
        constraints=constraints)
    optimal_params = result.x

    # move model based on optimal parameters
    transformed_model = moving_model.copy()
    transformed_model.points = rigid_transform(optimal_params, moving_model.points)

    return transformed_model


def registration_to_strocchi(fix_points, moving_points, electrode_points_array):
    '''registrer deformed electrodes to strocchi heart'''
    from scipy.optimize import minimize

    # Define the initial transformation parameters
    random_quaternion = Rotation.random().as_quat()
    initial_params = np.zeros(7)
    initial_params[:3] = np.random.rand(3)  # Random translation
    initial_params[3:] = random_quaternion

    # Constrain quaternion components
    constraints = ({'type': 'eq', 'fun': lambda params: 1.0 - np.sum(params[3:] ** 2)})

    # rigid transform function
    def rigid_transform(params, points):
        translation = params[:3]
        quaternion = params[3:]
        quaternion /= np.linalg.norm(quaternion)
        rotation_matrix = Rotation.from_quat(quaternion).as_matrix()
        transformed_points = np.dot(points - translation, rotation_matrix.T)
        return transformed_points

    # function to calculate the optimized parameter
    def objective_function(params, fixed_points, moving_points):
        transformed_points = rigid_transform(params, moving_points)
        distance = np.sum(np.square(transformed_points - fixed_points))
        return distance

    # Get the optimized parameters
    result = minimize(
        objective_function, 
        initial_params, 
        args=(fix_points, moving_points), 
        method='L-BFGS-B', 
        constraints=constraints)
    optimal_params = result.x

    transformed_model = electrode_points_array.copy()
    transformed_model = rigid_transform(optimal_params, transformed_model)

    return transformed_model

class PCA_model:
    def __init__(self):
        self.body = {}
        self.evectors = {}
        self.evalues = {}

        self.heart_case_id_of_reference: str = '01'
        # self.registration_matrix_file_path = r'D:\xuhu\humanshape\PCA_model_registration\combined_transform_matrix.csv'
        self.electrode_file_path = r"D:\Data\SSM_Input\input\PCA_electrodes.csv"
        
        # only if input data is meanshape: may be used later (_extract_torso), not needed at the moment
        # self.torso = None
        # self.torso_evectors = None
        
        self.electrodes = {}
        self.anatomical_points = {}

    def initialize_PCA_model(self, meanshape_file, eigenvalue_file, eignevectors_file):
        # READ pca DATA
        self._read_pca_data(meanshape_file, eigenvalue_file, eignevectors_file)
        self._define_electrodes()        
        # original_anatomical_points = self._extract_anatomical_points_of_reference()
        # self._caculate_heart_anatomical_points(
        #     self.registration_matrix_file_path, 
        #     original_anatomical_points
        # )
        return

    def _read_pca_data(self, path_to_vtk_file, path_to_evalues, path_to_evectors):
        # read/import 3 .mat file 
        # 1st step: you will use the 3 ones created manually and focusing on the torso (store it in PCA data folder)
        # Later we will manage all/specific one PCA from the website (matlab)
        # SAVE all the data in dictionary: rearranging data for evectors(x1 x2 x3....y1 y2 y3... z1 z2 z3...)
        body_data = loadmat(path_to_vtk_file)['points']
        evectors_data = loadmat(path_to_evectors)['evectors']
        evalues_data = loadmat(path_to_evalues)['evalues']

        self.body = {i: body_data[i, :] for i in range(body_data.shape[0])}

        # self.evectors = {i: evectors_data[:, i] for i in range(evectors_data.shape[1])}
        num_points = evectors_data.shape[0]
        num_cols = evectors_data.shape[1]

        # Rearranging data
        for i in range(num_points):
            x_values = evectors_data[i, :num_cols // 3]
            y_values = evectors_data[i, num_cols // 3: 2 * num_cols // 3]
            z_values = evectors_data[i, 2 * num_cols // 3:]
            self.evectors[i] = np.column_stack((x_values, y_values, z_values))

        self.evalues = {i: evalues_data[0, i] for i in range(evalues_data.size)}

    def _extract_torso(self, torso_file_path):
        # FOR LATER
        # This method should be implemented to extract torso data from the body data
        # For example, it could involve selecting a subset of points or applying a filter
        torso_model = self.body.copy()
        test_torso_model = np.hstack(torso_model.T.flatten())
        csv_points_df = pd.read_csv(torso_file_path)
        self.torso, self.torso_evectors = filter_data_relative_precision(test_torso_model, self.evectors, csv_points_df)

    def _define_electrodes(self):
        csv_data = pd.read_csv(self.electrode_file_path)
        electrode_positions = csv_data[['Points_0', 'Points_1', 'Points_2']].values
        
        # find the closest point -> define electrodes' references
        for i, position in enumerate(electrode_positions):
            min_distance = float('inf')
            closest_index = None

            for index, point in self.body.items():
                distance = np.sum((np.array(point) - position) ** 2)
                if distance < min_distance:
                    min_distance = distance
                    closest_index = index

            electrode_name = "V" + str(i + 1) if i < 6 else ["RA", "LA", "RL", "LL"][i - 6]
            self.electrodes[electrode_name] = closest_index

    def _extract_anatomical_points_of_reference(self):
        # import os
        # from pathlib import Path

        # import ansys.heart.preprocessor.models as models
        # import ansys.heart.writer.dynawriter as writers
        # from ansys.heart.misc.downloader import download_case, unpack_case
        # from ansys.heart.simulator.simulator import run_lsdyna
        # import pyvista
        # import numpy as np
        # import vtk

        # from vtkmodules.vtkCommonDataModel import vtkIterativeClosestPointTransform
        # from ansys.heart.simulator.simulator import DynaSettings, EPSimulator

        # os.environ["USE_OLD_HEART_MODELS"] = "1"

        # __file__ = r"c:\Users\xuhu\pyheart-lib\examples\preprocessor\doc_ECG_coordinates.py"

        # case_file = str(
        #     Path(Path(__file__).resolve().parents[2], "downloads", "Strocchi2020", self.heart_case_id_of_reference, "01.case")
        # )
        # download_folder = str(Path(Path(__file__).resolve().parents[2], "downloads"))
        # workdir = str(
        #     Path(Path(__file__).resolve().parents[2], "downloads", "Strocchi2020", self.heart_case_id_of_reference, "Biv")
        # )
        # path_to_model = str(Path(workdir, "heart_model.pickle"))


        # if not os.path.isfile(case_file):
        #     path_to_downloaded_file = download_case(
        #         "Strocchi2020", 1, download_folder=download_folder, overwrite=False
        #     )
        #     unpack_case(path_to_downloaded_file)

        # model: models.BiVentricle = models.BiVentricle.load_model(path_to_model)
        from ansys.heart.preprocessor.mesh.objects import Point

        template_anatomical_points = [
            [ 85.65327324,  86.50614136, 213.69885421],
            [ 77.74852222, 109.93585137, 217.75954732],
            [ 66.40093437,  52.36786898, 204.40208587],
            [127.88220618,  74.93271448, 229.93636546],
            [ 62.36470533,  71.99060696, 293.69408799],
            [ 59.75769555,  69.96465816, 298.15399064],
            [ 88.60548882,  80.38739545, 295.37047427],
            [ 89.12088316,  82.27674491, 298.51971067]
        ]

        return template_anatomical_points
        
    def _caculate_heart_anatomical_points(self, registration_matrix_file_path, anatomical_points_template):
        registration_matrix = np.loadtxt(registration_matrix_file_path, delimiter=',')

        ones = np.ones((anatomical_points_template.shape[0], 1))
        anatomical_points_homogeneous = np.hstack((anatomical_points_template, ones))
        # Apply transfer matrix to anatomical points
        transformed_anatomical_points = registration_matrix @ anatomical_points_homogeneous.T

        # Conversion from homogeneous coordinates back to 3D coordinates
        transformed_anatomical_points = transformed_anatomical_points.T[:, :3]
        # transformed_anatomical_points = pyvista.PolyData(transformed_anatomical_points
        self.anatomical_points = transformed_anatomical_points


    def deformation(self, mode_id, sigma, save: bool, filename: str, path_to_evectors: str):
        # return polydata, and if save == True, save the polydata as vtk using filename
        
        # reshape PCA model
        test_mean_model = np.array(list(self.body.values()))

        test_mean_model2 = np.ndarray.tolist(test_mean_model[:,0])
        test_mean_model2.extend(np.ndarray.tolist(test_mean_model[:,1]))
        test_mean_model2.extend(np.ndarray.tolist(test_mean_model[:,2]))
        test_mean_model2 = np.array(test_mean_model2)

        # caculate coef
        # coef = np.zeros(len(self.evalues.values()))
        eigenvalues = np.array(list(self.evalues.values()))
        coef = np.zeros(eigenvalues.shape)
        coef[mode_id-1] = eigenvalues[mode_id-1]*3*sigma

        # caculate evectors
        # path_to_evectors = r'D:\xuhu\humanshape\PCA_model_registration\input\Torso_evectors.mat'
        eigenvectors = loadmat(path_to_evectors)['evectors']
        e_evectors = np.dot(np.diag(coef),eigenvectors)

        deformed_model = np.transpose(test_mean_model2) + np.sum(e_evectors, axis=0)


        # reshape deformed PCA model
        deformed_x = deformed_model[0:int(deformed_model.size/3)-1]
        deformed_y = deformed_model[int(deformed_model.size/3): int(2*deformed_model.size/3)-1]
        deformed_z = deformed_model[int(2*deformed_model.size/3):-1]
        deformed_model_new = np.transpose([deformed_x, deformed_y, deformed_z])

        # PCA deformation
        centroid_mean = np.mean(test_mean_model, axis=0)
        centroid_deformed = np.mean(deformed_model_new, axis=0)

        displacement = centroid_mean - centroid_deformed
        moved_deformed_model_new = deformed_model_new + displacement

        if save == True:
            save_points_to_vtk(moved_deformed_model_new, filename)

        # update PCA_model
        for i in range(len(moved_deformed_model_new)):
            self.body[i] = moved_deformed_model_new[i]
    
        return


    def _rigid_transformation_scipy(self, fix_points, moving_points, moving_model):
        from scipy.optimize import minimize

        # Define the initial transformation parameters
        random_quaternion = Rotation.random().as_quat()
        initial_params = np.zeros(7)
        initial_params[:3] = np.random.rand(3)  # Random translation
        initial_params[3:] = random_quaternion

        # Constrain quaternion components
        constraints = ({'type': 'eq', 'fun': lambda params: 1.0 - np.sum(params[3:] ** 2)})

        # rigid transform function
        def rigid_transform(params, points):
            translation = params[:3]
            quaternion = params[3:]
            quaternion /= np.linalg.norm(quaternion)
            rotation_matrix = Rotation.from_quat(quaternion).as_matrix()
            transformed_points = np.dot(points - translation, rotation_matrix.T)
            return transformed_points

        # function to calculate the optimized parameter
        def objective_function(params, fixed_points, moving_points):
            transformed_points = rigid_transform(params, moving_points)
            distance = np.sum(np.square(transformed_points - fixed_points))
            return distance
        
        # Get the optimized parameters
        result = minimize(
            objective_function, 
            initial_params, 
            args=(fix_points, moving_points), 
            method='L-BFGS-B', 
            constraints=constraints)
        optimal_params = result.x

        transformed_model = moving_model.copy()
        transformed_model.points = rigid_transform(optimal_params, moving_model.points)

        return transformed_model

    def _rigid_transformation_landmark(self, fix_points, moving_points, moving_model):
        # 
        return

    def register_to_pca(self, model, registration_type: int, save: bool, filename: str):
        # register patient
        # return registered_mesh
        # use: new_mesh = PCA_model.register_to_pca 

        moving_points = [cap.centroid for cap in model.left_ventricle.caps]
        moving_points += [cap.centroid for cap in model.right_ventricle.caps]

        # Extract apex coordinates from the left and right ventricles
        moving_points += [apex.xyz for apex in model.left_ventricle.apex_points]
        moving_points += [apex.xyz for apex in model.right_ventricle.apex_points]

        fix_points = self.anatomical_points

        # moving_model = model.mesh
        
        if registration_type == 1:
            transformed_mesh = self._rigid_transformation_scipy(fix_points, moving_points, model)
        elif registration_type == 2:
            transformed_mesh = self._rigid_transformation_landmark(fix_points, moving_points, model)
        else:
            raise ValueError("Unsupported registration type")

        transformed_mesh
        if save:
            transformed_mesh.save(filename)
        return transformed_mesh
    
    def display_model(self, show_torso: bool, show_electrodes: bool, show_anatomical_points: bool):
        import pyvista as pv

        plotter = pv.Plotter()

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


        if show_torso:
            body_points = np.array(list(self.body.values()))
            body_polydata = pv.PolyData(body_points)
            plotter.add_mesh(body_polydata, color='white', label='Torso')

        if show_electrodes:
            for name, index in self.electrodes.items():
                electrode_point = np.array(self.body[index])
                electrode_polydata = pv.PolyData(electrode_point)
                color = electrode_colors.get(name, 'grey')
                plotter.add_mesh(electrode_polydata, color=color, point_size=10, label=name)

        if show_anatomical_points:
            anatomical_point_polydata = pv.PolyData(np.array(self.anatomical_points))
            plotter.add_mesh(anatomical_point_polydata, color='green', point_size=10, label=name)

        plotter.add_legend()
        plotter.camera_position = 'xz'
        plotter.background_color = 'grey'
        plotter.show()
    