"""Example to pre-process data from Strocchi2020 and Cristobal2021."""
import os
from pathlib import Path

import ansys.heart.preprocessor.models as models
import ansys.heart.writer.dynawriter as writers
from ansys.heart.misc.downloader import download_case, unpack_case
from ansys.heart.simulator.simulator import run_lsdyna
import pyvista
import numpy as np
import vtk

from vtkmodules.vtkCommonDataModel import vtkIterativeClosestPointTransform
from ansys.heart.simulator.simulator import DynaSettings, EPSimulator

os.environ["USE_OLD_HEART_MODELS"] = "1"

# __file__ = r"c:\Users\xuhu\pyheart-lib\examples\preprocessor\doc_ECG_coordinates.py"

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
    [81.90321388, 57.90000882, 205.76663367], # mitral-valve
    [94.35242091, 75.99022315, 213.31654731], # aortic-valve
    [67.14045655, 102.49380179, 216.20654707], # tricuspid-valve
    [121.58912558, 89.76497459, 223.29557159], # pulmonary-valve
    [70.87069056682236, 84.83837198547876, 295.6765864478138], # left endo
    [70.54655746919204, 84.50457846174797, 297.2737993295601], # left epi
    [76.04229182019685, 66.53094359081156, 297.7182142431582], # right endo
    [75.08606835375224, 66.33759424571653, 302.2811669120656], # right epi   
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


transformed_electrodes = model.define_ECG_coordinates(move_points, electrode_positions)

write_lsdyna_files = True


if write_lsdyna_files:
    for writer in (
        writers.ElectrophysiologyDynaWriter(model),
    ):

        writer.update()

        writer.export_databases(
            os.path.join(writer.model.info.workdir, "ECG")
        )


# specify LS-DYNA path
# lsdyna_path = r"C:\temporary\1\test\ls-dyna_smp_d_Dev_103633-gcc846f4b4e_winx64_ifort190.exe"
lsdyna_path = r"C:\Users\xuhu\lsdyna_smp_d_winx64\Other_version\mppdyna_d_winx64_msmpi\ls-dyna_mpp_d_Dev_104815-gc8c2d50328_winx64_ifort190_msmpi.exe"

if not os.path.isfile(lsdyna_path):
    raise FileExistsError(f"{lsdyna_path} not found.")

dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path,
    dynatype="smp",
    num_cpus=1,
)

run_lsdyna(
    path_to_input = os.path.join(
        writer.model.info.workdir, 
        "ECG", 
        "main.k"
    ),
    settings = dyna_settings,
    simulation_directory = os.path.join(workdir, "ECG"),
)

