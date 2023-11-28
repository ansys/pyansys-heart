import os
from pathlib import Path

from ansys.heart.misc.downloader import download_case, unpack_case
import ansys.heart.preprocessor.models as models
from ansys.heart.preprocessor.mesh.objects import Point

import pyvista
import numpy as np
import vtk
import sys

import ansys.heart.writer.dynawriter as writers
from vtkmodules.vtkCommonDataModel import vtkIterativeClosestPointTransform

from ansys.heart.simulator.simulator import DynaSettings, run_lsdyna


os.environ["USE_OLD_HEART_MODELS"] = "1"

# __file__ = r"c:\Users\xuhu\pyheart-lib\examples\preprocessor\doc_ECG_coordinates,py"

case_file = str(
    Path(Path(__file__).resolve().parents[2], "downloads", "Strocchi2020", "01", "01.case")
)
download_folder = str(Path(Path(__file__).resolve().parents[2], "downloads"))
workdir = str(
    Path(Path(__file__).resolve().parents[2], "downloads", "Strocchi2020", "01", "Biv")
)
# path_to_model = str(Path(workdir, str(sys.argv[1])+"heart_model.pickle"))
path_to_model = str(Path(workdir, "heart_model.pickle"))

print('11111111111111111111111')

if not os.path.isfile(case_file):
    path_to_downloaded_file = download_case(
        "Strocchi2020", 1, download_folder=download_folder, overwrite=False
    )
    unpack_case(path_to_downloaded_file)

# m_size = float(sys.argv[1])
# print('MYYYYYYYYYYY mesh : ', mesh_size)
print('22222222222222222222')

info = models.ModelInfo(
    database="Strocchi2020",
    path_to_case=case_file,
    work_directory=workdir,
    path_to_model=path_to_model,
    add_blood_pool=False,
    mesh_size=0.25,
)

print('3333333333333333333333')
# create the working directory
info.create_workdir()
# clean the working directory
info.clean_workdir(extensions_to_remove=[".stl", ".vtk", ".msh.h5"])
# dump information to stdout
info.dump_info()
print('4444444444444444')

# instantiate a four chamber model
model = models.FourChamber(info)
print('555555555555555555555')

# # extract the simulation mesh
# model.extract_simulation_mesh()

# # dump the model to disk for future use
# model.dump_model(path_to_model)
# # print the resulting information
# model.print_info()

# C:\Users\xuhu\pyheart-lib\examples\simulator\doc_ECG_scenario.py

# specify LS-DYNA path (last tested working versions is DEV-104399)
lsdyna_path = r"C:\Users\xuhu\lsdyna_smp_d_winx64\Other_version\mppdyna_d_winx64_msmpi\ls-dyna_mpp_d_Dev_104815-gc8c2d50328_winx64_ifort190_msmpi.exe"

if not os.path.isfile(lsdyna_path):
    raise FileExistsError(f"{lsdyna_path} not found.")
print('666666666666666666')

print('77777777777777777777')

# Define electrode positions and add them to model
electrodes = [
    Point(name="V1", xyz=[76.53798632905277, 167.67667039945263, 384.3139099410445]),
    Point(name="V2", xyz=[64.97540262482013, 134.94983038904573, 330.4783062379255]),
    Point(name="V3", xyz=[81.20629301587647, 107.06245851801455, 320.58645260857344]),
    Point(name="V4", xyz=[85.04956217691463, 59.54502732121309, 299.2838953724169]),
    Point(name="V5", xyz=[42.31377680589025, 27.997010728192166, 275.7620409440143]),
    Point(name="V6", xyz=[-10.105919604515957, -7.176987485426985, 270.46379012676135]),
    Point(name="RA", xyz=[-29.55095501940962, 317.12543912177983, 468.91891094294414]),
    Point(name="LA", xyz=[-100.27895839242505, 135.64520460914244, 222.56688206809142]),
    Point(name="RL", xyz=[203.38825799615842, 56.19020893502452, 538.5052677637375]),
    Point(name="LL", xyz=[157.56391664248335, -81.66615972595032, 354.17867264210076]),
]
model.electrodes = electrodes
print('8888888888888888888888:', model.electrodes)

if not isinstance(model, models.FourChamber):
    raise TypeError("Expecting a FourChamber heart model.")

# set base working directory
model.info.workdir = str(workdir)
print('99999999999999999999999')

write_lsdyna_files = True

print('1212121212')

if write_lsdyna_files:
    for writer in (
        writers.ElectrophysiologyDynaWriter(model),
        print('5552222')
    ):
        writer.update()
        print('writer.update() pass')

        writer.export_databases(
            os.path.join(writer.model.info.workdir, "ECG")
        )
        print('export_databases pass')


print('1313131313131')

# specify LS-DYNA path
# lsdyna_path = r"C:\temporary\1\test\ls-dyna_smp_d_Dev_103633-gcc846f4b4e_winx64_ifort190.exe"


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

# C:\Users\xuhu\pyheart-lib\examples\simulator\doc_ECG_scenario.py