import os
from pathlib import Path

from ansys.heart.preprocessor.mesh.objects import Point
import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.simulator import DynaSettings, EPSimulator, run_lsdyna
import pandas as pd
import pyvista as pv
import numpy as np
# __file__ = r"C:\Users\xuhu\pyheart-lib\examples\simulator\test.ipynb"

# set working directory and path to model.
workdir = Path(
    Path(__file__).resolve().parents[2], "downloads", "Strocchi2020", "01", "biventricle_scenario_2"
)

path_to_model = os.path.join(workdir, "heart_model.pickle")

if not os.path.isfile(path_to_model):
    raise FileExistsError(f"{path_to_model} not found")

# specify LS-DYNA path (last tested working versions is DEV-104399)
lsdyna_path = r"C:\Users\xuhu\lsdyna_smp_d_winx64\Other_version\mppdyna_d_winx64_msmpi\ls-dyna_mpp_d_Dev_104815-gc8c2d50328_winx64_ifort190_msmpi.exe"

if not os.path.isfile(lsdyna_path):
    raise FileExistsError(f"{lsdyna_path} not found.")

# load BiVentricle heart model.
model: models.BiVentricle = models.HeartModel.load_model(path_to_model)

# Define electrode positions and add them to model
electrodes = [
    Point(name="V1", xyz=[-29.893285751342773, 27.112899780273438, 373.30865478515625]),
    Point(name="V2", xyz=[33.68170928955078, 30.09606170654297, 380.5427551269531]),
    Point(name="V3", xyz=[56.33562469482422, 29.499839782714844, 355.533935546875]),
    Point(name="V4", xyz=[100.25729370117188, 43.61333465576172, 331.07635498046875]),
    Point(name="V5", xyz=[140.29800415039062, 81.36004638671875, 349.69970703125]),
    Point(name="V6", xyz=[167.9899139404297, 135.89862060546875, 366.18634033203125]),
    Point(name="RA", xyz=[-176.06332397460938, 57.632076263427734, 509.14202880859375]),
    Point(name="LA", xyz=[133.84518432617188, 101.44053649902344, 534.9176635742188]),
    Point(name="RL", xyz=[-103.60343933105469, 64.02100372314453, 160.0018310546875]),
    Point(name="LL", xyz=[128.9441375732422, 92.85327911376953, 173.07363891601562]),
]
model.electrodes = electrodes


if not isinstance(model, models.BiVentricle):
    raise TypeError("Expecting a BiVentricle heart model.")

# set base working directory
model.info.workdir = str(workdir)

###############################################################################
# Instantiate the simulator object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# instantiate the simulator and settings appropriately.

# instantaiate dyna settings of choice
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path,
    dynatype="msmpi",
    num_cpus=1,
)

simulator = EPSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "simulation-EP"),
)

simulator.compute_uvc()
simulator.compute_fibers()
simulator.compute_purkinje()
# simulator.compute_conduction_system()

simulator.model.dump_model(path_to_model)

param_file_path = Path(Path(__file__).resolve().parent, "two_parameter_combinations150.csv")
parameters_df = pd.read_csv(param_file_path)

for index, row_data in parameters_df.iterrows():
    '''Write LS-DYNA file'''
    my_simulation_directory = os.path.join(
        workdir, 
        "simulation-EP",
        str(row_data["sigmaX"])[:4]
        + "_"
        + str(row_data["Ratio2"])[:4],
    )
    '''
    sigmaX values: (0.2, 2)
    ratio2 (in keywords.EmMat001): [1,10] (such that: sigmaPurkinje=ratio2 * sigmaX)
    '''
    simulator.settings.electrophysiology.sigma11=float(row_data["sigmaX"])
    simulator.settings.purkinje.sigma=float(row_data["sigmaX"])*float(row_data["Ratio2"])

    simulator._write_main_simulation_files(my_simulation_directory)

# Change CPU number
simulator.dyna_settings.num_cpus=20

for index, row_data in parameters_df.iterrows():
    '''Run LS-DYNA file'''
    directory = os.path.join(
        workdir, 
        "simulation-EP",
        str(row_data["sigmaX"])[:4]
        + "_"
        + str(row_data["Ratio2"])[:4],
    )
    input_file = os.path.join(directory, "main.k")
    simulator._run_dyna(input_file)