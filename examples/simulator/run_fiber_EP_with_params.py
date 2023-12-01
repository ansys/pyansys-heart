import os
import pandas as pd

from pathlib import Path

from ansys.heart.preprocessor.mesh.objects import Point
import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.simulator import DynaSettings, EPSimulator
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


param_file_path = r'C:\Users\xuhu\pyheart-lib\examples\simulator\parameter_combinations.csv'
parameters_df = pd.read_csv(param_file_path)

for index, row in parameters_df.iterrows():
    # instantiate simulator. Change options where necessary.
    simulator = EPSimulator(
        model=model,
        dyna_settings=dyna_settings,
        simulation_directory=os.path.join(workdir, "simulation-EP"),
    )

    simulator.settings.load_defaults()

    simulator.settings.load_with_EP_params(
        Purkinje_edgelen=row['Purkinje_edgelen'], 
        sigmaX=row['SigmaX'], 
        ratio=row['Ratio'], 
        ratio2=row['Ratio2']
    )

    print(simulator.settings.electrophysiology.sigma11)
    print(simulator.settings.electrophysiology.sigma22)
    print(simulator.settings.electrophysiology.sigma33)
    print(simulator.settings.purkinje.edgelen)
    print(simulator.settings.purkinje.sigma)


    # # simulator.compute_fibers()
    # # simulator.model.plot_fibers(n_seed_points=2000)
    # # print('sucesss 5555555555')

    simulator.compute_purkinje()
    
    simulator.compute_conduction_system()

    # simulator.model.plot_purkinje()

    simulator.simulate()