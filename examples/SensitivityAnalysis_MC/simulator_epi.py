import os
from pathlib import Path

from ansys.heart.preprocessor.mesh.objects import Point
import ansys.heart.preprocessor.models.v0_1.models as models
from ansys.heart.simulator.simulator import DynaSettings, EPSimulator

import numpy as np
import pandas as pd
from ansys.dyna.keywords import keywords
import ansys.heart.writer.dynawriter as writers

# set working directory and path to model.
workdir = Path(
    Path(__file__).resolve().parents[2], "downloads", "Strocchi2020", "01", "Biv_SA_MC_15"
)

path_to_model = os.path.join(workdir, "heart_model.pickle")

if not os.path.isfile(path_to_model):
    raise FileExistsError(f"{path_to_model} not found")

# specify LS-DYNA path (last tested working versions is intelmpi-linux-DEV-106117)
lsdyna_path = r"C:\Users\xuhu\lsdyna_smp_d_winx64\Other_version\mppdyna_d_winx64_msmpi\ls-dyna_mpp_d_Dev_104815-gc8c2d50328_winx64_ifort190_msmpi.exe"

if not os.path.isfile(lsdyna_path):
    raise FileExistsError(f"{lsdyna_path} not found.")

# load four chamber heart model.
model: models.FourChamber = models.HeartModel.load_model(path_to_model)

# Define electrode positions and add them to model (correspond to patient 01 only)
# Positions were defined using a template torso geometry.
electrodes = [
    Point(name="V1", xyz=[-29.893285751342773, 27.112899780273438, 373.30865478515625]),
    Point(name="V2", xyz=[33.68170928955078, 30.09606170654297, 380.5427551269531]),
    Point(name="V3", xyz=[56.33562469482422, 29.499839782714844, 355.533935546875]),
    Point(name="V4", xyz=[100.25729370117188, 43.61333465576172, 331.07635498046875]),
    Point(name="V5", xyz=[140.29800415039062, 81.36004638671875, 349.69970703125]),
    Point(name="V6", xyz=[167.9899139404297, 135.89862060546875, 366.18634033203125]),
    Point(name="RA", xyz=[-176.06332397460938, 57.632076263427734, 509.14202880859375]),
    Point(name="LA", xyz=[133.84518432617188, 101.44053649902344, 534.9176635742188]),
    Point(name="RL", xyz=[203.38825799615842, 56.19020893502452, 538.5052677637375]),
    Point(name="LL", xyz=[128.9441375732422, 92.85327911376953, 173.07363891601562]),
]
model.electrodes = electrodes


if not isinstance(model, models.FourChamber):
    raise TypeError("Expecting a FourChamber heart model.")

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

# instantiate simulator. Change options where necessary.
simulator = EPSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "simulation-EP"),
)

###############################################################################
# Load simulation settings
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Here we load the default settings.

simulator.settings.load_defaults()

simulator.compute_uhc()
simulator.compute_fibers()
simulator.compute_purkinje()
simulator.compute_conduction_system()


simulator.simulate()

param_file_path = r'D:\xuhu\pyansys-heart\examples\SensitivityAnalysis_MC\scenarios\percent_epi.csv'
parameters_df = pd.read_csv(param_file_path)



for index, row_data in parameters_df.iterrows():
    '''Write LS-DYNA file'''
    my_simulation_directory = os.path.join(
        workdir, 
        "simulation-EP\\MC_epi",
        str(row_data["percent_epi"]),
    )

    percent_epi = row_data["percent_epi"]

    folder_name = my_simulation_directory

    directory = os.path.join(simulator.root_directory, folder_name)

    # simulator._write_main_simulation_files(folder_name)
    '''code in _write_main_simulation_files'''

    export_directory = os.path.join(simulator.root_directory, folder_name)
    simulator.directories["main-ep"] = export_directory

    dyna_writer = writers.ElectrophysiologyDynaWriter(simulator.model, simulator.settings)
    

    # dyna_writer.update()
    '''Code in update'''
    dyna_writer._update_main_db()
    dyna_writer._update_solution_controls()
    dyna_writer._update_export_controls()

    dyna_writer._update_node_db()
    dyna_writer._update_parts_db()
    dyna_writer._update_solid_elements_db(add_fibers=True)

    dyna_writer._update_dummy_material_db()
    dyna_writer._update_ep_material_db()

    dyna_writer._update_segmentsets_db()

    dyna_writer._update_nodesets_db()


    # dyna_writer._update_cellmodels()
    '''code in _update_cellmodels'''
    for part in dyna_writer.model.parts:
        partname = part.name.lower()
        if ("atrium" in partname) or ("ventricle" in partname) or ("septum" in partname):
            ep_mid = part.pid
            # One cell model for myocardium, default value is epi layer parameters
            dyna_writer._add_Tentusscher_keyword_epi(matid=ep_mid)

    # different cell models for endo/mid/epi layer
    # percent_endo=0.17, percent_mid=0.41, percent_epi=0.42
    if "transmural" in dyna_writer.model.mesh.point_data.keys():
        (
            endo_id,
            mid_id,
            epi_id,
        ) = dyna_writer._create_myocardial_nodeset_layers(
            percent_endo=0.17, 
            percent_mid=0.58-percent_epi, 
            percent_epi=percent_epi
            )
        dyna_writer._add_Tentusscher_keyword_endo(matid=-endo_id)
        dyna_writer._add_Tentusscher_keyword_mid(matid=-mid_id)
        dyna_writer._add_Tentusscher_keyword_epi(matid=-epi_id)

    if dyna_writer.model.beam_network:
        # with smcoupl=1, coupling is disabled
        dyna_writer.kw_database.ep_settings.append(keywords.EmControlCoupling(smcoupl=1))
        dyna_writer._update_use_Purkinje()

    # update ep settings
    dyna_writer._update_ep_settings()
    dyna_writer._update_stimulation()

    if dyna_writer.model.info.add_blood_pool == True:
        dyna_writer._update_blood_settings()

    if hasattr(dyna_writer.model, "electrodes") and len(dyna_writer.model.electrodes) != 0:
        dyna_writer._update_ECG_coordinates()

    dyna_writer._get_list_of_includes()
    dyna_writer._add_includes()

    dyna_writer.export(export_directory)


    input_file = os.path.join(directory, "main.k")


    simulator._run_dyna(input_file)





    # simulator.simulate(folder_name=my_simulation_directory)