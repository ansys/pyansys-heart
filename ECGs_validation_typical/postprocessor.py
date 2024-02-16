import os
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator
from math import ceil 

from ansys.heart.postprocessor.EPpostprocessor import EPpostprocessor

workdir = Path(
    Path(__file__).resolve().parents[2], "downloads", "Strocchi2020", "01", "Biv_typical_ecgs_15"
)

# Read DoE points
param_file_path = r'examples/SensitivityAnalysis_Sigma/scenarios/sigmaX.csv'
parameters_df = pd.read_csv(param_file_path)


def extract_lead_II(file_path):
    '''Extract lead II signal potential from ECG results'''
    ECG = np.loadtxt(file_path, comments=['$', '*'], skiprows=4)[:, 0:11]
    VRA = ECG[:, 9]
    VLF = ECG[:, 8]
    Lead_II = VLF - VRA

    return Lead_II


# Initialization for output dataset
columns = ['V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'sigmaX']
df = pd.DataFrame(columns=columns)

for index, row_data in parameters_df.iterrows():
    # scenario folder
    my_simulation_directory = os.path.join(
        workdir, 
        "simulation-EP\\Sigma",
        str(row_data["sigmaX"]),
    )
    # ep results folder
    ep_folder = os.path.join(
        my_simulation_directory,
        "d3plot",
    )
    # ecg results folder
    ecg_result_name = os.path.join(
        my_simulation_directory,
        "em_EKG_001.dat",
    )
    
    ECG = np.loadtxt(ecg_result_name, comments=['$', '*'], skiprows=4)[:, 0:11]
    # VRA = ECG[:, 7]
    # VLA = ECG[:, 8]
    # VRL = ECG[:, 9]
    # VLL = ECG[:, 10]

    # left -> right
    VLA = ECG[:, 7]
    VRA = ECG[:, 8]
    VLL = ECG[:, 9]
    VRL = ECG[:, 10]
    
    df = df.append(
        {
            'V1':ECG[:, 1], 
            'V2':ECG[:, 2], 
            'V3':ECG[:, 3], 
            'V4':ECG[:, 4], 
            'V5':ECG[:, 5], 
            'V6':ECG[:, 6], 
            'I': VLA - VRA,
            'II': VLL - VRA,
            'III': VLL - VLA,
            'aVR': (VRA - (VLA + VLL) / 2),
            'aVL': (VLA - (VRA + VLL) / 2),
            'aVF': (VLL - (VRA + VLA) / 2),
            'sigmaX': row_data["sigmaX"]
        }, 
        ignore_index=True
    )

dataset_save_name = os.path.join(
    workdir,
    'SA_sigma_12Lead.csv',
)
df.to_csv(dataset_save_name, index=False)



'''
VRA= ECG[:,9]
VLA= ECG[:,7]
VLF= ECG[:,8]
I  = VLA - VRA
II = VLF - VRA
III= VLF - VLA
aVR= VRA - (VLA+VLF)/2
aVL= VLA - (VLF+VRA)/2 
aVF= VLF - (VRA+VLA)/2
Vwct=(VLA+VRA+VLF)/3

'''