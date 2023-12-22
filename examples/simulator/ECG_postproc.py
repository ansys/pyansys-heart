"""
Post process EP simulation
--------------------------
This example shows you how to post process an EP simulation.
"""

###############################################################################
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules

# sphinx_gallery_start_ignore
# sphinx_gallery_thumbnail_path = '/_static/images/ep_post_activationtime.png'
# sphinx_gallery_end_ignore

import os
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator
from math import ceil 

from ansys.heart.postprocessor.EPpostprocessor import EPpostprocessor

workdir = Path(
    Path(__file__).resolve().parents[2], "downloads", "Strocchi2020", "01", "biventricle_scenario_2"
)

# Read DoE points
param_file_path = Path(Path(__file__).resolve().parent, "two_parameter_combinations150.csv")
parameters_df = pd.read_csv(param_file_path)

# ECG Visualization (code from Karim)
def _ax_plot(ax, x, y, secs=10, lwidth=0.5, amplitude_ecg=1.8, time_ticks=0.2):
    ax.set_xticks(np.arange(0, 11, time_ticks))    
    ax.set_yticks(np.arange(-ceil(amplitude_ecg), ceil(amplitude_ecg), 1.0))
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.set_ylim(-amplitude_ecg, amplitude_ecg)
    ax.set_xlim(0, secs)
    ax.grid(which='major', linestyle='-', linewidth='0.5', color='red')
    ax.grid(which='minor', linestyle='-', linewidth='0.5', color=(1, 0.7, 0.7))
    ax.plot(x, y, linewidth=lwidth)

def plot_1(ecg, sample_rate=500, title='ECG', fig_width=6, fig_height=2, line_w=1, ecg_amp=1.8, timetick=0.2, save_path=None):
    """Plot multi lead ECG chart.
    # Arguments
        ecg        : m x n ECG signal data, which m is number of leads and n is length of signal.
        sample_rate: Sample rate of the signal.
        title      : Title which will be shown on top off chart
        fig_width  : The width of the plot
        fig_height : The height of the plot
    """
    plt.figure(figsize=(fig_width,fig_height))
    plt.suptitle(title)
    plt.subplots_adjust(
        hspace = 0, 
        wspace = 0.04,
        left   = 0.04,  # the left side of the subplots of the figure
        right  = 0.98,  # the right side of the subplots of the figure
        bottom = 0.2,   # the bottom of the subplots of the figure
        top    = 0.88
        )
    seconds = len(ecg)/sample_rate

    ax = plt.subplot(1, 1, 1)
    #plt.rcParams['lines.linewidth'] = 5
    step = 1.0/sample_rate
    _ax_plot(ax,np.arange(0,len(ecg)*step,step),ecg, seconds, line_w, ecg_amp,timetick)

    if save_path:
        plt.savefig(save_path, dpi=300)

def extract_lead_II(file_path):
    '''Extract lead II signal potential from ECG results'''
    ECG = np.loadtxt(file_path, comments=['$', '*'], skiprows=4)[:, 0:11]
    VRA = ECG[:, 9]
    VLF = ECG[:, 8]
    Lead_II = VLF - VRA

    return Lead_II


# Initialization for output dataset
columns = ['LeadII', 'sigmaX', 'ratio2']
df = pd.DataFrame(columns=columns)

for index, row_data in parameters_df.iterrows():
    # scenario folder
    my_simulation_directory = os.path.join(
        workdir, 
        "simulation-EP",
        str(row_data["sigmaX"])[:4]
        + "_"
        + str(row_data["Ratio2"])[:4],
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
    # ecg result image
    lead_II_save_name = os.path.join(
        my_simulation_directory,
        "Lead_II.png",
    )

    # get activation times
    postproc = EPpostprocessor(results_path=ep_folder)
    activation_time_field = postproc.get_activation_times()
    # activation_time_field.plot(show_edges=False)
    
    activation_time_data = activation_time_field.data_as_list
    total_acctivation_time = max(activation_time_data) - min(activation_time_data)
    
    print(
        str(row_data["sigmaX"])[:4]
        + " _ "
        + str(row_data["Ratio2"])[:4]
        +"  Total activation time: " 
        + str(total_acctivation_time)
        +" ms")
    
    if total_acctivation_time>=40:
        # keep scenario
        # ECGs = postproc.read_ECGs(ecg_result_name)
        lead_II = extract_lead_II(ecg_result_name)
        plot_1(
            10*lead_II, 
            sample_rate = 1000, 
            line_w=1, 
            title = 'Lead II', 
            save_path=lead_II_save_name
        )
        sigmaX = str(row_data["sigmaX"])[:4]
        ratio2 = str(row_data["Ratio2"])[:4]
        
        # Create dataset for training model
        # data = {'lead_II': lead_II, 'sigmaX': sigmaX, 'ratio2': ratio2}
        df = df.append({'LeadII': lead_II, 'sigmaX': sigmaX, 'ratio2': ratio2}, ignore_index=True)

dataset_save_name = os.path.join(
    workdir,
    'dataset.csv',
)
df.to_csv(dataset_save_name, index=False)