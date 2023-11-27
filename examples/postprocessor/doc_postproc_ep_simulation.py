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
import pathlib

from ansys.heart.postprocessor.EPpostprocessor import EPpostprocessor

# set ep results folder
ep_folder = os.path.join(
    pathlib.Path(__file__).absolute().parents[2],
    "downloads\\Strocchi2020\\01\FourChamber\\simulation-EP\\main-ep\\d3plot",
)
###############################################################################
# Instantiate the Postprocessor
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# instantiate Postprocessor

postproc = EPpostprocessor(results_path=ep_folder)

###############################################################################
# Activation times
# ~~~~~~~~~~~~~~~~
# Get activation times and plot the field

activation_time_field = postproc.get_activation_times()
activation_time_field.plot(show_edges=False)
###############################################################################
# .. image:: /_static/images/ep_post_activationtime.png
#   :width: 300pt
#   :align: center

# Compute total activation time
activation_time_data = activation_time_field.data_as_list
total_acctivation_time = max(activation_time_data) - min(activation_time_data)
print("Total activation time: " + str(total_acctivation_time) + " ms")

###############################################################################
# Transmembrane potentials
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Get transmembrane potentials on list of nodes and plot
vm, times = postproc.get_transmembrane_potential(node_id=[0, 1, 100], plot=True)
###############################################################################
# .. image:: /_static/images/ep_tm.png
#   :width: 300pt
#   :align: center

# Animate and export in vtk format
postproc.export_transmembrane_to_vtk()
# postproc.animate_transmembrane()

import numpy as np

electrodes = np.array(
    [
        [76.53798632905277, 167.67667039945263, 384.3139099410445],
        [64.97540262482013, 134.94983038904573, 330.4783062379255],
        [81.20629301587647, 107.06245851801455, 320.58645260857344],
        [85.04956217691463, 59.54502732121309, 299.2838953724169],
        [42.31377680589025, 27.997010728192166, 275.7620409440143],
        [-10.105919604515957, -7.176987485426985, 270.46379012676135],
        [-29.55095501940962, 317.12543912177983, 468.91891094294414],
        [-100.27895839242505, 135.64520460914244, 222.56688206809142],
        [203.38825799615842, 56.19020893502452, 538.5052677637375],
        [157.56391664248335, -81.66615972595032, 354.17867264210076],
    ]
)
ECGs, times = postproc.compute_ECGs(electrodes=electrodes)
print("finished")
