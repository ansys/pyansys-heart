# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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
# sphinx_gallery_thumbnail_path = '_static/images/ep_post_activationtime.png'
# sphinx_gallery_end_ignore

import os
import pathlib

from ansys.heart.postprocessor.EPpostprocessor import EPpostprocessor

os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"
# set ep results folder
ep_folder = os.path.join(
    r"D:\REPOS\pyansys-heart\downloads\Strocchi2020\01\FourChamber\simulation-EP-origin\main-ep-ReactionEikonal\d3plot"
)
###############################################################################
# Instantiate the Postprocessor
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# instantiate Postprocessor

postproc = EPpostprocessor(results_path=ep_folder)


###############################################################################
# 12-LEAD ECGs
# ~~~~~~~~~~~~~~~~
# Plot 12-Lead ECGs


ECGsdyna, times = postproc.read_ECGs(
    r"D:\REPOS\pyansys-heart\downloads\Strocchi2020\01\FourChamber\simulation-EP-origin\main-ep-ReactionEikonal\em_EKG_001.dat"
)
import numpy as np

electrodes = np.array(
    [
        [-29.893285751342773, 27.112899780273438, 373.30865478515625],
        [33.68170928955078, 30.09606170654297, 380.5427551269531],
        [56.33562469482422, 29.499839782714844, 355.533935546875],
        [100.25729370117188, 43.61333465576172, 331.07635498046875],
        [140.29800415039062, 81.36004638671875, 349.69970703125],
        [167.9899139404297, 135.89862060546875, 366.18634033203125],
        [-176.06332397460938, 57.632076263427734, 509.14202880859375],
        [133.84518432617188, 101.44053649902344, 534.9176635742188],
        [203.38825799615842, 56.19020893502452, 538.5052677637375],
        [128.9441375732422, 92.85327911376953, 173.07363891601562],
    ]
)
ECGs12dyna = postproc.compute_12_lead_ECGs(ECGs=ECGsdyna, times=times, plot=True)
ECGs, times = postproc.compute_ECGs(electrodes)
ECGs12 = postproc.compute_12_lead_ECGs(ECGs=ECGs, times=times, plot=True)

###############################################################################
# .. image:: /_static/images/ep_post_12LeadECGs.png
#   :width: 300pt
#   :align: center

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
postproc.animate_transmembrane()
