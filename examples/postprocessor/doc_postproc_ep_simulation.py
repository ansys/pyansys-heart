# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
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

from ansys.heart.postprocessor.ep_postprocessor import EPpostprocessor

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
# 12-LEAD ECGs
# ~~~~~~~~~~~~~~~~
# Plot 12-Lead ECGs

ECGs, times = postproc.read_ECGs(
    os.path.join(
        pathlib.Path(__file__).absolute().parents[2],
        "downloads\\Strocchi2020\\01\FourChamber\\simulation-EP\\main-ep\\em_EKG_001.dat",
    )
)


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
