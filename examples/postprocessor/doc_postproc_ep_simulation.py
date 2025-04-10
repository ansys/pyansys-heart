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
Postprocess an EP simulation
-----------------------------
This example shows how to postprocess an EP simulation.
"""

###############################################################################
# Perform required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~

# sphinx_gallery_start_ignore
# sphinx_gallery_thumbnail_path = '_static/images/ep_post_activationtime.png'
# sphinx_gallery_end_ignore

import pathlib

from ansys.health.heart.post.dpf_utils import EPpostprocessor

# set ep results folder
ep_folder = (
    pathlib.Path.home()
    / "pyansys-heart"
    / "downloads"
    / "Strocchi2020"
    / "01"
    / "FourChamber"
    / "simulation-EP"
    / "main-ep"
    / "d3plot"
)
###############################################################################
# Instantiate the postprocessor
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

postproc = EPpostprocessor(results_path=ep_folder)


###############################################################################
# Plot 12-LEAD ECGs
# ~~~~~~~~~~~~~~~~~

path_to_ecg_file = ep_folder.parent / "em_EKG_001.dat"

ECGs, times = postproc.read_ECGs(path_to_ecg_file)


ECGs12 = postproc.compute_12_lead_ECGs(ECGs=ECGs, times=times, plot=True)

###############################################################################
# .. image:: /_static/images/ep_post_12LeadECGs.png
#   :width: 300pt
#   :align: center

###############################################################################
# Get activation times
# ~~~~~~~~~~~~~~~~~~~~
# Get activation times and plot the field.

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
# Get transmembrane potentials
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get transmembrane potentials on a list of nodes and plot.

vm, times = postproc.get_transmembrane_potential(node_id=[0, 1, 100], plot=True)
###############################################################################
# .. image:: /_static/images/ep_tm.png
#   :width: 300pt
#   :align: center

# Animate and export in vtk format
postproc.export_transmembrane_to_vtk()
postproc.animate_transmembrane()
