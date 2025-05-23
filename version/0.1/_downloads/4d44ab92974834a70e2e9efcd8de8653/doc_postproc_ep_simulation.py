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
postproc.animate_transmembrane()
