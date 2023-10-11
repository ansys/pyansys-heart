"""

Post process EP simulation
--------------------------
This example shows you how to use post process an EP simulation.
"""

###############################################################################
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules

# sphinx_gallery_start_ignore
# Note that we need to put the thumbnail here to avoid weird rendering in the html page.
# sphinx_gallery_thumbnail_path = '_static/images/thumbnails/klotz.png'
# sphinx_gallery_end_ignore
# import os

from ansys.heart.postprocessor.EPpostprocessor import EPpostprocessor

# from ansys.heart.postprocessor.dpf_utils import D3plotReader
# import ansys.heart.preprocessor.models as models

###############################################################################
# Set relevant paths
# ~~~~~~~~~~~~~~~~~~

# path_to_model = r"D:\REPOS\pyheart-lib\downloads\Strocchi2020\01\FourChamber\heart_model.pickle"

# if not os.path.isfile(path_to_model):
#     raise FileExistsError(f"{path_to_model} not found")

# load heart model.
# model: models.FourChamber = models.FourChamber.load_model(path_to_model)

# set zerop simulation path
ep_folder = (
    r"D:\REPOS\pyheart-lib\downloads\Strocchi2020\01\FourChamber\simulation-EP\main-ep\d3plot"
)


postproc = EPpostprocessor(results_path=ep_folder)
activation_time_field = postproc.get_activation_times()
activation_time_field.plot()

start_time = 0
endtime = 100


###############################################################################
# Run default process scripts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a folder "post" will be created with key simulation results (json, png, vtk...)
# postprocessor = D3plotReader(ep_folder)
# fields = postprocessor.get_epfields()
# fields[10].plot()
# print("finished")
