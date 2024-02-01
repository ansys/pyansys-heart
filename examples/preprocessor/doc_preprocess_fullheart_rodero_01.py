"""

Create a full heart model
---------------------------------
This example shows you how to download a case from the Rodero et al (2021) database
and process that case into a simulation-ready four chamber heart model.
"""

###############################################################################
# Example setup
# -------------
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory and generated model

# sphinx_gallery_start_ignore
# sphinx_gallery_thumbnail_path = '/_static/images/full_heart_mesh.png'
# sphinx_gallery_end_ignore

import os
from pathlib import Path

from ansys.heart.misc.downloader import download_case, unpack_case
import ansys.heart.preprocessor.models.v0_1.models as models

# sphinx_gallery_start_ignore
os.environ["USE_OLD_HEART_MODELS"] = "1"
# sphinx_gallery_end_ignore

# specify necessary paths.
# Note that we need to cast the paths to strings to facilitate serialization.
case_file = str(
    Path(Path(__file__).resolve().parents[2], "downloads", "Rodero2021", "01", "01.vtk")
)
download_folder = str(Path(Path(__file__).resolve().parents[2], "downloads"))
workdir = str(
    Path(Path(__file__).resolve().parents[2], "downloads", "Rodero2021", "01", "FullHeart2.0")
)

path_to_model = str(Path(workdir, "heart_model.pickle"))

###############################################################################
# Download the case
# ~~~~~~~~~~~~~~~~~
# Download and unpack the case from the public repository of full hearts if it is
# not yet available. This model will be used as input for the preprocessor.
if not os.path.isfile(case_file):
    path_to_downloaded_file = download_case(
        "Rodero2021", 1, download_folder=download_folder, overwrite=False
    )
    unpack_case(path_to_downloaded_file)

###############################################################################
# Set required information
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Set the right database to which this case belongs, and set other relevant
# information such as the desired mesh size.
info = models.ModelInfo(
    database="Rodero2021",
    path_to_case=case_file,
    work_directory=workdir,
    path_to_model=path_to_model,
    mesh_size=2.0,
)

# create the working directory
info.create_workdir()
# clean the working directory
info.clean_workdir(extensions_to_remove=[".stl", ".vtk", ".msh.h5"])
# dump information to stdout
info.dump_info()

###############################################################################
# Initialize the heart model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize the desired heart model, and invoke the main method to
# extract the simulation mesh and dump the model to disk. This will extract
# the relevant parts and features from the original model and remesh
# the entire surface and volume. Relevant anatomical features that are extracted include
# the endocardium, epicardium, apex, and valve annuli

# instantiate a four chamber model
model = models.FullHeart(info)

# extract the simulation mesh
model.extract_simulation_mesh()

# dump the model to disk for future use
model.dump_model(path_to_model)

# print information regarding the generated model
model.print_info()

# print part names
print(model.part_names)
# print volumes of some cavities:
print(f"Volume of LV cavity: {model.left_ventricle.cavity.volume} mm^3")
print(f"Volume of LV cavity: {model.left_atrium.cavity.volume} mm^3")

###############################################################################
# Visualize results
# ~~~~~~~~~~~~~~~~~
# You can visualize and inspect the components of the model by accessing
# various properties/attributes and invoke methods.

# plot the remeshed model
model.plot_mesh(show_edges=False)

###############################################################################
# .. image:: /_static/images/full_heart_mesh.png
#   :width: 400pt
#   :align: center

# plot the endocardial surface of the left ventricle.
model.left_ventricle.endocardium.plot(show_edges=True, color="r")

###############################################################################
# .. image:: /_static/images/full_heart_lv_endocardium.png
#   :width: 400pt
#   :align: center

# loop over all cavities and plot these in a single window.
import pyvista as pv

cavities = pv.PolyData()
for c in model.cavities:
    cavities += c.surface
cavities.plot(show_edges=True)

###############################################################################
# .. image:: /_static/images/full_heart_cavities.png
#   :width: 400pt
#   :align: center

# sphinx_gallery_start_ignore
# Generate static images for docs.
#
docs_images_folder = Path(Path(__file__).resolve().parents[2], "doc", "source", "_static", "images")

# Full mesh
filename = Path(docs_images_folder, "full_heart_mesh.png")
plotter = pv.Plotter(off_screen=True)
model.mesh.set_active_scalars("tags")
plotter.add_mesh(model.mesh, show_edges=False)
# plotter.show()
plotter.camera.roll = -60
plotter.screenshot(filename)

# Clipped full mesh

# left ventricle endocardium
filename = Path(docs_images_folder, "full_heart_lv_endocardium.png")
plotter = pv.Plotter(off_screen=True)
model.mesh.set_active_scalars(None)
plotter.add_mesh(model.left_ventricle.endocardium, color="r", show_edges=True)
plotter.camera.roll = -60
plotter.screenshot(filename)

# Cavities
filename = Path(docs_images_folder, "full_heart_cavities.png")
plotter = pv.Plotter(off_screen=True)
for c in model.cavities:
    plotter.add_mesh(c.surface, show_edges=True)
# plotter.show()
plotter.camera.roll = -60
plotter.screenshot(filename)
# sphinx_gallery_end_ignore
