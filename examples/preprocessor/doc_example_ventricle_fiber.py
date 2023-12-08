"""
Ventricle fiber
------------
This examples shows how to generate ventricle fibers with Laplace-Dirichlet-Rule-Based-Method
"""
###############################################################################
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory, model, and ls-dyna executable (uses > DEV-104373-g6d20c20aee).

# sphinx_gallery_start_ignore
# Note that we need to put the thumbnail here to avoid weird rendering in the html page.
# sphinx_gallery_thumbnail_path = '_static/images/lsdyna_fibers.png'
# sphinx_gallery_end_ignore

import os
from pathlib import Path

import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.simulator import BaseSimulator, DynaSettings
import numpy as np
import pyvista as pv

# set working directory and path to model.
workdir = Path(
    Path(__file__).resolve().parents[2], "downloads", "Strocchi2020", "01", "FourChamber"
)

path_to_model = os.path.join(workdir, "heart_model.pickle")

if not os.path.isfile(path_to_model):
    raise FileExistsError(f"{path_to_model} not found")

# specify LS-DYNA path
lsdyna_path = r"."

if not os.path.isfile(lsdyna_path):
    raise FileExistsError(f"{lsdyna_path} not found.")

# load heart model.
model: models.FourChamber = models.HeartModel.load_model(path_to_model)

# set base working directory
model.info.workdir = str(workdir)

###############################################################################
# Instantiate the simulator object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# instantiate simulator. Change options where necessary.

# instantaiate dyna settings of choice
dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path, dynatype="intelmpi", num_cpus=4, platform="wsl"
)

simulator = BaseSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "simulation"),
)

# remove fiber/sheet information if already exists
model.mesh.cell_data["fiber"] = np.zeros((model.mesh.n_cells, 3))
model.mesh.cell_data["sheet"] = np.zeros((model.mesh.n_cells, 3))

###############################################################################
# Compute ventricle fibers use LS-DYNA implemented method
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This method is a simplified version of doi:10.1007/s10439-012-0593-5 without interpolation


simulator.compute_fibers(method="LSDYNA")  # this is also the default option

plotter = pv.Plotter()
mesh = model.mesh.ctp()
streamlines = mesh.streamlines(vectors="fiber", source_radius=100, n_points=50000)
tubes = streamlines.tube()
plotter.add_mesh(mesh, opacity=0.5, color="white")
plotter.add_mesh(tubes, color="red")
plotter.show()

###############################################################################
# .. image:: /_static/images/lsdyna_fibers.png
#   :width: 400pt
#   :align: center


model.dump_model(os.path.join(workdir, "with_lsdyna_fibers.pickle"))

###############################################################################
# Compute ventricle fibers use D-RBM
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This method is described in https://doi.org/10.1016/j.cma.2020.113468

# Define rotation angle for fibers, from endo to epi
rotation_angles = {
    "left alpha": [-60, 60],
    "right alpha": [90, -25],
    "left beta": [-20, 20],
    "right beta": [0, 20],
    "outflow alpha": [90, 0],  # special parameters for D-RBM
    "outflow beta": [0, 0],  # special parameters for D-RBM
}
simulator.compute_fibers(method="D-RBM")

plotter = pv.Plotter()
mesh = model.mesh.ctp()
streamlines = mesh.streamlines(vectors="fiber", source_radius=100, n_points=50000)
tubes = streamlines.tube()
plotter.add_mesh(mesh, opacity=0.5, color="white")
plotter.add_mesh(tubes, color="red")
plotter.show()

###############################################################################
# .. image:: /_static/images/d-rbm_fibers.png
#   :width: 400pt
#   :align: center

model.dump_model(os.path.join(workdir, "with_drbm_fibers.pickle"))
