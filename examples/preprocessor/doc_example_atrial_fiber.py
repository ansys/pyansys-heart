"""
Atrial fiber
------------
This examples shows how to generate fibers with Laplace-Dirichlet-Rule-Based-Method
"""
###############################################################################
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory, model, and ls-dyna executable (uses DEV-104373-g6d20c20aee).

# sphinx_gallery_start_ignore
# Note that we need to put the thumbnail here to avoid weird rendering in the html page.
# sphinx_gallery_thumbnail_path = '_static/images/thumbnails/atrial_fiber.png'
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
lsdyna_path = r"ls-dyna_smp"

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
    lsdyna_path=lsdyna_path,
    dynatype="smp",
    num_cpus=1,
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
# Compute atrial fibers
# ~~~~~~~~~~~~~~~~~~~~~

# Compute left atrium fiber
la = simulator.compute_left_atrial_fiber()

# Appendage apex point should be manually given to compute right atrium fiber
appendage_apex = [-50, 106, 425]
ra = simulator.compute_right_atrial_fiber(appendage_apex)

###############################################################################
# Plot bundle selection results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
la.set_active_scalars("bundle")
la.plot()

###############################################################################
# .. image:: /_static/images/la_bundle.png
#   :width: 400pt
#   :align: center

ra.set_active_scalars("bundle")
ra.plot()

###############################################################################
# .. image:: /_static/images/ra_bundle.png
#   :width: 400pt
#   :align: center

###############################################################################
# Plot fibers
# ~~~~~~~~~~~
plotter = pv.Plotter()
mesh = la.ctp()
streamlines = mesh.streamlines(vectors="e_l", source_radius=50, n_points=50000)
tubes = streamlines.tube()
plotter.add_mesh(mesh, opacity=0.5, color="white")
plotter.add_mesh(tubes, color="red")
plotter.show()

###############################################################################
# .. image:: /_static/images/la_fiber.png
#   :width: 400pt
#   :align: center

plotter = pv.Plotter()
mesh = ra.ctp()
streamlines = mesh.streamlines(vectors="e_l", source_radius=50, n_points=50000)
tubes = streamlines.tube()
plotter.add_mesh(mesh, opacity=0.5, color="white")
plotter.add_mesh(tubes, color="red")
plotter.show()

###############################################################################
# .. image:: /_static/images/ra_fiber.png
#   :width: 400pt
#   :align: center

###############################################################################

# Atrial fibers are automatically assigned to heart model after computation.

plotter = pv.Plotter()
mesh = model.mesh.ctp()
streamlines = mesh.streamlines(vectors="fiber", source_radius=100, n_points=50000)
tubes = streamlines.tube()
plotter.add_mesh(mesh, opacity=0.5, color="white")
plotter.add_mesh(tubes, color="red")
plotter.show()

###############################################################################
# .. image:: /_static/images/atrial_fiber_assign.png
#   :width: 400pt
#   :align: center
