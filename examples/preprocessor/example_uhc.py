"""

UHC example
--------------------
This example shows how to compute universal heart coordinate for a BiVentricle heart model.
"""
###############################################################################
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory, model, and ls-dyna executable.

# sphinx_gallery_start_ignore
# Note that we need to put the thumbnail here to avoid weird rendering in the html page.
# sphinx_gallery_thumbnail_path = '_static/images/thumbnails/uvc.png'
# sphinx_gallery_end_ignore

import copy
import os
from pathlib import Path

import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.simulator import BaseSimulator
import pyvista as pv

# set working directory and path to model.
workdir = Path(Path(__file__).parents[2], "downloads", "Strocchi2020", "01", "FourChamber")

path_to_model = os.path.join(workdir, "heart_model.pickle")

if not os.path.isfile(path_to_model):
    raise FileExistsError(f"{path_to_model} not found")

# specify LS-DYNA path
lsdyna_path = Path(Path(__file__).parents[4], "dyna-versions", "ls-dyna_smp_d.exe")

if not os.path.isfile(lsdyna_path):
    raise FileExistsError(f"{lsdyna_path} not found.")

# load heart model.
model: models.BiVentricle = models.HeartModel.load_model(path_to_model)

# set base working directory
model.info.workdir = str(workdir)

###############################################################################
# Instantiate the simulator object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# instantiate simulator. Change options where necessary.

simulator = BaseSimulator(
    model=model,
    lsdynapath=lsdyna_path,
    dynatype="smp",
    num_cpus=1,
    simulation_directory=os.path.join(workdir, "simulation"),
    platform="windows",
)

###############################################################################
# Compute UHC
# ~~~~~~~~~~~
# Compute UHC using Laplace Dirichlet method.

simulator.compute_uhc()

###############################################################################
# .. note::
#    There are several definitions for UHC (see https://github.com/KIT-IBT/Cobiveco).
#    Here, a simple approach is taken and the
#    Dirichlet conditions are shown below. At rotational direction, the start (pi), end (-pi)
#    and middle (0) points are defined from 4CV long axis cut view.

###############################################################################
# .. image:: /_static/images/uvc_bc.png
#   :width: 600pt
#   :align: center

###############################################################################
# Visualization
# ~~~~~~~~~~~~~

data = pv.read(os.path.join(workdir, "simulation", "uhc", "uhc.vtk"))

plotter = pv.Plotter(shape=(1, 3))

plotter.subplot(0, 0)
plotter.add_mesh(data, scalars="longitudinal")

plotter.subplot(0, 1)
plotter.add_mesh(copy.copy(data), scalars="transmural")

plotter.subplot(0, 2)
plotter.add_mesh(copy.copy(data), scalars="rotational")
plotter.show()

###############################################################################
# .. image:: /_static/images/uvc_result.png
#   :width: 600pt
#   :align: center
