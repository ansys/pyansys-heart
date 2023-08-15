"""Example for running a FourChamber EP simulation from an existing preprocessed model."""
# import required modules
import os
from pathlib import Path
import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.simulator import EPSimulator


# set working directory and path to model.
workdir = Path(Path(__file__).parents[3], "downloads", "Strocchi2020", "01", "FourChamber")
path_to_model = os.path.join(workdir, "heart_model.pickle")

if not os.path.isfile(path_to_model):
    raise FileExistsError(f"{path_to_model} not found")

# specify LS-DYNA path
lsdyna_path = Path(Path(__file__).parents[5], "dyna-versions", "ls-dyna_smp_d.exe")

if not os.path.isfile(lsdyna_path):
    raise FileExistsError(f"{lsdyna_path} not found.")

# load four chamber heart model.
model: models.FourChamber = models.HeartModel.load_model(path_to_model)

if not isinstance(model, models.FourChamber):
    raise TypeError("Expecting a FourChamber heart model.")

# set base working directoy
model.info.workdir = str(workdir)

## instantiate simulator, please change the dynatype accordingly
simulator = EPSimulator(
    model=model,
    lsdynapath=lsdyna_path,
    dynatype="smp",
    num_cpus=1,
    simulation_directory=os.path.join(workdir, "simulation-EP"),
    platform="windows",
)

# load default simulation settings
simulator.settings.load_defaults()
# compute fiber orientation
simulator.compute_fibers()
# visualize computed fibers
simulator.model.plot_fibers(n_seed_points=2000)

# save image fibers
import pyvista as pv

doc_root = Path(Path(__file__).parents[3], "doc")
filename = Path(doc_root, "source\\images\\fibers.png")

mesh = simulator.model.mesh.ctp()
streamlines = mesh.streamlines(vectors="fiber", source_radius=75, n_points=2000)
tubes = streamlines.tube()
plotter = pv.Plotter(off_screen=True)
plotter.add_mesh(mesh, opacity=0.5, color="white")
plotter.add_mesh(tubes, color="white")
plotter.camera.roll = -60
plotter.screenshot(filename, transparent_background=False)

# compute purkinje network
simulator.compute_purkinje()
# compute the conduction system
simulator.compute_conduction_system()
# visualize purkinje
simulator.model.plot_purkinje()

# save purkinje image
import numpy as np

filename = Path(doc_root, "source\\images\\purkinje.png")
plotter = pv.Plotter(off_screen=True)
plotter.add_mesh(simulator.model.mesh, color="w", opacity=0.3)
for beams in simulator.model.mesh.beam_network:
    plotter.add_mesh(beams, color=np.random.uniform(size=3), line_width=3)
plotter.camera.roll = -60
plotter.screenshot(filename)

# start main ep-simulation
simulator.simulate()
