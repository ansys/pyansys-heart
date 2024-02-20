from ansys.heart.preprocessor.models.v0_2.models import HeartModel
import os

# environment variable indicating what model version to use.
# if the preprocessor generated a model with v0.2, then 
# this env variable needs to be set explicitely.
os.environ["ANSYS_HEART_MODEL_VERSION"]="v0.2"

# import relevant classes
from ansys.heart.simulator.settings.settings import DynaSettings
from ansys.heart.simulator.simulator import EPSimulator

# set path to ls-dyna, model and working directory.
lsdyna_path = "mppdyna_d_sse2_linux86_64_intelmmpi_105630"
lsdyna_path = r"D:\development\dyna-versions\daily_builds\daily_build_26102023\mppdyna_d_sse2_linux86_64_intelmmpi_105630"
path_to_model = r"D:\development\pyansys-heart\downloads\Rodero2021\01\FullHeart_v0_2\heart_model.pickle"
workdir = r"D:\development\pyansys-heart\downloads\Rodero2021\01\FullHeart_v0_2"

# initialize settings
dyna_settings = DynaSettings(lsdyna_path=lsdyna_path, dynatype="intelmpi", num_cpus=6, platform="wsl")

# load heart model
model = HeartModel.load_model(path_to_model)

# initialize simulator object
simulator = EPSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory=os.path.join(workdir, "ep-only"),
)

# load default simulation settings
simulator.settings.load_defaults()

# compute fiber orientation in the ventricles and atria
simulator.compute_fibers()

# Compute the conduction system
simulator.compute_purkinje()
# simulator.model.plot_purkinje()
simulator.compute_conduction_system()
simulator.model.plot_purkinje()


###############################################################################
# Start EP simulation
# ~~~~~~~~~~~~~~~~~~~

# store model with fibers and conduction system
simulator.model.dump_model(os.path.join(workdir, "heart_fib_beam.pickle"))

# start main simulation
simulator.simulate()