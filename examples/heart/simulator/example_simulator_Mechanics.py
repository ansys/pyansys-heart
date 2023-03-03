"""Example to preprocess and run a biventruclar mechanics simulation."""
import os
import pathlib

import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.simulator import MechanicsSimulator
from ansys.heart.simulator.support import run_preprocessor

# import ansys.heart.writer.dynawriter as writers

if __name__ == "__main__":
    """BiVentricle example.

    1. Extracts simulation mesh
    2. Use simulator class to launch Mechanics simulation
            - compute fibers
            - compute stress-free configuration
            - launch main simulation

    Please change paths according to your workspace
    """
    # extract simulation mesh(es)
    path_to_case = os.path.join(
        pathlib.Path(__file__).parents[3], "downloads\\Strocchi2020\\01\\01.case"
    )
    workdir = os.path.join(pathlib.Path(path_to_case).parent, "BiVentricle")

    path_to_model = os.path.join(workdir, "heart_model.pickle")

    use_preprocessor = False

    if use_preprocessor:
        model = run_preprocessor(
            model_type=models.BiVentricle,
            database="Strocchi2020",
            path_original_mesh=path_to_case,
            work_directory=workdir,
            path_to_model=path_to_model,
            mesh_size=2.0,
        )

    # specify LS-DYNA path
    lsdyna_path = r"D:\my_path_to_ls_dyna\lsdyna_executable.exe"

    # BASIC examples with defaults:

    # Load model (e.g. when you skip the preprocessor):
    model: models.BiVentricle = models.HeartModel.load_model(path_to_model)
    ## instantiate simulator, please change the dynatype accordingly
    simulator = MechanicsSimulator(
        model=model,
        lsdynapath=lsdyna_path,
        dynatype="smp",
        num_cpus=1,
        simulation_directory=os.path.join(model.info.workdir, "simulation-mechanics"),
    )

    # load default settings.
    simulator.settings.load_defaults()
    # compute the fiber orientation
    simulator.compute_fibers()
    # visualize computed fibers
    simulator.model.plot_fibers(n_seed_points=2000)
    # compute the stress free configuration
    simulator.compute_stress_free_configuration()
    # do the main simulation
    simulator.simulate()
    print("done")
