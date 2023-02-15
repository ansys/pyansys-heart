"""Example to pre-process data from Strocchi2020 and Rodero2021."""
import os
import pathlib

import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.simulator import EPSimulator, MechanicsSimulator
from ansys.heart.simulator.settings.settings import SimulationSettings
from ansys.heart.simulator.support import run_preprocessor

# import ansys.heart.writer.dynawriter as writers

if __name__ == "__main__":
    """Full Heart example.

    1. Extracts simulation mesh
    2. Use simulator class to launch
        1. EP-simulation
            - compute fibers
            - compute purkinje network
            - launch main simulation
        2. Mechanics-simulation
            - compute fibers
            - compute stress-free configuration
            - launch main simulation

    Please change paths
    """
    # extract simulation mesh(es)
    path_to_case = os.path.join(
        pathlib.Path(__file__).parents[3], "downloads\\Strocchi2020\\01\\01.case"
    )
    workdir = os.path.join(pathlib.Path(path_to_case).parent, "BiVentricleTestSimulator1")

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
    lsdyna_path = (
        r"D:\development\dyna-versions\ls-dyna_smp_d_Dev_97584-g1b99fd817b_winx64_ifort190.exe"
    )

    # BASIC examples with defaults:

    # Load model (e.g. when you skip the preprocessor):
    model: models.BiVentricle = models.HeartModel.load_model(path_to_model)
    ## instantiate simulator.
    simulator = EPSimulator(model, lsdyna_path, "smp", num_cpus=4)
    # load default settings
    simulator.settings.load_defaults()
    # compute the fiber orientation
    simulator.compute_fibers()
    # visualize computed fibers
    simulator.model.plot_fibers(n_seed_points=2000)
    # compute purkinje network
    simulator.compute_purkinje()
    # visualize purkinje
    simulator.model.plot_purkinje()
    # start simulation
    simulator.simulate()

    # ** Mechanics simulation (simplified EP and fluids)

    # Load model (e.g. when you skip the preprocessor):
    model: models.BiVentricle = models.HeartModel.load_model(path_to_model)
    ## instantiate simulator.
    simulator = MechanicsSimulator(model, lsdyna_path, "smp", num_cpus=4)
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

    # ADVANCED example: create array of different settings and use these in simulator.
    settings = SimulationSettings()
    settings.load_defaults()
    from pint import Quantity
    import copy
    from typing import List

    # prepare list of settings.
    list_of_settings: List[SimulationSettings] = []
    k1_list = [
        Quantity(0.00336, "newton / millimeter"),
        Quantity(0.00436, "newton / millimeter"),
        Quantity(0.00536, "newton / millimeter"),
        Quantity(0.00636, "newton / millimeter"),
    ]
    for k1 in k1_list:
        s = copy.deepcopy(settings)
        s.mechanics.material.myocardium["isotropic"]["k1"] = k1
        list_of_settings.append(s)

    # Load model (e.g. when you skip the preprocessor):
    model: models.BiVentricle = models.HeartModel.load_model(path_to_model)

    # loop over all settings and run simulation.
    base_directory = os.path.join(model.info.workdir, "simulations")
    for ii, settings in enumerate(list_of_settings):
        ## create simulation directory.
        sim_directory = os.path.join(base_directory, "simulation-{:03d}".format(ii))
        ## instanatiate simulator.
        simulator = MechanicsSimulator(
            model, lsdyna_path, "smp", num_cpus=4, simulation_directory=sim_directory
        )
        ## load default settings.
        simulator.settings = settings
        simulator.settings.mechanics.analysis.end_time = Quantity(20, "ms")

        ## compute the fiber orientation
        # simulator.compute_fibers()

        ## compute the stress free configuration
        # simulator.compute_stress_free_configuration()

        ## do the main simulation
        simulator.simulate()
