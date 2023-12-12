import argparse
import os
from pathlib import Path

import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.simulator import (
    DynaSettings,
    EPMechanicsSimulator,
    EPSimulator,
    MechanicsSimulator,
)

#############################################################
os.environ["USE_OLD_HEART_MODELS"] = "1"
root_folder = r"D:\ansysdev"

#############################################################
##  instantiate dyna settings object
dyna_settings = DynaSettings(
    lsdyna_path=r"D:\ansysdev\lsdyna_mpp\mppdyna_d_sse2_linux86_64_intelmmpi_105630",
    dynatype="intelmpi",
    num_cpus=6,
    platform="wsl",
)
# dyna_settings._set_env_variables()
#############################################################


def main(args):
    # input
    database = args.database
    cases = args.case
    types = [args.type]
    size = args.meshsize
    simulation = args.simulation

    #
    simu_folder = "end2end"

    #
    for id in cases:
        case_str = "{0:02d}".format(id)
        root = os.path.join(root_folder, database, case_str)

        for model_type in types:
            workdir = os.path.join(root, f"{model_type}_{size}")
            path_to_model = str(Path(workdir, "heart_model.pickle"))

            model = models.HeartModel.load_model(path_to_model)

            # instantiate simulator object
            if simulation == "meca":
                simulator = MechanicsSimulator(
                    model=model,
                    dyna_settings=dyna_settings,
                    simulation_directory=os.path.join(workdir, simu_folder),
                )
                simulator.settings.load_defaults()

                simulator.compute_fibers()
                simulator.compute_stress_free_configuration()
                simulator.simulate()

            elif simulation == "ep":
                simulator = EPSimulator(
                    model=model,
                    dyna_settings=dyna_settings,
                    simulation_directory=os.path.join(workdir, simu_folder),
                )
                simulator.settings.load_defaults()
                #
                simulator.compute_fibers()
                simulator.compute_purkinje()
                if isinstance(model, models.FourChamber):
                    simulator.compute_conduction_system()
                simulator.simulate()

            elif simulation == "epmeca":
                simulator = EPMechanicsSimulator(
                    model=model,
                    dyna_settings=dyna_settings,
                    simulation_directory=os.path.join(workdir, simu_folder),
                )
                simulator.settings.load_defaults()
                #
                if isinstance(model, models.FourChamber):
                    raise NotImplementedError("We need right appendage apex point...")
                    simulator.compute_left_atrial_fiber()
                    simulator.compute_right_atrial_fiber(appendage=[39, 29, 98])

                    simulator.model.left_atrium.has_fiber = True
                    simulator.model.left_atrium.is_active = True
                    simulator.model.right_atrium.has_fiber = True
                    simulator.model.right_atrium.is_active = True

                simulator.compute_stress_free_configuration()

                simulator.compute_purkinje(after_zerop=True)
                if isinstance(model, models.FourChamber):
                    simulator.compute_conduction_system(after_zerop=True)

                simulator.simulate()


if __name__ == "__main__":
    # Create an argument parser
    parser = argparse.ArgumentParser(description="EndToEnd Test: Batch run simulation")

    # Define command-line arguments
    parser.add_argument(
        "--database",
        help="Cristobal2021 or Strocchi2020",
        choices=["Cristobal2021", "Strocchi2020"],
        default="Cristobal2021",
    )

    # Use nargs='+' to accept a list of float values for option2
    parser.add_argument(
        "--case",
        help="e.g. 1 for case 1, 1,2 for case 1 and 2",
        type=lambda x: [int(val) for val in x.split(",")],
        default=[1],
    )
    parser.add_argument("--type", help="Heart model type", choices=["lv", "bv", "fh"], default="fh")

    parser.add_argument("--meshsize", help="Mesh Size", type=float, default=2.0)

    parser.add_argument(
        "--simulation", help="simulaton type", choices=["meca", "ep", "epmeca"], default="ep"
    )
    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(args)
