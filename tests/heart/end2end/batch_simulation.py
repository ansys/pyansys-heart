# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import argparse
import os
from pathlib import Path

# right atrium appendage apex is manually selected
right_appendage_apex = {
    "Rodero2021": {
        "1": [39, 29, 98],
        "2": [40, 38, 101],
        "3": [40, 57, 98],
        "4": [35, 34, 101],
        "5": [30, 32, 86],
        # ...
    },
    "Strocchi2020": {
        "1": [-33, 82, 417],
        "2": [-35, 84, -169],
        "3": [-63, 76, 184],
        "4": [-27, 107, 279],
        "5": [-23, 67, 231],
        # ...
    },
}


def main(args):
    # input
    database = args.database
    cases = args.case
    folder = args.folder
    simulation = args.simulation
    root_folder = args.root

    #############################################################
    # package import
    import ansys.heart.preprocessor.models.v0_1.models as models
    from ansys.heart.simulator.simulator import (
        DynaSettings,
        EPMechanicsSimulator,
        EPSimulator,
        MechanicsSimulator,
    )

    #############################################################
    ##  instantiate dyna settings object
    dyna_settings = DynaSettings(
        lsdyna_path=args.lsdynapath,
        dynatype="intelmpi",
        num_cpus=6,
        platform="wsl",
    )
    # dyna_settings._set_env_variables()
    #############################################################

    #
    simu_folder = "end2end"

    #
    for id in cases:
        case_str = "{0:02d}".format(id)
        root = os.path.join(root_folder, database, case_str)

        workdir = os.path.join(root, f"{folder}")
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
            simulator.compute_fibers()

            if isinstance(model, models.FourChamber):
                try:
                    coord = right_appendage_apex[database][str(id)]
                except KeyError:
                    print(f"raa coordinate need to be defined for case {id} of {database}")
                    exit()

                simulator.compute_left_atrial_fiber()
                simulator.compute_right_atrial_fiber(appendage=coord)

                simulator.model.left_atrium.has_fiber = True
                simulator.model.left_atrium.is_active = True
                simulator.model.right_atrium.has_fiber = True
                simulator.model.right_atrium.is_active = True

            _ = simulator.create_stiff_ventricle_base()
            _ = simulator.model._create_atrial_stiff_ring()

            simulator.compute_stress_free_configuration()

            simulator.compute_purkinje()
            if isinstance(model, models.FourChamber):
                simulator.compute_conduction_system()

            simulator.simulate()


if __name__ == "__main__":
    # Create an argument parser
    parser = argparse.ArgumentParser(description="EndToEnd Test: Batch run simulation")

    # Define command-line arguments
    parser.add_argument(
        "--heartversion",
        help="Heart model version. 0: Uses HeartModels from old version of models.py,"
        + "1: Uses HeartModels from new version of models.py (models_new.py)",
        type=str,
        default="0",
    )

    parser.add_argument(
        "--root",
        help="Root folder. The script will look for cases relative to this folder.",
        default="D:\\ansysdev",
    )

    parser.add_argument(
        "--lsdynapath",
        help="Path to LS-DYNA",
        default=r"D:\ansysdev\lsdyna_mpp\mppdyna_d_sse2_linux86_64_intelmmpi_105630",
    )

    parser.add_argument(
        "--database",
        help="Rodero2021 or Strocchi2020",
        choices=["Rodero2021", "Strocchi2020"],
        default="Rodero2021",
    )

    # Use nargs='+' to accept a list of float values for option2
    parser.add_argument(
        "--case",
        help="e.g. 1 for case 1, 1,2 for case 1 and 2",
        type=lambda x: [int(val) for val in x.split(",")],
        default=[1],
    )
    parser.add_argument("--folder", help="Folder of heart model type", type=str, required=True)

    parser.add_argument(
        "--simulation", help="simulaton type", choices=["meca", "ep", "epmeca"], default="ep"
    )
    # Parse the command-line arguments
    args = parser.parse_args()

    # set right environment variable
    if args.heartversion == "0":
        os.environ["USE_OLD_HEART_MODELS"] = "1"
    else:
        os.environ["USE_OLD_HEART_MODELS"] = "0"

    # Call the main function with parsed arguments
    main(args)
