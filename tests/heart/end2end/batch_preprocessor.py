import argparse
import os
from pathlib import Path

import ansys.heart.preprocessor.models as models

#############################################################
os.environ["USE_OLD_HEART_MODELS"] = "1"


def main(args):
    # input
    database = args.database
    cases = args.case
    types = [args.type]
    size = args.meshsize
    root_folder = args.root

    #
    if database == "Strocchi2020":
        extension = "case"
    elif database == "Cristobal2021":
        extension = "vtk"

    for id in cases:
        case_str = "{0:02d}".format(id)
        root = os.path.join(root_folder, database, case_str)
        case_file = os.path.join(root, case_str + "." + extension)

        print("Path to case file: %s" % case_file)

        for model_type in types:
            workdir = os.path.join(root, f"{model_type}_{size}")
            path_to_model = str(Path(workdir, "heart_model.pickle"))

            info = models.ModelInfo(
                database=database,
                path_to_case=case_file,
                work_directory=workdir,
                path_to_model=path_to_model,
                add_blood_pool=False,
                mesh_size=size,
            )

            # create the working directory
            info.create_workdir()
            # clean the working directory
            info.clean_workdir(extensions_to_remove=[".stl", ".vtk", ".msh.h5"])
            # dump information to stdout
            info.dump_info()
            if model_type == "lv":
                model = models.LeftVentricle(info)
            elif model_type == "bv":
                model = models.BiVentricle(info)
            elif model_type == "fh":
                model = models.FullHeart(info)

            # extract the simulation mesh
            model.extract_simulation_mesh()

            # dump the model to disk for future use
            model.dump_model(path_to_model)
            # print the resulting information
            model.print_info()


if __name__ == "__main__":
    # Create an argument parser
    parser = argparse.ArgumentParser(description="EndToEnd Test: Batch run preprocessor")

    # Define command-line arguments
    parser.add_argument(
        "--root",
        help="Root folder. The script will look for cases relative to this part.",
        default="D:\\ansysdev",
    )

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

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(args)
