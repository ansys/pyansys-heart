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


def main(args):
    # input
    database = args.database
    cases = args.case
    types = [args.type]
    size = args.meshsize
    root_folder = args.root

    #############################################################
    # package import
    import ansys.heart.preprocessor.models as models

    #
    if database == "Strocchi2020":
        extension = "case"
    elif database == "Rodero2021":
        extension = "vtk"

    for id in cases:
        case_str = "{0:02d}".format(id)
        root = os.path.join(root_folder, database, case_str)
        case_file = os.path.join(root, case_str + "." + extension)

        print("Path to case file: %s" % case_file)

        for model_type in types:
            workdir = os.path.join(root, f"{model_type}_{size}")
            path_to_model = str(Path(workdir, "heart_model.pickle"))

            # TODO: this will not work with v0.2
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
            elif model_type == "4c":
                model = models.FourChamber(info)

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
        "--database",
        help="Database to use.",
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
    parser.add_argument(
        "--type",
        help="Heart model type: lv: left-ventricular model, bv: biventricular model, "
        + "fh: full heart model",
        choices=["lv", "bv", "4c", "fh"],
        default="fh",
    )

    parser.add_argument(
        "--meshsize", help="Uniform Mesh size used for remeshing ", type=float, default=2.0
    )

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(args)
