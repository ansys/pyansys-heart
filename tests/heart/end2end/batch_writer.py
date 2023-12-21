import argparse
import os
from pathlib import Path

#############################################################
os.environ["USE_OLD_HEART_MODELS"] = "1"
root_folder = r"D:\ansysdev"


# right atrium appendage apex is manually selected
right_appendage_apex = {
    "Cristobal2021": {
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
    import ansys.heart.preprocessor.models as models
    import ansys.heart.writer.dynawriter as writers

    # input
    database = args.database
    cases = args.case
    folder = args.folder
    writer_type = args.writer

    #
    model_name = "heart_model.pickle"
    simu_folder = "writer_" + writer_type

    #
    for id in cases:
        case_str = "{0:02d}".format(id)
        root = os.path.join(root_folder, database, case_str)

        workdir = os.path.join(root, f"{folder}")
        path_to_model = str(Path(workdir, model_name))

        model = models.HeartModel.load_model(path_to_model)

        #
        if writer_type == "meca":
            w = writers.MechanicsDynaWriter(model)

        elif writer_type == "zerop":
            w = writers.ZeroPressureMechanicsDynaWriter(model)

        elif writer_type == "epmeca":
            w = writers.ElectroMechanicsDynaWriter(model)

        elif writer_type == "fiber":
            w = writers.FiberGenerationDynaWriter(model)

        elif writer_type == "uvc":
            w = writers.UHCWriter(model, type=writer_type)

        elif writer_type == "la_fiber":
            w = writers.UHCWriter(model, type=writer_type)
        elif writer_type == "ra_fiber":
            try:
                coord = right_appendage_apex[database][str(id)]
            except KeyError:
                print(f"raa coordinate need to be defined for case {id} of {database}")
                exit()
            w = writers.UHCWriter(model, type=writer_type, raa=coord)

        else:
            print(f"Unknown writer type: {writer_type}")
            exit()

        w.update()
        w.export(os.path.join(workdir, simu_folder))


if __name__ == "__main__":
    # Create an argument parser
    parser = argparse.ArgumentParser(description="EndToEnd Test: Batch run writer")

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
    parser.add_argument(
        "--folder",
        help="Folder of heart model type",
        type=str,
        # required=True,
        default="fh_1.5",
    )

    parser.add_argument(
        "--writer",
        help="Writer type",
        choices=["zerop", "fiber", "la_fiber", "ra_fiber", "meca", "ep", "epmeca", "uvc"],
        # required=True,
        default="meca",
    )
    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(args)
