"""Example to pre-process data from Strocchi2020 and Cristobal2021."""
import glob as glob
import os
import pathlib

import ansys.heart.preprocessor.models as models

DOWNLOAD_DIR = os.path.join(pathlib.Path(__file__).parents[3], "downloads")

if __name__ == "__main__":
    """Loop over all available cases and databases."""

    databases_to_run = ["Cristobal2021"]
    # databases_to_run = ["Strocchi2020"]

    for database in databases_to_run:
        database_dir = os.path.join(DOWNLOAD_DIR, database)
        case_folders = glob.glob(os.path.join(database_dir, "*"))

        for case_folder in case_folders:

            if database == "Strocchi2020":
                case_path = os.path.join(case_folder, pathlib.Path(case_folder).name + ".case")
            elif database == "Cristobal2021":
                case_path = os.path.join(case_folder, pathlib.Path(case_folder).name + ".vtk")

            case_num = os.path.split(pathlib.Path(case_path).parent)[-1]
            print("Database: {} | case: {}".format(database, case_num))

            workdir = os.path.join(pathlib.Path(case_path).parent, "ExtractEndocardium")

            info = models.ModelInfo(
                database=database, work_directory=workdir, path_to_case=case_path
            )
            info.create_workdir()

            model = models.FourChamber(info)
            model.read_input_mesh()
            model._remove_unused_tags()
            model._prepare_for_meshing()

            save_directory = os.path.join(
                DOWNLOAD_DIR, "extracted_atrial_endocardia", database, case_num
            )
            # for interface in model.mesh_raw.interfaces:
            #     print("interface: %s " % interface.name )

            # for boundary in model.mesh_raw.boundaries:
            #     print("boundary: %s " % boundary.name )
            if not os.path.isdir(save_directory):
                os.makedirs(save_directory)

            for boundary in model.mesh_raw.boundaries:
                if "atrium" in boundary.name:
                    boundary.write_to_stl(os.path.join(save_directory, boundary.name + ".stl"))
            for interface in model.mesh_raw.interfaces:
                if "atrium" in interface.name:
                    interface.write_to_stl(os.path.join(save_directory, interface.name + ".stl"))

    pass
