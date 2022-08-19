"""Example to pre-process data from Strocchi2020 and Cristobal2021"""
import os
import pathlib
import glob as glob


import ansys.heart.preprocessor.models as models
from ansys.heart.workflow.support import run_preprocessor
import ansys.heart.writer.dynawriter as writers

DOWNLOAD_DIR = os.path.join(pathlib.Path(__file__).parents[3], "downloads")

if __name__ == "__main__":

    """
    Loops over all available cases, databases, and available models
    """

    use_preprocessor = True
    write_lsdyna_files = False

    databases_to_run = ["Strocchi2020", "Cristobal2021"]
    models_to_run = (models.LeftVentricle, models.BiVentricle, models.FourChamber, models.FullHeart)

    for database in databases_to_run:
        database_dir = os.path.join(DOWNLOAD_DIR, database)
        case_folders = glob.glob(os.path.join(database_dir, "*"))

        for case_folder in case_folders:

            if database == "Strocchi2020":
                case_path = os.path.join(case_folder, pathlib.Path(case_folder).name + ".case")
            elif database == "Cristobal2021":
                case_path = os.path.join(case_folder, pathlib.Path(case_folder).name + ".vtk")

            for model_type in models_to_run:
                workdir = os.path.join(
                    pathlib.Path(case_path).parent, model_type.__name__ + "Refactored"
                )
                path_to_model = os.path.join(workdir, "heart_model.pickle")

                if os.path.isfile(path_to_model):
                    print("Skipping: {0}  -  {1}".format(case_path, model_type.__name__))
                    continue

                try:
                    if use_preprocessor:
                        model = run_preprocessor(
                            model_type=model_type,
                            database="Strocchi2020",
                            path_original_mesh=case_path,
                            work_directory=workdir,
                            path_to_model=path_to_model,
                            mesh_size=2.0,
                        )
                        model.mesh.write_to_vtk(os.path.join(workdir, "volume_mesh.vtk"))
                except:
                    print(">> Failed to run case: {0} | {1}".format(case_path, model_type.__name__))
                    continue

                # write LS-DYNA files
                # Load model (e.g. when you skip the preprocessor):
                model = models.HeartModel.load_model(path_to_model)

                if write_lsdyna_files:
                    for writer in (
                        writers.MechanicsDynaWriter(model, "ConstantPreloadWindkesselAfterload"),
                        writers.ZeroPressureMechanicsDynaWriter(model),
                        writers.FiberGenerationDynaWriter(model),
                        writers.PurkinjeGenerationDynaWriter(model),
                    ):
                        exportdir = os.path.join(
                            writer.model.info.workdir,
                            writer.__class__.__name__.lower().replace("dynawriter", ""),
                        )

                        writer.model.mesh.write_to_vtk(
                            os.path.join(writer.model.info.workdir, "volume_mesh.vtk")
                        )
                        writer.update()
                        writer.export(exportdir)

    pass
