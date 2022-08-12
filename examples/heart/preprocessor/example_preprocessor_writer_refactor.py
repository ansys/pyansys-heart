"""Example to pre-process data from Strocchi2020 and Cristobal2021"""
import os
import pathlib


import ansys.heart.preprocessor.models as models
import ansys.heart.writer.dynawriter as writers


if __name__ == "__main__":

    """Full Heart example

    1. Extracts simulation mesh
    2. Writes files for mechanics, zero-pressure, and fiber generation
    """

    path_to_case = "D:\\development\\pyheart-lib\\pyheart-lib\\downloads\\Strocchi2020\\02\\02.case"
    workdir = os.path.join(pathlib.Path(path_to_case).parent, "FullHeartRefactored")
    info = models.ModelInfo(
        database="Strocchi2020", work_directory=workdir, path_to_case=path_to_case
    )
    info.mesh_size = 1.5

    path_to_model = os.path.join(info.workdir, "heart_model.pickle")

    model = models.FullHeart(info)
    run_preprocessor = True

    if run_preprocessor:
        info.clean_workdir(remove_all=True)
        info.create_workdir()
        info.dump_info()
        model.extract_simulation_mesh()
        model.dump_model(path_to_model)
        model.print_info()
        model.info.clean_workdir([".stl", ".vtk", ".jou", ".log", ".trn"])

    # write LS-DYNA files
    # Load model (e.g. when you skip the preprocessor):
    model = models.HeartModel.load_model(path_to_model)

    for writer in (
        writers.MechanicsDynaWriter(model, "ConstantPreloadWindkesselAfterload"),
        writers.ZeroPressureMechanicsDynaWriter(model),
        writers.FiberGenerationDynaWriter(model),
    ):
        exportdir = os.path.join(
            writer.model.info.workdir, writer.__class__.__name__.lower().replace("dynawriter", "")
        )

        writer.update()
        writer.export(exportdir)

    pass
