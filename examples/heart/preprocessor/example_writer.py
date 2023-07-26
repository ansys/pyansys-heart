"""Example to pre-process data from Strocchi2020 and Cristobal2021."""
import os
import pathlib 

import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.support import run_preprocessor
import ansys.heart.writer.dynawriter as writers

PROJECT_DIRECTORY =  pathlib.Path(__file__).absolute().parents[3]
PATH_TO_CASE =  os.path.join(PROJECT_DIRECTORY, "downloads\\Strocchi2020\\01\\01.case")
WORKING_DIRECTORY = os.path.join(pathlib.Path(PATH_TO_CASE).parent, "BiVentricle")

if __name__ == "__main__":
    """BiVentricle example.
    model_type can be changed to BiVentricle or FullHeart or FourChamber model based on user
    requirements

    1. Extracts simulation mesh and creates "heart_model.pickle"
    2. Writes files for mechanics, zero-pressure, fiber generation, and purkinje using "heart_model.pickle" file

    Please change paths according to your workspace
    """
    path_to_model = os.path.join(WORKING_DIRECTORY, "heart_model.pickle")

    use_preprocessor = True
    write_lsdyna_files = True

    # Preprocessing geometry, remeshing with mesh_size=2
    if use_preprocessor:
        model = run_preprocessor(
            model_type=models.BiVentricle,
            database="Strocchi2020",
            path_original_mesh=PATH_TO_CASE,
            work_directory=WORKING_DIRECTORY,
            path_to_model=path_to_model,
            mesh_size=2.0,
        )

    # Load model (e.g. when you skip the preprocessor):
    model = models.HeartModel.load_model(path_to_model)
    if not isinstance(model, models.HeartModel):
        exit()
    model.info.workdir = os.path.join(WORKING_DIRECTORY, "simulation_files")

    # Write LS-DYNA k files for mechanics, zero-pressure, fiber generation, and purkinje
    # generation
    if write_lsdyna_files:
        for writer in (
            writers.ElectrophysiologyDynaWriter(model),
            writers.MechanicsDynaWriter(model),
            writers.ZeroPressureMechanicsDynaWriter(model),
            writers.FiberGenerationDynaWriter(model),
            writers.PurkinjeGenerationDynaWriter(model),
        ):
            exportdir = os.path.join(
                writer.model.info.workdir,
                writer.__class__.__name__.lower().replace("dynawriter", ""),
            )

            writer.update()
            writer.export(exportdir)

            writer.model.mesh.write_to_vtk(
                os.path.join(writer.model.info.workdir, "volume_mesh.vtk")
            )
    print("done")
