"""Example to pre-process data from Strocchi2020 and Cristobal2021"""
import os
import pathlib

import ansys.heart.preprocessor.models as models
from ansys.heart.workflow.support import run_preprocessor
import ansys.heart.writer.dynawriter as writers

if __name__ == "__main__":

    """Full Heart example

    1. Extracts simulation mesh
    2. Writes files for mechanics, zero-pressure, fiber generation, and purkinje

    Please change paths
    """

    path_to_case = os.path.join(
        pathlib.Path(__file__).parents[3], "downloads\\Strocchi2020\\01\\01.case"
    )
    workdir = os.path.join(pathlib.Path(path_to_case).parent, "BiVentricle")
    path_to_model = os.path.join(workdir, "heart_model.pickle")

    use_preprocessor = True
    write_lsdyna_files = True

    if use_preprocessor:
        model = run_preprocessor(
            # model_type=models.FullHeart,
            model_type=models.BiVentricle,
            database="Strocchi2020",
            path_original_mesh=path_to_case,
            work_directory=workdir,
            path_to_model=path_to_model,
            mesh_size=1.5,
            add_blood_pool=False,
        )

    # write LS-DYNA files
    # Load model (e.g. when you skip the preprocessor):
    model = models.HeartModel.load_model(path_to_model)

    if write_lsdyna_files:
        for writer in (
            writers.ElectrophysiologyDynaWriter(model),
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
