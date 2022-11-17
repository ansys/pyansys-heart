import os.path
import pathlib
import shutil

from ansys.heart.calibration.passive_calibration import create_calibration_folder
import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.support import run_preprocessor
import ansys.heart.writer.dynawriter as writers

if __name__ == "__main__":

    os.environ["REMOTING_SERVER_ADDRESS"] = "localhost"
    use_preprocessor = False
    write_lsdyna_files = False

    for i in range(2, 3):
        id = "{0:02d}".format(i)
        path_to_case = r"D:\Heart20\healthy20\{0:s}.vtk".format(id)
        workdir = os.path.join(pathlib.Path(path_to_case).parent, "health{0:s}_BV_2mm".format(id))

        path_to_model = os.path.join(workdir, "heart_model.pickle")

        if use_preprocessor:
            model = run_preprocessor(
                model_type=models.BiVentricle,
                database="Strocchi2020",
                path_original_mesh=path_to_case,
                work_directory=workdir,
                path_to_model=path_to_model,
                mesh_size=2,
            )

        if write_lsdyna_files:
            # write LS-DYNA files
            # Load model (e.g. when you skip the preprocessor):
            model = models.HeartModel.load_model(path_to_model)

            for writer in (
                writers.FiberGenerationDynaWriter(model),
                writers.ZeroPressureMechanicsDynaWriter(model),
            ):
                exportdir = os.path.join(
                    workdir,
                    writer.__class__.__name__.lower().replace("dynawriter", "_calibration"),
                )
                # writer.model.mesh.write_to_vtk(os.path.join(workdir, "volume_mesh.vtk"))
                writer.update()
                writer.export(exportdir)

        # Run fiber rule
        from ansys.heart.general import run_lsdyna

        os.chdir(os.path.join(workdir, "fibergeneration_calibration"))
        run_lsdyna("main.k")

        # Copy elements.k
        shutil.copy2(
            os.path.join(workdir, "fibergeneration_calibration", "element_solid_ortho.k"),
            os.path.join(workdir, "zeropressuremechanics_calibration", "solid_elements.k"),
        )
        # create calibration directory
        create_calibration_folder(os.path.join(workdir, "zeropressuremechanics_calibration"))
        # run lsopt.exe PassiveCalibration.lsopt
