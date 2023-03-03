"""Example to pre-process (morphed) data from Casis.

Note
----
Uses 'protected' methods from the HeartModel class to chain
various operations. Moreover, some non-generic steps are required such
as compute uvc, providing a list of surfaces instead of an entire model, etc

"""
import os

import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.support import run_preprocessor
import ansys.heart.writer.dynawriter as writers

if __name__ == "__main__":
    """BiVentricle LabeledSurface data."""

    os.environ["PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION"] = "python"
    os.environ["REMOTING_SERVER_ADDRESS"] = "localhost"

    preproc = True
    if preproc:
        model = run_preprocessor(
            model_type=models.LeftVentricle,
            database="LabeledSurface",
            path_original_mesh=os.path.join("data", "lv_surface.vtk"),
            work_directory=os.path.join("data", "LeftVentricle"),
            mesh_size=1.5,
        )
    else:
        # Load model (e.g. when you skip the preprocessor):
        path_to_model = os.path.join("data", "BiVentricle", "heart_model.pickle")
        model = models.HeartModel.load_model(path_to_model)

    # write LS-DYNA files
    write_lsdyna_files = True
    if write_lsdyna_files:
        for writer in (
            writers.FiberGenerationDynaWriter(model),
            writers.ElectrophysiologyDynaWriter(model),
            writers.MechanicsDynaWriter(model),
            writers.ZeroPressureMechanicsDynaWriter(model),
            writers.PurkinjeGenerationDynaWriter(model),
        ):
            exportdir = os.path.join(
                writer.model.info.workdir,
                "new",
                writer.__class__.__name__.lower().replace("dynawriter", ""),
            )

            writer.model.mesh.write_to_vtk(
                os.path.join(writer.model.info.workdir, "volume_mesh.vtk")
            )
            writer.update()
            writer.export(exportdir)

    pass
