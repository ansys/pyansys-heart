"""For information only, example not yet stable"""
from ansys.heart.postprocessor.EPpostprocessor import EPpostprocessor
from ansys.heart.preprocessor.models import HeartModel
import pyvista as pv

if __name__ == "__main__":
    model: HeartModel
    model = HeartModel.load_model(r"path_to_model.pickle")

    results_path = r"path_to_simulation_results"

    postproc = EPpostprocessor(results_path=results_path, model=model)
    postproc.read_EP_results()

    pl = pv.Plotter()
    pl.add_mesh(postproc.mesh, scalars="activation_time")
    pl.show()
