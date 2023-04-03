"""For information only, example not yet stable"""
import pathlib as Path

from ansys.heart.postprocessor.EPpostprocessor import EPpostprocessor
from ansys.heart.preprocessor.models import HeartModel
import pyvista as pv

if __name__ == "__main__":
    model: HeartModel
    model = HeartModel.load_model(
        r"D:\REPOS\pyheart-lib\downloads\Strocchi2020\01\BiVentricle\simulation-EP\model_with_fiber.pickle"
    )

    results_path = (
        r"D:\REPOS\pyheart-lib\downloads\Strocchi2020\01\BiVentricle\simulation-EP\main-ep"
    )

    postproc = EPpostprocessor(results_path=results_path, model=model)
    postproc.read_EP_results()

    pl = pv.Plotter()
    pl.add_mesh(postproc.mesh, scalars="activation_time")
    pl.show()
