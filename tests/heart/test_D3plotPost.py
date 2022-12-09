"""unit test for D3plotPost."""
import os
import pathlib

from ansys.heart.postprocessor.D3plotPost import D3plotExporter, LVContourExporter
from ansys.heart.preprocessor.models import HeartModel
import numpy as np
import pytest

from .conftest import get_assets_folder

#
# todo: only work with Paraview

d3plot: pathlib.Path
model: HeartModel


@pytest.fixture(autouse=True, scope="module")
def get_data():
    global model, d3plot
    assets_folder = get_assets_folder()
    path_to_reference_model = os.path.join(
        assets_folder,
        "reference_models",
        "strocchi2020/01",
        "BiVentricle",
    )
    model = HeartModel.load_model(os.path.join(path_to_reference_model, "heart_model.pickle"))
    model.compute_left_ventricle_anatomy_axis()
    d3plot = os.path.join(path_to_reference_model, "post", "d3plot")


def test_convert_to_vtk():
    exporter = D3plotExporter(d3plot, model)
    exporter.convert_to_vtk(model.part_names, "temp")
    target = os.path.join(exporter.work_dir, "temp")
    assert os.listdir(target) == ["model_0.vtk", "model_1.vtk"]


def test_contour_exporter():
    exporter = LVContourExporter(d3plot, model)
    contours = exporter.export_contour_to_vtk("l4cv", model.l4cv_axis)
    # exporter.export_contour_to_vtk('l2cv', model.l2cv_axis)
    #
    # normal = model.short_axis["normal"]
    #
    # p_start = model.short_axis["center"]
    # for ap in model.left_ventricle.apex_points:
    #     if ap.name == "apex epicardium":
    #         p_end = ap.xyz
    # for icut in range(2):
    #     p_cut = p_start + (p_end - p_start) * icut / 2
    #     cutter = {"center": p_cut, "normal": normal}
    #     exporter.export_contour_to_vtk(f"shor_{icut}",cutter)
    print()
    print(contours[0].GetPoint(0))
    print(contours[1].GetPoint(0))
    assert np.allclose(
        np.array(contours[0].GetPoint(0)),
        np.array([7.273597903198675, 128.824467818063, 387.9997416596942]),
    )
    assert np.allclose(
        np.array(contours[1].GetPoint(0)),
        np.array([7.254521621126342, 128.83922762291593, 388.01083324324804]),
        rtol=0.001,
    )


def test_lvls():
    exporter = LVContourExporter(d3plot, model)
    p1, p2 = exporter.compute_lvls()
    print()
    assert np.allclose(
        p1,
        np.array(
            [[14.84525339, 138.85722581, 381.68174194], [14.88606218, 138.8574518, 381.66178657]]
        ),
        rtol=0.001,
    )
    assert np.allclose(
        p2,
        np.array([[70.7569, 72.7259, 351.93], [70.79526416, 72.73520587, 351.97354634]]),
        rtol=0.001,
    )
