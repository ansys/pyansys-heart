"""unit test for post-processing."""
import os
import pathlib as Path

from ansys.heart.postprocessor.SystemModelPost import SystemModelPost
from ansys.heart.postprocessor.aha17_strain import AhaStrainCalculator
from ansys.heart.postprocessor.auto_process import mech_post, zerop_post
from ansys.heart.postprocessor.exporter import LVContourExporter
from ansys.heart.preprocessor.models import HeartModel
import pytest

model: HeartModel
test_dir: str


pytestmark = pytest.mark.local


@pytest.fixture(autouse=True, scope="module")
def get_data():
    global test_dir, model

    # TODO: test case in locally saved, need to upload to Github

    test_dir = r"D:\PyAnsys-Heart\test_case\test_lv"
    model = HeartModel.load_model(Path.Path(test_dir) / "model_with_fiber.pickle")
    model.compute_left_ventricle_anatomy_axis()
    model.compute_left_ventricle_aha17()


@pytest.mark.xfail(reason="Test requires local data.")
def test_compute_myocardial_strain():
    d3plot = Path.Path(test_dir) / "main-mechanics" / "d3plot"

    s = AhaStrainCalculator(model, d3plot)
    element_lrc, aha_lrc, element_lrc_averaged = s._compute_myocardial_strain(1)
    assert aha_lrc[-1, -1] == pytest.approx(0.08878163)


@pytest.mark.xfail(reason="Test requires local data.")
def test_compute_aha_strain():
    d3plot = Path.Path(test_dir) / "main-mechanics" / "d3plot"

    s = AhaStrainCalculator(model, d3plot)
    aha_lrc = s.compute_aha_strain(".")
    assert aha_lrc[1, -1] == pytest.approx(0.08878163)


@pytest.mark.xfail(reason="Test requires local data.")
def test_mech_post():
    dct = mech_post(Path.Path(test_dir) / "main-mechanics", model)
    assert os.path.exists(Path.Path(test_dir) / "main-mechanics" / "post")


@pytest.mark.xfail(reason="Test requires local data.")
def test_zerop_post():
    dct = zerop_post(Path.Path(test_dir) / "zeropressure", model)
    assert dct["True left ventricle volume (mm3)"] == pytest.approx(288876.8)
    assert os.path.exists(Path.Path(test_dir) / "zeropressure" / "post")


@pytest.mark.xfail(reason="Test requires local data.")
def test_contour_exporter():
    d3plot = Path.Path(test_dir) / "main-mechanics" / "d3plot"
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
    # assert np.allclose(
    #     np.array(contours[0].GetPoint(0)),
    #     np.array([7.273597903198675, 128.824467818063, 387.9997416596942]),
    # )
    # assert np.allclose(
    #     np.array(contours[1].GetPoint(0)),
    #     np.array([7.254521621126342, 128.83922762291593, 388.01083324324804]),
    #     rtol=0.001,
    # )


@pytest.mark.xfail(reason="Test requires local data.")
def test_lvls():
    d3plot = Path.Path(test_dir) / "main-mechanics" / "d3plot"
    exporter = LVContourExporter(d3plot, model)
    p1, p2 = exporter._compute_lvls()
    exporter.export_lvls_to_vtk(folder="lvls")
    print()
    # assert np.allclose(
    #     p1,
    #     np.array(
    #         [[14.84525339, 138.85722581, 381.68174194], [14.88606218, 138.8574518, 381.66178657]]
    #     ),
    #     rtol=0.001,
    # )
    # assert np.allclose(
    #     p2,
    #     np.array([[70.7569, 72.7259, 351.93], [70.79526416, 72.73520587, 351.97354634]]),
    #     rtol=0.001,
    # )


class TestSystemModelPost:
    @pytest.fixture
    def system_model(self):
        return SystemModelPost(Path.Path(test_dir) / "main-mechanics")

    @pytest.mark.xfail(reason="Test requires local data.")
    def test_plot_pv_loop(self, system_model):
        ef = system_model.get_ejection_fraction(t_start=2, t_end=3)
        fig = system_model.plot_pv_loop(ef=ef)
        fig.savefig("pv_{0:d}.png".format(0))
        assert os.path.isfile("pv_0.png")
        # ffmpeg -f image2 -i pv_%d.png output.mp4
