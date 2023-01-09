from ansys.heart.postprocessor.D3plotPost import LVContourExporter
from ansys.heart.preprocessor.models import HeartModel
import matplotlib.pyplot as plt

if __name__ == "__main__":
    """
    Show exporting LV motion in vtk
    """
    model: HeartModel
    model = HeartModel.load_model(r"D:\Heart20\last_model\heart_model.pickle")
    model.compute_left_ventricle_anatomy_axis(first_cut_short_axis=0.2)

    d3plot_file = "D:/Heart20/last_model/main/d3plot"

    exporter = LVContourExporter(d3plot_file, model)

    exporter.export_contour_to_vtk("l4cv", model.l4cv_axis)

    exporter.export_contour_to_vtk("l2cv", model.l2cv_axis)

    normal = model.short_axis["normal"]
    p_start = model.short_axis["center"]
    for ap in model.left_ventricle.apex_points:
        if ap.name == "apex epicardium":
            p_end = ap.xyz

    for icut in range(2):
        p_cut = p_start + (p_end - p_start) * icut / 2
        cutter = {"center": p_cut, "normal": normal}
        exporter.export_contour_to_vtk(f"shor_{icut}", cutter)

    # x,y= b.compute_lvls()
    exporter.export_lvls_to_vtk("lvls")
    exporter.plot_lvls()
    plt.show()
