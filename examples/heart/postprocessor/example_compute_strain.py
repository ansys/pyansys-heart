"""Example show how to compute myocardium strain"""
from ansys.heart.preprocessor.models import HeartModel
from matplotlib import pyplot as plt

from heart.postprocessor.aha17_strain import AhaStrainCalculator

if __name__ == "__main__":
    model: HeartModel
    model = HeartModel.load_model(r".")
    model.compute_left_ventricle_anatomy_axis(first_cut_short_axis=0.2)
    model.compute_left_ventricle_aha17()

    d3plot_file = r"."

    aha_calculator = AhaStrainCalculator(model, d3plot_file)

    # get LRC strain at time=0 and export to vtk
    aha_strain0 = aha_calculator.compute_aha_strain_once(frame=0, save_vtk=True)

    # bulleye plot
    fig, ax = plt.subplots(figsize=(24, 16), nrows=1, ncols=3, subplot_kw=dict(projection="polar"))
    fig.canvas.manager.set_window_title("Left Ventricle Bulls Eyes (AHA)")
    for i in range(3):
        aha_calculator.bullseye_plot(ax[i], aha_strain0[:, i])
    ax[0].set_title("longitudinal")
    ax[1].set_title("radial")
    ax[2].set_title("circumferential")
    plt.show()

    # export all strain into csv (Long)
    # aha_strain = aha_calculator.compute_aha_strain(out_file='strain.csv')
