import os

from ansys.heart.postprocessor.binout_helper import Elout
from ansys.heart.postprocessor.bulleye import bullseye_plot
from ansys.heart.postprocessor.compute_strain import (
    compute_AHA17_segment_strain,
    compute_myocardial_strain,
)
from ansys.heart.preprocessor.models import HeartModel
from matplotlib import pyplot as plt
import numpy as np

if __name__ == "__main__":
    os.chdir("..\\data")
    model = HeartModel.load_model("heart_model.pickle")
    res = Elout("binout0000")
    strain = compute_myocardial_strain(model, res, refrence_time=0)
    np.save("strain", strain)
    strain = np.load("strain.npy")
    aha_strain = compute_AHA17_segment_strain(model, res, strain)

    # plot longitudinal strain for all segment
    plt.plot(res.time, aha_strain[0, :, :], label=f"segment")
    plt.title("longitudinal strain")
    plt.legend()
    plt.show()

    # plot bulleye
    fig, ax = plt.subplots(figsize=(24, 16), nrows=1, ncols=3, subplot_kw=dict(projection="polar"))
    fig.canvas.manager.set_window_title("Left Ventricle Bulls Eyes (AHA)")
    for i in range(3):
        bullseye_plot(ax[i], aha_strain[i, 3, :])  # at second time step

    ax[0].set_title("longitudinal")
    ax[1].set_title("radial")
    ax[2].set_title("circumferential")
    plt.show()

    # write strain into vtk
    import meshio

    ele_id = np.where(model.aha_ids != 0)[0]
    elems = model.mesh.tetrahedrons[ele_id]
    nodes = model.mesh.nodes[np.unique(elems.ravel())]
    _, a = np.unique(elems, return_inverse=True)
    connect = a.reshape(elems.shape)

    for i_time in range(len(res.time)):
        meshio.write_points_cells(
            f"strain_{i_time}.vtk",
            nodes,
            [("tetra", connect)],
            cell_data={
                "lrc_strain": [strain[:, i_time, :].T],
            },
        )
