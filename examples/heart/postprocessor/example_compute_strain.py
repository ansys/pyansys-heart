from ansys.heart.postprocessor.bulleye import bullseye_plot
from ansys.heart.postprocessor.compute_strain import (
    compute_AHA17_segment_strain,
    compute_myocardial_strain,
)
from ansys.heart.postprocessor.dpf_d3plot import D3plotReader
from ansys.heart.preprocessor.models import HeartModel
from matplotlib import pyplot as plt

if __name__ == "__main__":
    model: HeartModel
    model = HeartModel.load_model(r"D:\Heart20\healthy20\health03_BV_2mm\heart_model_bv.pickle")
    model.compute_left_ventricle_anatomy_axis(first_cut_short_axis=0.2)
    model.compute_left_ventricle_aha17()

    d3plot_file = r"D:\Heart20\healthy20\health03_BV_2mm\simulation\main-mechanics\d3plot"
    data = D3plotReader(d3plot_file)
    df = data.get_history_variable(hv_index=list(range(9)), at_frame=1)
    strain = compute_myocardial_strain(model, df.T, reference=None)

    # write strain into vtk
    import meshio
    import numpy as np

    ele_id = np.where(~np.isnan(model.aha_ids))[0]
    elems = model.mesh.tetrahedrons[ele_id]
    nodes = model.mesh.nodes[np.unique(elems.ravel())]
    _, a = np.unique(elems, return_inverse=True)
    connect = a.reshape(elems.shape)
    meshio.write_points_cells(
        f"strain_{0}.vtk",
        nodes,
        [("tetra", connect)],
        cell_data={
            "lrc_strain": [strain],
            "aha17": [model.aha_ids[model.aha_ids > 0]],
        },
    )

    aha_strain = compute_AHA17_segment_strain(model, strain)

    # plot bulleye
    fig, ax = plt.subplots(figsize=(24, 16), nrows=1, ncols=3, subplot_kw=dict(projection="polar"))
    fig.canvas.manager.set_window_title("Left Ventricle Bulls Eyes (AHA)")
    for i in range(3):
        bullseye_plot(ax[i], aha_strain[:, i])  # at second time step

    ax[0].set_title("longitudinal")
    ax[1].set_title("radial")
    ax[2].set_title("circumferential")
    plt.show()
