"""Compute myocardial strain."""
from ansys.heart.postprocessor.binout_helper import Elout
from ansys.heart.postprocessor.bulleye import bullseye_plot
from ansys.heart.preprocessor.models import HeartModel
import matplotlib.pyplot as plt
import numpy as np


def compute_myocardial_strain(model: HeartModel, elout_res: Elout, refrence_time=0):
    """
    Compute left ventricle myocardial strain.

    Parameters
    ----------
    model: HeartModel
    elout_res: Elout object
    refrence_time: reference time of strain

    Returns
    -------
    strain 3 * time * element
    """
    # model info
    model.compute_left_ventricle_axis()
    model.compute_left_ventricle_AHA17()
    e_l, e_r, e_c = model.compute_left_ventricle_element_cs()

    # check
    ele_id = np.where(model.aha_ids != 0)[0] + 1
    if not np.array_equal(np.sort(ele_id), np.sort(elout_res.ids)):
        Exception("model id cannot match with elout file.")

    indices = np.where(np.in1d(elout_res.ids, ele_id))[0]
    strain = np.zeros((3, len(elout_res.time), len(ele_id)))
    def_grad = elout_res.get_history_variable()[:, indices, 0:9]

    # get strain at reference time
    if refrence_time != 0:
        ref_index = np.argmin(np.abs(elout_res.time - 1))
        ref_def_grad = def_grad[ref_index, :, :].reshape(-1, 3, 3).swapaxes(1, 2)
        ref_def_grad_inv = np.linalg.inv(ref_def_grad)

    # todo: vectorization
    for i_time in range(len(elout_res.time)):
        for i_ele in range(len(ele_id)):
            if refrence_time != 0:
                # absolute deformation gradient
                f = def_grad[i_time, i_ele, :].reshape(3, 3).T
                # deformation gradient based on the reference time
                ff = np.matmul(ref_def_grad_inv[i_ele], f)
                right_cauchy_green = np.matmul(ff.T, ff)

            else:
                right_cauchy_green = np.matmul(
                    def_grad[i_time, i_ele, :].reshape(3, 3),
                    def_grad[i_time, i_ele, :].reshape(3, 3).T,
                )

            # Green Lagrangian strain: E = 0.5*(lambda**2-1)
            # lambda = sqrt(e*right_cauchy_green*e)
            strain[0, i_time, i_ele] = 0.5 * (
                np.matmul(np.matmul(e_l[i_ele].T, right_cauchy_green), e_l[i_ele]) - 1
            )
            strain[1, i_time, i_ele] = 0.5 * (
                np.matmul(np.matmul(e_r[i_ele].T, right_cauchy_green), e_r[i_ele]) - 1
            )
            strain[2, i_time, i_ele] = 0.5 * (
                np.matmul(np.matmul(e_c[i_ele].T, right_cauchy_green), e_c[i_ele]) - 1
            )

    return strain


def compute_AHA17_segment_strain(model: HeartModel, elout_res: Elout, element_strain):
    """
    Average elemental strain for AHA17 segments.

    Parameters
    ----------
    model
    elout_res
    element_strain

    Returns
    -------
    strain 3 * time * 17
    """
    aha_strain = np.zeros((3, len(elout_res.time), 17))
    for i in range(1, 18):
        ele_id = np.where(model.aha_ids == i)[0]
        indices = np.where(np.in1d(elout_res.ids, ele_id))[0]
        strain_avg = np.mean(element_strain[:, :, indices], axis=2)
        aha_strain[:, :, i - 1] = strain_avg

    return aha_strain


if __name__ == "__main__":
    model = HeartModel.load_model("heart_model.pickle")
    res = Elout("binout0000")

    # strain = compute_myocardial_strain(model, res,refrence_time=1)
    # np.save('strain',strain);exit()
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
        bullseye_plot(ax[i], aha_strain[i, 2, :])  # at second time step

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
