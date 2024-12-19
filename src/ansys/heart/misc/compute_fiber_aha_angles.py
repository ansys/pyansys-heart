# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
Compute average fiber orientations with respect to UHCs in each AHA region in the LV.

from davoks
"""

import os
import pathlib

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
import pandas as pd
import pyvista as pv

import ansys.heart.core.models as models
from ansys.heart.simulator.simulator import DynaSettings, EPSimulator


def angle_btw_vectors(x, y):
    """Computes angles between N number of M-dimensional vectors
    Input:  shape(x): (N,M)
            shape(y): (N,M)
    Output: shape(): (N,)
    """  # noqa
    return (
        np.arccos(np.clip(np.sum(x * y, axis=1) / (norm(x, axis=1) * norm(y, axis=1)), -1, 1))
        * 180
        / np.pi
    )


def signed_angle_btw_vectors(x, y, n):
    """Computes signed angles between N number of M-dimensional vectors
        Input:  shape(x): (N,M)
                shape(y): (N,M)
                shape(n): (N,M)
        Output: shape(): (N,)

    Va . Vb == |Va| * |Vb| * cos(alpha)    (by definition)
            == |Va| * |Vb| * cos(beta)     (cos(alpha) == cos(-alpha) == cos(360° - alpha)


    Va x Vb == |Va| * |Vb| * sin(alpha) * n1
        (by definition; n1 is a unit vector perpendicular to Va and Vb with
         orientation matching the right-hand rule)

    Therefore (again assuming Vn is normalized):
       n1 . Vn == 1 when beta < 180
       n1 . Vn == -1 when beta > 180

    ==>  (Va x Vb) . Vn == |Va| * |Vb| * sin(beta)

    """  # noqa
    n /= norm(n, axis=1)[:, None]  # normalize normal vectors

    return np.arctan2(np.sum(np.cross(x, y) * n, axis=1), np.sum(x * y, axis=1)) * 180 / np.pi


def compute_aha_fiber_angles(mesh: pv.UnstructuredGrid, out_dir):
    """Computes average fiber inclination in each AHA region
    as a function of the transmural depth
    """  # noqa
    assert "aha17" in mesh.cell_data
    assert "fiber" in mesh.cell_data
    assert "sheet" in mesh.cell_data
    assert "transmural" in mesh.point_data
    assert "rotational" in mesh.point_data

    aha_ids = mesh["aha17"]
    aha_elements = np.where(~np.isnan(aha_ids))[0]
    aha_model = mesh.extract_cells(aha_elements)
    aha_model.cell_data["AHA"] = aha_ids[aha_elements]

    # print("Nmbr of points:\t{:d}".format(mesh.n_points))
    # print("Nmbr of cells:\t{:d}".format(mesh.n_cells))
    # print("Nmbr of AHA cells:\t{:d}".format(len(aha_elements)))

    # load fibers and sheets at cells
    el_fibers = mesh.cell_data["fiber"]
    el_sheets = mesh.cell_data["sheet"]
    el_fibers = el_fibers[aha_elements]
    el_sheets = el_sheets[aha_elements]
    el_fibers /= np.linalg.norm(el_fibers, axis=1)[:, None]  # make sure fibers are normalized

    # interpolate transmural depth from points to cells
    el_depths = mesh.point_data_to_cell_data()["transmural"][aha_elements]
    el_depths = 2.0 * el_depths - 1.0  # map from [0,1] -> [-1,1]
    # interpolate rotational coordinates from points to cells
    # el_rotat = mesh.point_data_to_cell_data()["rotational"][aha_elements]

    # compute transmural vector from derivative of transmural depth
    # TODO : interpolate t instead of grad_t
    pt_grad_t = mesh.compute_derivative(scalars="transmural", preference="cell")["gradient"]
    mesh.point_data.set_scalars(name="grad_t", scalars=pt_grad_t)
    el_grad_t = mesh.point_data_to_cell_data()["grad_t"][aha_elements]
    el_grad_t = el_grad_t / np.linalg.norm(el_grad_t, axis=1)[:, None]

    # compute transmural vector from derivative of rotational coordinate
    pt_grad_r = mesh.compute_derivative(scalars="rotational", preference="cell")["gradient"]
    mesh.point_data.set_scalars(name="grad_r", scalars=pt_grad_r)
    el_grad_r = mesh.point_data_to_cell_data()["grad_r"][aha_elements]
    # set elements at the rotational coordinate discontinuity to NaN
    # TODO: prescribe to average gradient from nearest neighbors
    norm_grad_r = norm(el_grad_r, axis=1)
    id_discont = np.where(norm_grad_r > np.average(norm_grad_r) + 2 * np.std(norm_grad_r))[0]
    nans = np.empty((3,))
    nans[:] = np.nan
    el_grad_r[id_discont, :] = nans
    el_grad_r = el_grad_r / np.linalg.norm(el_grad_r, axis=1)[:, None]

    # compute angle between fibers and transmural vectors
    el_angles_t = signed_angle_btw_vectors(el_fibers, el_grad_t, el_grad_r)

    # compute angle between fibers and rotational vectors
    el_angles_r = signed_angle_btw_vectors(el_grad_r, el_fibers, el_grad_t)

    # get aha17 label for left ventricle elements
    aha17_label = aha_ids[aha_elements]

    # average fiber angles wrt to transmural and rotational coords in each AHA and depth
    ndepths = 9
    depth_bin_edges = np.linspace(-1.0, 1.0, ndepths + 1)
    depth_bin_ctrs = 0.5 * (depth_bin_edges[:-1] + depth_bin_edges[1:])
    el_angles_t_avg = np.zeros((len(aha_elements)))
    aha_angles_t = np.zeros((17, ndepths))
    el_angles_r_avg = np.zeros((len(aha_elements)))
    aha_angles_r = np.zeros((17, ndepths))
    for i in range(1, 18):
        idx_aha = aha17_label == i
        for j in range(ndepths):
            depth_min = depth_bin_edges[j]
            depth_max = depth_bin_edges[j + 1]
            idx_depth = (el_depths >= depth_min) & (el_depths < depth_max)
            idx = np.where(idx_aha & idx_depth)[0]
            aha_angles_t[i - 1, j] = np.nanmean(el_angles_t[idx])
            el_angles_t_avg[idx] = aha_angles_t[i - 1, j]
            aha_angles_r[i - 1, j] = np.nanmean(el_angles_r[idx])
            el_angles_r_avg[idx] = aha_angles_r[i - 1, j]

    # save to vtk
    aha_model.cell_data["transmural_angles"] = el_angles_t
    aha_model.cell_data["transmural_angles_aha"] = el_angles_t_avg
    aha_model.cell_data["rotational_angles"] = el_angles_r
    aha_model.cell_data["rotational_angles_aha"] = el_angles_r_avg
    aha_model.cell_data["transmural_depth"] = el_depths
    aha_model.cell_data.set_vectors(el_fibers, "fibers", deep_copy=True)
    aha_model.cell_data.set_vectors(el_sheets, "sheets", deep_copy=True)
    aha_model.cell_data.set_vectors(el_grad_t, "grad_t", deep_copy=True)
    aha_model.cell_data.set_vectors(el_grad_r, "grad_r", deep_copy=True)
    aha_model.save(os.path.join(pathlib.Path(out_dir), "aha_averaged_angles.vtk"))

    # save to csv
    cols = ["{:1.2f}".format(x) for x in depth_bin_ctrs]
    rows = ["{:d}".format(x) for x in range(1, 18)]
    df_r = pd.DataFrame(data=aha_angles_r, index=rows, columns=cols)
    df_r.to_csv(os.path.join(out_dir, "AHA_fiber_angles_r.csv"), index=True)
    df_t = pd.DataFrame(data=aha_angles_t, index=rows, columns=cols)
    df_t.to_csv(os.path.join(out_dir, "AHA_fiber_angles_t.csv"), index=True)

    # print(df_r)
    # print(df_t)

    return df_r, df_t


def plot_fiber_aha_angles(data: pd.DataFrame | str):
    """
    Plot average fiber helical orientation in each AHA as a function of transmural depth

    """  # noqa
    if isinstance(data, str):
        df = pd.read_csv(os.path.join(data, "AHA_fiber_angles_t.csv"))
    elif isinstance(data, pd.DataFrame):
        df = data

    aha_names = [
        "Basal anterior",
        "Basal anteroseptal",
        "Basal inferoseptal",
        "Basal inferior",
        "Basal inferolateral",
        "Basal anterolateral",
        "Mid anterior",
        "Mid anteroseptal",
        "Mid inferoseptal",
        "Mid inferior",
        "Mid inferolateral",
        "Mid anterolateral",
        "Apical anterior",
        "Apical inferior",
        "Apical septal",
        "Apical lateral",
        "Apex",
    ]

    fig, axs = plt.subplots(3, 6)
    fig.set_size_inches(31, 18)
    # print(df.columns.tolist()[1:])
    depths = [float(x) for x in df.columns.tolist()[1:]]
    for iaha in range(17):
        i = iaha // 6
        j = iaha % 6
        alphas = df.iloc[iaha].tolist()[1:]
        # print(alphas)
        axs[i, j].plot(depths, alphas, "o-b")
        axs[i, j].plot([-1, 1], [-60, -60], "--k")
        axs[i, j].plot([-1, 1], [60, 60], "--k")
        axs[i, j].set_title(aha_names[iaha])
        axs[i, j].set_xlim(xmin=-1, xmax=1)
        axs[i, j].set_ylim(ymin=-100, ymax=100)

    for ax in axs.flat:
        ax.set(xlabel="Transmural Depth", ylabel="$\\alpha_h$")

    for ax in axs.flat:
        ax.label_outer()

    axs[2, 5].remove()

    # plt.savefig("fiber_helical_angles.png", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    """Demo
    """
    # assumed LS-DYNA
    lsdyna_path = r"D:\ansysdev\lsdyna\ls-dyna_smp_d_DEV_112901-gcbb8e36701_winx64_ifort190.exe"
    dyna_settings = DynaSettings(
        lsdyna_path=lsdyna_path, dynatype="smp", num_cpus=4, platform="windows"
    )

    # assumed case
    workdir = os.path.abspath(os.path.join("downloads", "Rodero2021", "01", "BV"))
    path_to_model = os.path.join(workdir, "heart_model.vtu")

    model: models.BiVentricle = models.BiVentricle(workdir)
    model.load_model_from_mesh(path_to_model, path_to_model.replace(".vtu", ".partinfo.json"))

    # get fields data
    simulator = EPSimulator(
        model=model,
        dyna_settings=dyna_settings,
        simulation_directory=os.path.join(workdir, "simulation-EP"),
    )
    simulator.settings.load_defaults()
    simulator.compute_uhc()
    simulator.compute_fibers()

    # get AHA labels
    from ansys.heart.core.helpers.landmarks import compute_aha17

    aha_ids = compute_aha17(simulator.model, simulator.model.short_axis, simulator.model.l4cv_axis)
    grid = simulator.model.mesh.extract_cells_by_type(10)
    grid.cell_data["aha17"] = aha_ids

    # compute angles
    df_r, df_t = compute_aha_fiber_angles(grid, workdir)
    plot_fiber_aha_angles(df_t)
