# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
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

"""Compute average fiber orientations with respect to UHCs in each AHA region in the LV."""

import os

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
import pandas as pd
import pyvista as pv

from ansys.health.heart.utils.misc import angle_between_vectors


def compute_gradient(mesh: pv.UnstructuredGrid, scalar: str, interpolate_to_cell: bool = False):
    """Compute gradient of a scalar field on the mesh.

    Parameters
    ----------
    mesh : pv.UnstructuredGrid
        Input mesh.
    scalar : str
        Name of the scalar field.

    Returns
    -------
    gradient : np.ndarray
        Gradient of the scalar field at each point.
    """
    if scalar not in mesh.point_data:
        raise ValueError(f"Scalar '{scalar}' not found in mesh point data.")

    # compute gradient at points
    mesh_with_grad = mesh.compute_derivative(scalars=scalar, preference="point")

    if interpolate_to_cell:
        grad = mesh_with_grad.point_data_to_cell_data().cell_data["gradient"]
    else:
        grad = mesh_with_grad.point_data["gradient"]
    return grad


def _get_angles_from_mesh(grid: pv.UnstructuredGrid) -> pd.DataFrame:
    """Retrieve the fiber angles from the mesh and store in pandas dataframe."""
    expected_arrays = ["aha17", "helix_angle", "transverse_angle", "depth_bin"]
    for arr in expected_arrays:
        if arr not in grid.cell_data and arr not in grid.point_data:
            raise ValueError(f"Expected array '{arr}' not found in mesh data.")

    # Retrieve angles for each AHA segment and depth bin
    ndepths = np.unique(grid["depth_bin"]).shape[0]
    n_ahas = np.unique(grid["aha17"]).shape[0]
    aha_ids = np.unique(grid["aha17"])

    data_helix = np.full((ndepths, n_ahas), np.nan)
    data_transverse = np.full((ndepths, n_ahas), np.nan)

    # Retrieve helix and transverse angles for each AHA segment and depth bin
    for ii in range(1, 18):
        for jj in range(1, ndepths + 1):
            mask = grid["aha17"] == ii
            mask &= grid["depth_bin"] == jj

            helix_angle = np.nanmean(grid.extract_cells(mask)["helix_angle"])
            transverse_angle = np.nanmean(grid.extract_cells(mask)["transverse_angle"])

            data_helix[jj - 1, ii - 1] = helix_angle  # depth in row, AHA in column
            data_transverse[jj - 1, ii - 1] = transverse_angle  # depth in row, AHA in column

    legend_entries = [f"AHA_{i}" for i in aha_ids]

    return data_helix, data_transverse, legend_entries


def compute_aha_fiber_angles1(mesh: pv.UnstructuredGrid, out_dir: str):
    """Compute the average fiber helix and transverse angles."""
    expected_arrays = ["aha17", "fiber", "sheet", "transmural", "rotational"]
    for arr in expected_arrays:
        if arr not in mesh.cell_data and arr not in mesh.point_data:
            raise ValueError(f"Expected array '{arr}' not found in mesh data.")

    aha_ids = mesh["aha17"]
    aha_elements = np.where(~np.isnan(aha_ids))[0]
    aha_model: pv.UnstructuredGrid = mesh.extract_cells(aha_elements)
    # aha_model.cell_data["fiber"] = aha_model.cell_data["fiber"] * -1  # flip fiber direction
    aha_ids = aha_model.cell_data["aha17"]

    # flip fibers
    # aha_model.cell_data["fiber"] = 1 * aha_model.cell_data["fiber"]

    # load fibers and sheets at cells
    el_fibers = aha_model.cell_data["fiber"]
    el_sheets = aha_model.cell_data["sheet"]  # noqa N841

    # TODO : interpolate t instead of grad_t
    aha_model.cell_data["grad_transmural"] = compute_gradient(
        aha_model, "transmural", interpolate_to_cell=True
    )
    el_grad_t = aha_model.cell_data["grad_transmural"]

    # compute rotational vector from derivative of rotational coordinate
    aha_model.cell_data["grad_rotational"] = compute_gradient(
        aha_model, "rotational", interpolate_to_cell=True
    )
    el_grad_r = aha_model.cell_data["grad_rotational"]

    # set elements at the rotational coordinate discontinuity to NaN
    # TODO: prescribe to average gradient from nearest neighbors
    norm_grad_r = norm(el_grad_r, axis=1)
    id_discont = np.where(norm_grad_r > np.average(norm_grad_r) + 2 * np.std(norm_grad_r))[0]
    nans = np.empty((3,))
    nans[:] = np.nan
    el_grad_r[id_discont, :] = nans
    el_grad_r = el_grad_r / np.linalg.norm(el_grad_r, axis=1)[:, None]

    aha_model.cell_data["grad_rotational"] = el_grad_r

    # visualize: rotational vector
    # plotter = pv.Plotter()
    # aha_model.set_active_scalars("aha17")
    # glyph = aha_model.glyph("grad_transmural", scale=False)
    # plotter.add_mesh(aha_model, opacity=0.1, color="white")
    # plotter.add_mesh(glyph, color="blue")
    # plotter.show()

    # compute longitudinal vector as cross product of transmural and rotational vectors
    el_grad_l = np.cross(el_grad_t, el_grad_r)
    aha_model.cell_data["grad_longitudinal"] = el_grad_l

    # compute angle between rotational and fiber vectors in plane defined
    # by longitudinal vector (transverse angle)

    el_angles_t = angle_between_vectors(el_grad_r, el_fibers, el_grad_l, "2-quadrant-signed")
    aha_model.cell_data["transverse_angle"] = el_angles_t

    # compute angle between fibers and rotational vectors in plane defined
    # by transmural vector (helix angle)
    el_angles_r = angle_between_vectors(el_grad_r, el_fibers, el_grad_t, "2-quadrant-signed")
    aha_model.cell_data["helix_angle"] = el_angles_r

    # Add depth bins to the mesh
    ndepths = 9
    depth_lower_limit = np.linspace(-1.01, 1.0, ndepths + 1)[0:-1]
    depth_upper_limit = np.append(depth_lower_limit[1:], 1.01)
    depth_bin_ctrs = np.mean(np.vstack([depth_lower_limit, depth_upper_limit]), axis=0)

    # Compute the depth of each cell in the mesh, normalized to [-1 (endocardium), 1 (epicardium)]
    depth_cells = aha_model.point_data_to_cell_data()["transmural"]
    aha_model.cell_data["depth"] = 2 * depth_cells - 1.0

    # Split in 9 transmural bins and store bin number in cell data
    aha_model.cell_data["depth_bin"] = 0.0
    depth_bin = 1
    for lower, upper in zip(depth_lower_limit, depth_upper_limit):
        mask = aha_model["depth"] > lower
        mask &= aha_model["depth"] < upper

        aha_model.cell_data["depth_bin"][mask] = float(depth_bin)
        depth_bin += 1

    if np.any(aha_model.cell_data["depth_bin"] == 0):
        raise ValueError("Some cells were not assigned to a depth bin - check tolerances.")

    data_helix, data_transverse, legend_entries = _get_angles_from_mesh(aha_model)

    # save to csv
    cols = ["{:1.2f}".format(x) for x in depth_bin_ctrs]
    rows = ["{:d}".format(x) for x in range(1, 18)]
    df_helix_angle = pd.DataFrame(data=data_helix.T, index=rows, columns=cols)
    df_helix_angle.to_csv(os.path.join(out_dir, "AHA_fiber_angles_r.csv"), index=True)
    df_transverse_angle = pd.DataFrame(data=data_transverse.T, index=rows, columns=cols)
    df_transverse_angle.to_csv(os.path.join(out_dir, "AHA_fiber_angles_t.csv"), index=True)

    # fig, axs = plt.subplots(1, 2)
    # axs[0].plot(depth_bin_ctrs, data_helix, "o-")
    # axs[0].set_title(r"Helix Angle alpha_h")
    # axs[0].set_xlim([-1, 1])
    # axs[0].set_xlabel("transmural depth")
    # axs[0].set_ylabel("angle [deg]")
    # axs[1].plot(depth_bin_ctrs, data_transverse, "o-")
    # axs[1].set_title(r"Transverse Angle alpha_t")
    # axs[1].set_xlabel("transmural depth [-]")
    # axs[1].set_xlim([-1, 1])
    # axs[0].legend(legend_entries)
    # fig.show()

    # aha_model.save(os.path.join(out_dir, "aha_model.vtu"))

    # sub_model: pv.UnstructuredGrid = aha_model.threshold([1, 1], "aha17").threshold(
    #     [1, 1], "depth_bin"
    # )

    # # plot components: transmural, rotational and longitudinal vectors
    # plotter = pv.Plotter()
    # plotter.add_mesh(sub_model.glyph("grad_transmural", scale=False), color="b")
    # plotter.add_mesh(sub_model.glyph("grad_rotational", scale=False), color="r")
    # plotter.add_mesh(sub_model.glyph("grad_longitudinal", scale=False), color="g")
    # plotter.add_mesh(sub_model.glyph("fiber", scale=False), color="white")
    # plotter.add_mesh(aha_model, opacity=0.5, color="white")
    # plotter.show()

    return df_helix_angle, df_transverse_angle


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
