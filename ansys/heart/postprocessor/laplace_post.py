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

"""Post process script related to Laplace solving (UHC, fibers)."""

import copy
import os
from typing import List

from ansys.heart.postprocessor.dpf_utils import D3plotReader
import numpy as np
import pyvista as pv

from ansys.heart.core import LOG as LOGGER
from ansys.heart.simulator.settings.settings import AtrialFiber


def read_uvc(directory):
    """Read UVC from d3plot files."""
    data = D3plotReader(os.path.join(directory, "apico-basal.d3plot"))
    grid: pv.UnstructuredGrid = data.model.metadata.meshed_region.grid

    for field_name in ["apico-basal", "transmural", "rotational"]:
        data = D3plotReader(os.path.join(directory, field_name + ".d3plot"))
        t = data.model.results.temperature.on_last_time_freq.eval()[0].data
        if len(t) == grid.n_points:
            t = t
        elif len(t) == 3 * grid.n_points:
            t = t[::3]
        else:
            LOGGER.error("Failed to read d3plot in UVC")
            exit()
        grid.point_data[field_name] = np.array(t, dtype=float)

    return grid


def orthogonalization(grid) -> pv.UnstructuredGrid:
    """Gram–Schmidt orthogonalization."""
    norm = np.linalg.norm(grid["grad_trans"], axis=1)
    bad_cells = np.argwhere(norm == 0).ravel()

    LOGGER.debug(
        f"{len(bad_cells)} cells have null gradient in transmural direction."
        f" This should only be at valve regions and can be checked from the vtk file."
    )

    # grad_trans_new = grid["grad_trans"]
    # grad_trans_new[bad_cells] = np.array([[1, 0, 0]]).repeat(len(bad_cells), axis=0)
    # norm = np.linalg.norm(grad_trans_new, axis=1)
    # grid["e_t"] = grad_trans_new / norm[:, None]

    norm = np.where(norm != 0, norm, 1)
    grid["e_t"] = grid["grad_trans"] / norm[:, None]

    k_e = np.einsum("ij,ij->i", grid["k"], grid["e_t"])
    en = grid["k"] - np.einsum("i,ij->ij", k_e, grid["e_t"])
    norm = np.linalg.norm(en, axis=1)
    norm = np.where(norm != 0, norm, 1)
    grid["e_n"] = en / norm[:, None]

    grid["e_l"] = np.cross(grid["e_n"], grid["e_t"])
    return grid


def get_gradient(directory, field_list: List[str]) -> pv.UnstructuredGrid:
    """Read thermal fields from d3plot and compute gradient."""
    data = D3plotReader(os.path.join(directory, field_list[0] + ".d3plot"))
    grid: pv.UnstructuredGrid = data.model.metadata.meshed_region.grid

    for name in field_list:
        data = D3plotReader(os.path.join(directory, name + ".d3plot"))
        t = data.model.results.temperature.on_last_time_freq.eval()[0].data
        if len(t) == grid.n_points:
            t = t
        elif len(t) == 3 * grid.n_points:
            t = t[::3]
        else:
            LOGGER.error("Failed to read d3plot.")
            exit()

        # force a deep copy
        grid.point_data[name] = copy.deepcopy(t)

    # note vtk gradient method shows warning/error for some cells
    grid2 = grid.point_data_to_cell_data()
    for name in field_list:
        derivative = grid2.compute_derivative(scalars=name, preference="cell")
        res = derivative["gradient"]
        grid2["grad_" + name] = res

    grid2.save("gradient.vtk")

    return grid2


def _update_trans_by_normal(grid: pv.UnstructuredGrid, surface: pv.PolyData):
    """Use surface normal to replace gradient of transmural direction."""
    with_normals = surface.clean().compute_normals()

    from scipy import spatial

    tree = spatial.cKDTree(with_normals.cell_centers().points)

    cell_center = grid.cell_centers().points
    d, t = tree.query(cell_center, 1)
    # print(max(d))
    grid["grad_trans"] = with_normals.cell_data["Normals"][t]

    return grid


def compute_la_fiber_cs(
    directory: str, settings: AtrialFiber, endo_surface: pv.PolyData = None
) -> pv.UnstructuredGrid:
    """Compute left atrium fibers coordinate system.

    Parameters
    ----------
    directory : str
        directory of d3plot files.
    settings : AtrialFiber
        Atrial fiber settings.
    endo_surface : pv.PolyData, optional
        _description_, by default None
        If given, normal direction will be updated by surface normal instead of Laplace solution.

    Notes
    -----
    Method descrbed in https://doi.org/10.1016/j.cma.2020.113468

    Returns
    -------
    pv.UnstructuredGrid
        pv object with fiber coordinates system.
    """

    def bundle_selection(grid) -> pv.UnstructuredGrid:
        """Left atrium bundle selection."""
        # grid = pv.read(os.path.join(directory, "res.vtk"))

        # bundle selection
        tau_mv = settings.tau_mv  # 0.65
        tau_lpv = settings.tau_lpv  # 0.65
        tau_rpv = settings.tau_rpv  # 0.1

        grid["k"] = np.zeros((grid.n_cells, 3))
        grid["bundle"] = np.zeros(grid.n_cells, dtype=int)

        mask_mv = grid["r"] >= tau_mv
        grid["k"][mask_mv] = grid["grad_r"][mask_mv]
        grid["bundle"][mask_mv] = 1

        mask = np.invert(mask_mv) & (grid["v"] < tau_lpv)
        grid["k"][mask] = grid["grad_v"][mask]
        grid["bundle"][mask] = 2

        mask = np.invert(mask_mv) & (grid["v"] > tau_rpv)
        grid["k"][mask] = grid["grad_v"][mask]
        grid["bundle"][mask] = 3

        mask = grid["bundle"] == 0
        grid["k"][mask] = grid["grad_ab"][mask]

        # grid.save(os.path.join(directory, "res2.vtk"))

        return grid

    grid = get_gradient(directory, field_list=["trans", "ab", "v", "r"])
    # TODO sometimes, pv object broken when pass directly

    grid = pv.read("gradient.vtk")
    if endo_surface is not None:
        grid = _update_trans_by_normal(grid, endo_surface)

    grid = bundle_selection(grid)

    grid.save("la_fiber.vtk")
    grid = pv.read("la_fiber.vtk")

    grid = orthogonalization(grid)
    grid.save("la_fiber.vtk")

    return grid


def compute_ra_fiber_cs(
    directory: str, settings: AtrialFiber, endo_surface: pv.PolyData = None
) -> pv.UnstructuredGrid:
    """Compute right atrium fibers coordinate system.

    Parameters
    ----------
    directory : str
        directory of d3plot files.
    settings : AtrialFiber
        Atrial fiber settings.
    endo_surface : pv.PolyData, optional
        _description_, by default None
        If given, normal direction will be updated by surface normal instead of Laplace solution.

    Notes
    -----
    Method descrbed in https://doi.org/10.1016/j.cma.2020.113468

    Returns
    -------
    pv.UnstructuredGrid
        pv object with fiber coordinates system.
    """

    def bundle_selection(grid) -> pv.UnstructuredGrid:
        """Left atrium bundle selection."""
        tao_tv = settings.tau_tv  # 0.9
        tao_raw = settings.tau_raw  # 0.55
        tao_ct_minus = settings.tau_ct_minus  # -0.18
        tao_ct_plus = settings.tau_ct_plus  # -0.1
        tao_icv = settings.tau_icv  # 0.9
        tao_scv = settings.tau_scv  # 0.1
        tao_ib = settings.tau_ib  # 0.35
        tao_ras = settings.tau_ras  # 0.135

        trans = grid["trans"]
        ab = grid["ab"]
        v = grid["v"]
        r = grid["r"]
        w = grid["w"]

        trans_grad = grid["grad_trans"]
        ab_grad = grid["grad_ab"]
        v_grad = grid["grad_v"]
        r_grad = grid["grad_r"]
        w_grad = grid["grad_w"]
        tag = np.zeros(ab.shape)
        k = np.zeros(ab_grad.shape)

        tv = 1
        icv = 2
        scv = 3
        raw = 4
        ct = 5
        ib = 6
        ras_top = 7
        ras_center = 9
        ras_bottom = 10
        raw_ist_raa = 8

        for i in range(grid.n_cells):
            if r[i] >= tao_tv:
                k[i] = r_grad[i]
                tag[i] = tv
            else:
                if r[i] < tao_raw:
                    if tao_ct_minus <= w[i] <= tao_ct_plus:
                        k[i] = w_grad[i]
                        tag[i] = ct
                    elif w[i] < tao_ct_minus:
                        if v[i] >= tao_icv or v[i] <= tao_scv:
                            k[i] = v_grad[i]
                            if v[i] >= tao_icv:
                                tag[i] = icv
                            if v[i] <= tao_scv:
                                tag[i] = scv
                        else:
                            k[i] = ab_grad[i]
                            tag[i] = raw
                    else:
                        if v[i] >= tao_icv or v[i] <= tao_scv:
                            k[i] = v_grad[i]
                            if v[i] >= tao_icv:
                                tag[i] = icv
                            if v[i] <= tao_scv:
                                tag[i] = scv
                        else:
                            if w[i] < tao_ib:
                                k[i] = v_grad[i]
                                tag[i] = ib
                            elif w[i] > tao_ras:
                                k[i] = r_grad[i]
                                tag[i] = ras_center
                            else:
                                k[i] = w_grad[i]
                                tag[i] = ras_top
                else:
                    if v[i] >= tao_icv or v[i] <= tao_scv:
                        k[i] = v_grad[i]
                        if v[i] >= tao_icv:
                            tag[i] = icv
                        if v[i] <= tao_scv:
                            tag[i] = scv
                    else:
                        if w[i] >= 0:
                            k[i] = r_grad[i]
                            tag[i] = ras_bottom
                        else:
                            k[i] = ab_grad[i]
                            tag[i] = raw_ist_raa

        grid["k"] = k
        grid["bundle"] = tag.astype(int)

        return grid

    grid = get_gradient(directory, field_list=["trans", "ab", "v", "r", "w"])
    grid = pv.read("gradient.vtk")

    if endo_surface is not None:
        grid = _update_trans_by_normal(grid, endo_surface)

    grid = bundle_selection(grid)

    grid.save("ra_fiber.vtk")
    grid = pv.read("ra_fiber.vtk")

    grid = orthogonalization(grid)
    grid.save("ra_fiber.vtk")

    return grid
