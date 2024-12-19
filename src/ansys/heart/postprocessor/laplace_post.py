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

import os

import numpy as np
import pyvista as pv

from ansys.heart.core import LOG as LOGGER
from ansys.heart.postprocessor.dpf_utils import D3plotReader
from ansys.heart.simulator.settings.settings import AtrialFiber


def read_temperature_field(directory: str, field_list: list[str]) -> pv.UnstructuredGrid:
    """Read thermal fields from d3plot.

    Parameters
    ----------
    directory : str
        directory of d3plot files
    field_list : list[str]
        name of each d3plot file/field

    Returns
    -------
    pv.UnstructuredGrid
        grid with point data of each field
    """
    data = D3plotReader(os.path.join(directory, field_list[0] + ".d3plot"))
    grid: pv.UnstructuredGrid = data.model.metadata.meshed_region.grid

    for name in field_list:
        data = D3plotReader(os.path.join(directory, name + ".d3plot"))
        t = data.model.results.temperature.on_last_time_freq.eval()[0].data
        if len(t) == grid.n_points:
            t = t
        elif len(t) == 3 * grid.n_points:
            LOGGER.warning(
                "DPF reads temperature as a vector field, but expecting a scalar field.\
                Consider updating the DPF server."
            )
            t = t[::3]
        else:
            LOGGER.error("Failed to read d3plot.")
            exit()

        grid.point_data[name] = np.array(t, dtype=float)

    return grid


def compute_cell_gradient(grid: pv.UnstructuredGrid) -> pv.UnstructuredGrid:
    """Compute cell gradient.

    Parameters
    ----------
    grid : pv.UnstructuredGrid
        grid with point data associated

    Returns
    -------
    pv.UnstructuredGrid
        grid with gradient vectors in cell data
    """
    grid2 = grid.point_data_to_cell_data()
    for name in grid2.cell_data.keys():
        derivative = grid2.compute_derivative(scalars=name, preference="cell")
        res = derivative["gradient"]
        grid2["grad_" + name] = res

    return grid2


def update_transmural_by_normal(grid: pv.UnstructuredGrid, surface: pv.PolyData) -> np.ndarray:
    """Use surface normal for transmural direction.

    Note
    ----
    Assume mesh is coarse compared to the thinkness, solid cell normal
    is interpolated from closest surface normal

    Parameters
    ----------
    grid : pv.UnstructuredGrid
        atrium grid
    surface : pv.PolyData
        atrium endocardium surface

    Returns
    -------
    np.ndarray
        cell transmural direction vector
    """
    surface_normals = surface.clean().compute_normals()

    from scipy import spatial

    tree = spatial.cKDTree(surface_normals.cell_centers().points)

    cell_center = grid.cell_centers().points
    d, t = tree.query(cell_center, 1)

    grad_trans = surface_normals.cell_data["Normals"][t]

    return grad_trans


def orthogonalization(
    grad_trans: np.ndarray, k: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Orthogonalization.

    Parameters
    ----------
    grad_trans : np.ndarray
        transmural vector
    k : np.ndarray
        Bundle selection vector

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        local coordinate system e_l,e_n,e_t
    """
    norm = np.linalg.norm(grad_trans, axis=1)
    bad_cells = np.argwhere(norm == 0).ravel()

    LOGGER.debug(
        f"{len(bad_cells)} cells have null gradient in transmural direction."
        f" This should only be at valve regions and can be checked from the vtk file."
    )

    norm = np.where(norm != 0, norm, 1)
    e_t = grad_trans / norm[:, None]

    k_e = np.einsum("ij,ij->i", k, e_t)
    en = k - np.einsum("i,ij->ij", k_e, e_t)
    norm = np.linalg.norm(en, axis=1)
    norm = np.where(norm != 0, norm, 1)
    e_n = en / norm[:, None]

    e_l = np.cross(e_n, e_t)

    return e_l, e_n, e_t


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

    def bundle_selection(grid):
        """Left atrium bundle selection.

        Add two cell data to grid.
        - 'k' is unit vector from different gradient fields.
        - 'bundle' labels regions of selection.

        """
        # bundle selection
        tau_mv = settings.tau_mv  # 0.65
        tau_lpv = settings.tau_lpv  # 0.65
        tau_rpv = settings.tau_rpv  # 0.1

        grid["k"] = np.zeros((grid.n_cells, 3))
        grid["bundle"] = np.zeros(grid.n_cells, dtype=int)

        # MV region
        mask_mv = grid["r"] >= tau_mv
        grid["k"][mask_mv] = grid["grad_r"][mask_mv]
        grid["bundle"][mask_mv] = 1
        # LPV region
        mask = np.invert(mask_mv) & (grid["v"] < tau_lpv)
        grid["k"][mask] = grid["grad_v"][mask]
        grid["bundle"][mask] = 2
        # RPV region
        mask = np.invert(mask_mv) & (grid["v"] > tau_rpv)
        grid["k"][mask] = grid["grad_v"][mask]
        grid["bundle"][mask] = 3

        # rest and assign to grad_ab
        mask = grid["bundle"] == 0
        grid["k"][mask] = grid["grad_ab"][mask]

        return

    data = read_temperature_field(directory, field_list=["trans", "ab", "v", "r"])
    grid = compute_cell_gradient(data)

    if endo_surface is not None:
        grid.cell_data["grad_trans"] = update_transmural_by_normal(grid, endo_surface)

    bundle_selection(grid)

    el, en, et = orthogonalization(grid["grad_trans"], grid["k"])

    grid.cell_data["e_l"] = el
    grid.cell_data["e_n"] = en
    grid.cell_data["e_t"] = et

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

    def bundle_selection(grid):
        """Right atrium bundle selection.

        Add two cell data to grid.
        - 'k' is unit vector from different gradient fields.
        - 'bundle' labels regions of selection.

        """
        tao_tv = settings.tau_tv  # 0.9
        tao_raw = settings.tau_raw  # 0.55
        tao_ct_minus = settings.tau_ct_minus  # -0.18
        tao_ct_plus = settings.tau_ct_plus  # -0.1
        tao_icv = settings.tau_icv  # 0.9
        tao_scv = settings.tau_scv  # 0.1
        tao_ib = settings.tau_ib  # 0.35
        tao_ras = settings.tau_ras  # 0.135

        ab = grid["ab"]
        v = grid["v"]
        r = grid["r"]
        w = grid["w"]

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

        return

    data = read_temperature_field(directory, field_list=["trans", "ab", "v", "r", "w"])
    grid = compute_cell_gradient(data)

    if endo_surface is not None:
        grid.cell_data["grad_trans"] = update_transmural_by_normal(grid, endo_surface)

    bundle_selection(grid)

    el, en, et = orthogonalization(grid["grad_trans"], grid["k"])

    grid.cell_data["e_l"] = el
    grid.cell_data["e_n"] = en
    grid.cell_data["e_t"] = et

    return grid
