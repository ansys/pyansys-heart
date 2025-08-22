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

"""Postprocessing script related to Laplace solving (UHC, fibers)."""

import os

from deprecated import deprecated
import numpy as np
import pyvista as pv

from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.exceptions import D3PlotNotSupportedError
from ansys.health.heart.post.dpf_utils import D3plotReader
from ansys.health.heart.settings.settings import AtrialFiber


def read_laplace_solution(
    directory: str, field_list: list[str], read_heatflux: bool = False
) -> pv.UnstructuredGrid:
    """Read laplace fields from d3plot files.

    Parameters
    ----------
    directory : str
        Directory of d3plot files.
    field_list : list[str]
        Name of each d3plot file/field.
    read_heatflux : bool, default: False
        Whether to read heatflux.

    Returns
    -------
    pv.UnstructuredGrid
        Grid with point data of each field.
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
                "DPF reads temperature as a vector field but is expecting a scalar field.\
                Consider updating the DPF server."
            )
            t = t[::3]
        else:
            LOGGER.error("Failed to read d3plot.")
            raise D3PlotNotSupportedError("Failed to read d3plot.")

        grid.point_data[name] = np.array(t, dtype=float)

        if read_heatflux:
            last_step = data.model.metadata.time_freq_support.n_sets
            grid.point_data["grad_" + name] = -data.get_heatflux(last_step)

    return grid.copy()


@deprecated(reason="Transmural direction can be automatically read by d3plot heat flux.")
def update_transmural_by_normal(grid: pv.UnstructuredGrid, surface: pv.PolyData) -> np.ndarray:
    """Use surface normal for transmural direction.

    Notes
    -----
    Assume mesh is coarse compared to the thickness. Solid cell normal
    is interpolated from closest surface normal.

    Parameters
    ----------
    grid : pv.UnstructuredGrid
        Atrium grid.
    surface : pv.PolyData
        Atrium endocardium surface.

    Returns
    -------
    np.ndarray
        Cell transmural direction vector.
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
        Transmural vector.
    k : np.ndarray
        Bundle selection vector.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        Local coordinate system ``e_l, e_n, e_t``.
    """
    norm = np.linalg.norm(grad_trans, axis=1)
    bad_cells = np.argwhere(norm == 0).ravel()

    LOGGER.debug(
        f"{len(bad_cells)} cells have null gradient in transmural direction."
        f" This should only be at valve regions and can be checked from the VTK file."
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
        Directory of d3plot files.
    settings : AtrialFiber
        Atrial fiber settings.
    endo_surface : pv.PolyData, default: None
        ``_description_``. If given, normal direction is updated by the surface
        normal instead of the Laplace solution.

    Notes
    -----
    This method is described in `Modeling cardiac muscle fibers in ventricular and
    atrial electrophysiology simulations <https://doi.org/10.1016/j.cma.2020.113468>`_.

    Returns
    -------
    pv.UnstructuredGrid
        PV object with fiber coordinates system.
    """

    def bundle_selection(grid):
        """Left atrium bundle selection.

        Add two-cell data to grid.
        - 'k' is the unit vector from different gradient fields.
        - 'bundle' labels the regions of selection.

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

    solutions = ["trans", "ab", "v", "r"]
    data = read_laplace_solution(directory, field_list=solutions, read_heatflux=True)
    grid = data.point_data_to_cell_data()

    if endo_surface is not None:
        grid.cell_data["grad_trans"] = update_transmural_by_normal(grid, endo_surface)

    bundle_selection(grid)

    el, en, et = orthogonalization(grid["grad_trans"], grid["k"])

    grid.cell_data["e_l"] = el
    grid.cell_data["e_n"] = en
    grid.cell_data["e_t"] = et

    return grid.copy()


def compute_ra_fiber_cs(
    directory: str, settings: AtrialFiber, endo_surface: pv.PolyData = None
) -> pv.UnstructuredGrid:
    """Compute right atrium fibers coordinate system.

    Parameters
    ----------
    directory : str
        Directory of d3plot files.
    settings : AtrialFiber
        Atrial fiber settings.
    endo_surface : pv.PolyData, default: None
        ``_description_``. If given, normal direction is updated by the surface normal
        instead of the Laplace solution.

    Notes
    -----
    This method is described in `Modeling cardiac muscle fibers in ventricular and
    atrial electrophysiology simulations <https://doi.org/10.1016/j.cma.2020.113468>`_.

    Returns
    -------
    pv.UnstructuredGrid
        PV object with the fiber coordinates system.
    """

    def bundle_selection(grid):
        """Right atrium bundle selection.

        Add two-cell data to grid.
        - 'k' is the unit vector from different gradient fields.
        - 'bundle' labels the regions of selection.

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

    solution = ["trans", "ab", "v", "r", "w"]
    data = read_laplace_solution(directory, field_list=solution, read_heatflux=True)
    grid = data.point_data_to_cell_data()

    if endo_surface is not None:
        grid.cell_data["grad_trans"] = update_transmural_by_normal(grid, endo_surface)

    bundle_selection(grid)

    el, en, et = orthogonalization(grid["grad_trans"], grid["k"])

    grid.cell_data["e_l"] = el
    grid.cell_data["e_n"] = en
    grid.cell_data["e_t"] = et

    return grid.copy()


def set_rotation_bounds(
    w: np.ndarray, endo: float, epi: float, outflow_tracts: list[float, float] = None
) -> tuple[np.ndarray, np.ndarray]:
    """Define rotation bounds from input parameters.

    Parameters
    ----------
    w : np.ndarray
        Intra-ventricular interpolation weight if ``outflow_tracts`` is not ``None``.
    endo : float
        Rotation angle at endocardium.
    epi : float
        Rotation angle at epicardium.
    outflow_tracts : list[float, float], default: None
        Rotation angle of enendocardium do and epicardium on outflow tract.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Cell-wise rotation bounds for endocardium and epicardium.
    """

    def _sigmoid(z):
        """Sigmoid function."""
        return 1 / (1 + np.exp(-z))

    if outflow_tracts is not None:
        # rescale w with sigmoid function so it affects only on outflow tracts region
        c0 = 0.3  # clip point
        c1 = 20  # steepness
        w = _sigmoid((w - c0) * c1)
        ro_endo = w * endo + (1 - w) * outflow_tracts[0]
        ro_epi = w * epi + (1 - w) * outflow_tracts[1]
    else:
        # constant rotation angle
        ro_endo = np.ones(len(w)) * endo
        ro_epi = np.ones(len(w)) * epi

    return ro_endo, ro_epi


def compute_rotation_angle(
    grid: pv.UnstructuredGrid,
    w: np.ndarray,
    rotation: list[float, float],
    outflow_tracts: list[float, float] = None,
) -> np.ndarray:
    """Rotate by alpha and beta angles.

    Parameters
    ----------
    grid : pv.UnstructuredGrid
        Mesh grid.
    w : np.ndarray
        Intral ventricular interpolation weight.
    rotation : list[float, float]
        Rotation angles in degrees at endocardium and epicardium.
    outflow_tracts : list[float, float], default: None
        Rotation angle of enendocardium do and epicardium on outflow tract.

    Returns
    -------
    np.ndarray
        Cell-wise rotation angles.

    Notes
    -----
    Compute for all cells, but filter by left/right mask outside of this function.
    """
    rot_endo, rot_epi = set_rotation_bounds(w, rotation[0], rotation[1], outflow_tracts)

    # interpolate along transmural direction
    angle = np.zeros(grid.n_cells)
    angle = rot_epi * (np.ones(grid.n_cells) - grid["d"]) + rot_endo * grid["d"]
    return angle


def compute_ventricle_fiber_by_drbm(
    directory: str,
    settings: dict = {
        "alpha_left": [-60, 60],
        "alpha_right": [-60, 60],
        "alpha_ot": None,
        "beta_left": [-65, 25],
        "beta_right": [-65, 25],
        "beta_ot": None,
    },
    left_only: bool = False,
) -> pv.UnstructuredGrid:
    """Compute the fiber coordinate system from Laplace solving.

    Parameters
    ----------
    directory : str
        Directory of d3plot/tprint files.
    settings : dict, optional
        Rotation angles. By default: ``{ "alpha_left": [-60, 60], "alpha_right": [-60, 60],
        "alpha_ot": None, "beta_left": [-65, 25], "beta_right": [-65, 25], "beta_ot": None, }``.
    left_only : bool, default: False
        Whether to only compute fibers on the left ventricle.

    Notes
    -----
    The D-RBM method is described in `Modeling cardiac muscle fibers in ventricular and
    atrial electrophysiology simulations <https://doi.org/10.1016/j.cma.2020.113468>`_.

    Returns
    -------
    pv.UnstructuredGrid
        Grid containing ``fiber``, ``cross-fiber``, and ``sheet`` vectors.
    """
    solutions = ["trans", "ab_l", "ot_l", "w_l"]
    if not left_only:
        solutions.extend(["ab_r", "ot_r", "w_r", "lr"])

    data = read_laplace_solution(directory, field_list=solutions, read_heatflux=True)
    grid = data.point_data_to_cell_data()

    if left_only:
        # label to 1 for all cells
        left_mask = np.ones(grid.n_cells, dtype=bool)
        grid.cell_data["label"] = np.ones(grid.n_cells, dtype=int)
    else:
        # label to 1 for left ventricle, 2 for right ventricle
        left_mask = grid["lr"] >= 0
        right_mask = grid["lr"] < 0
        label = np.zeros(grid.n_cells, dtype=int)
        label[left_mask] = 1
        label[right_mask] = 2
        grid.cell_data["label"] = label

    # normal direction
    k = np.zeros((grid.n_cells, 3))
    w_l = np.tile(grid["w_l"], (3, 1)).T
    result = w_l * grid["grad_ab_l"] + (np.ones((grid.n_cells, 3)) - w_l) * grid["grad_ot_l"]
    k[left_mask] = result[left_mask]

    if not left_only:
        w_r = np.tile(grid["w_r"], (3, 1)).T
        result = w_r * grid["grad_ab_r"] + (np.ones((grid.n_cells, 3)) - w_r) * grid["grad_ot_r"]
        k[right_mask] = result[right_mask]

    grid.cell_data["k"] = k

    # build local coordinate system
    if not left_only:
        grid.cell_data["grad_trans"][right_mask] *= -1.0  # both LV & RV point to inside

    el, en, et = orthogonalization(grid["grad_trans"], k)

    # normalized transmural distance
    if left_only:
        grid["d"] = grid["trans"]
    else:
        d_l = grid["trans"] / 2
        d_r = np.absolute(grid["trans"])
        grid["d"] = np.zeros(grid.n_cells, dtype=float)
        grid["d"][left_mask] = d_l[left_mask]
        grid["d"][right_mask] = d_r[right_mask]

    # rotation angles for each cell
    alpha = np.zeros(grid.n_cells)
    beta = np.zeros(grid.n_cells)

    alpha[left_mask] = compute_rotation_angle(
        grid, grid["w_l"], settings["alpha_left"], settings["alpha_ot"]
    )[left_mask]
    beta[left_mask] = compute_rotation_angle(
        grid, grid["w_l"], settings["beta_left"], settings["beta_ot"]
    )[left_mask]

    if not left_only:
        alpha[right_mask] = compute_rotation_angle(
            grid, grid["w_r"], settings["alpha_right"], settings["alpha_ot"]
        )[right_mask]
        beta[right_mask] = compute_rotation_angle(
            grid, grid["w_r"], settings["beta_right"], settings["beta_ot"]
        )[right_mask]

    # save data for inspection
    grid.cell_data["alpha"] = alpha
    grid.cell_data["beta"] = beta

    #
    grid.cell_data["fiber"] = np.zeros((grid.n_cells, 3))

    # use f,n,s in Quateroni, it's n, cross fiber
    # use FTS in Bayer, it's S, sheet normal
    grid.cell_data["cross-fiber"] = np.zeros((grid.n_cells, 3))

    # use f,n,s in Quateroni, it's s, sheet
    # use FTS in Bayer, it's T, transverse
    grid.cell_data["sheet"] = np.zeros((grid.n_cells, 3))

    # apply rotation
    for i in range(grid.n_cells):
        q = np.array([el[i], en[i], et[i]]).T
        # rotate alpha around e_t
        a = alpha[i] * np.pi / 180
        rot1 = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
        # rotate beta around e_l
        b = beta[i] * np.pi / 180
        rot2 = np.array([[1, 0, 0], [0, np.cos(b), np.sin(b)], [0, -np.sin(b), np.cos(b)]])
        # apply rotation
        qq = np.matmul(np.matmul(q, rot1), rot2)

        grid.cell_data["fiber"][i] = qq[:, 0]
        grid.cell_data["cross-fiber"][i] = qq[:, 1]
        grid.cell_data["sheet"][i] = qq[:, 2]

    return grid.copy()
