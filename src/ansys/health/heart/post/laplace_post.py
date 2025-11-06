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
from scipy.spatial.transform import Rotation as R  # noqa N817

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


def orthogonalization(e1: np.ndarray, e2: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create a orthonormal coordinate system.

    Parameters
    ----------
    e1 : np.ndarray
        First unit (N,M) vector of the coordinate system.
    e2 : np.ndarray
        Second unit (N,M) vector of the coordinate system.

    Notes
    -----
    e3 is orthogonal to the plane spanned by e1 and e2 following the right hand rule.
    Project e2 onto e1, and subtract to ensure orthogonality.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        Local orthonormal coordinate system ``e1, e2, e3``.
    """
    e1_norm = np.linalg.norm(e1, axis=1)
    bad_vectors = np.argwhere(e1_norm == 0).ravel()

    LOGGER.debug(f"{len(bad_vectors)} vectors have length zero.")

    e1_norm = np.where(e1_norm != 0, e1_norm, 1)
    e1 = e1 / e1_norm[:, None]

    # Ensure e1 and e2 are orthogonal
    dot_prod = np.einsum("ij,ij->i", e1, e2)
    e2 = e2 - dot_prod[:, None] * e1

    # Normalize
    e2 /= np.linalg.norm(e2, axis=1)[:, None]

    # Use right hand rule to compute e3
    e3 = np.cross(e1, e2)

    return e1, e2, e3


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
        Endocardium surface mesh. If provided, normal direction is updated by the
        surface normal instead of the Laplace solution.

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

    et, en, _ = orthogonalization(grid["grad_trans"], grid["k"])
    el = np.cross(en, et)

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
        Endocardium surface mesh. If provided, normal direction is updated by the
        surface normal instead of the Laplace solution.

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

    et, en, _ = orthogonalization(grid["grad_trans"], grid["k"])
    el = np.cross(en, et)

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


@deprecated(reason="Use _compute_rotation_angles instead.")
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
    # follow definition in Doste et al:
    # α = α_endo(w) · (1 − d) + α_epi(w) · d

    angle = np.zeros(grid.n_cells)
    # angle = rot_epi * (np.ones(grid.n_cells) - grid["d"]) + rot_endo * grid["d"]
    angle = rot_endo * (np.ones(grid.n_cells) - grid["d"]) + rot_epi * grid["d"]
    return angle


def _compute_rotation_angle(
    transmural_distance: float | list | np.ndarray,
    rotation_endocardium: float,
    rotation_epicardium: float,
):
    """Compute the rotation angle a for a given transmural depth and weight factor."""
    # follow definition in Doste et al:
    # α = α_endo(w) · (1 − d) + α_epi(w) · d
    angle = (
        rotation_endocardium * (1 - transmural_distance) + rotation_epicardium * transmural_distance
    )
    return angle


def compute_ventricle_fiber_by_drbm(
    directory: str,
    settings: dict = {
        "alpha_left": [60, -60],
        "alpha_right": [90, -25],
        "alpha_ot": None,
        "beta_left": [-65, 25],
        "beta_right": [-65, 25],
        "beta_ot": None,
    },
    left_only: bool = False,
) -> pv.UnstructuredGrid:
    """Compute the fiber coordinate system from the Laplace solutions.

    Parameters
    ----------
    directory : str
        Directory of d3plot/tprint files.
    settings : dict[str, list[float] | None]
        Rotation angles for fiber generation. The defaults are
        ``{"alpha_left": [60, -60], "alpha_right": [90, -25], "alpha_ot": None,
        "beta_left": [-65, 25], "beta_right": [-65, 25], "beta_ot": None}``.
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
        solutions.extend(["ab_r", "ot_r", "w_r"])

    data = read_laplace_solution(directory, field_list=solutions, read_heatflux=True)
    grid = data.point_data_to_cell_data()

    if left_only:
        # label to 1 for all cells
        left_mask = np.ones(grid.n_cells, dtype=bool)
        grid.cell_data["label"] = np.ones(grid.n_cells, dtype=int)
        right_mask = np.invert(left_mask)
    else:
        # label to 1 for left ventricle, 2 for right ventricle
        left_mask = grid["trans"] <= 0
        right_mask = grid["trans"] > 0
        grid.cell_data["label"] = np.zeros(grid.n_cells, dtype=int)
        grid.cell_data["label"][left_mask] = 1
        grid.cell_data["label"][right_mask] = 2

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

    # Build local coordinate system:
    # The right ventricle transmural gradient is flipped to ensure
    # a consistent coordinate system:
    # e_t points from endocardium to epicardium
    # e_n points from apex to base
    # e_c = e_n x e_t
    grid.cell_data["grad_trans"][right_mask] *= -1.0  # both LV & RV point to inside

    # Create orthonormal coordinate system
    en, et, ec = orthogonalization(k, grid["grad_trans"])

    # Add (unrotated) local coordinate system
    grid.cell_data["e_c"] = ec  # circumferential direction
    grid.cell_data["e_n"] = en  # normal/longitudinal direction
    grid.cell_data["e_t"] = et  # transmural direction

    if left_only:
        grid["d"] = grid["trans"]
    else:
        # normalize transmural distance to [0,1) in each ventricle
        # where 0 is endocardium, and 1 is epicardium
        d_l = np.absolute(grid["trans"][left_mask] / 2)
        d_r = np.absolute(grid["trans"][right_mask])
        grid["d"] = np.zeros(grid.n_cells, dtype=float)
        grid["d"][left_mask] = d_l
        grid["d"][right_mask] = d_r
        grid["d"] = grid["d"] * -1 + 1

    # rotation angles alpha and beta for each cell
    alpha = np.zeros(grid.n_cells)
    beta = np.zeros(grid.n_cells)

    alpha[left_mask] = _compute_rotation_angle(
        grid["d"][left_mask], settings["alpha_left"][0], settings["alpha_left"][1]
    )
    alpha[right_mask] = _compute_rotation_angle(
        grid["d"][right_mask], settings["alpha_right"][0], settings["alpha_right"][1]
    )

    beta[left_mask] = _compute_rotation_angle(
        grid["d"][left_mask], settings["beta_left"][0], settings["beta_left"][1]
    )
    beta[right_mask] = _compute_rotation_angle(
        grid["d"][right_mask], settings["beta_right"][0], settings["beta_right"][1]
    )

    # save rotation angles
    grid.cell_data["alpha"] = alpha
    grid.cell_data["beta"] = beta

    # 1) rotate vector ec counterclockwise around et by an angle alpha
    rot_alpha = R.from_rotvec(alpha[:, None] * et, degrees=True)

    fibers = rot_alpha.apply(ec)
    cross_fibers = rot_alpha.apply(en)
    sheets = et

    # 2) rotate vector ec counterclockwise around el or fibers by an angle beta
    rot_beta = R.from_rotvec(beta[:, None] * fibers, degrees=True)

    cross_fibers = rot_beta.apply(cross_fibers)
    sheets = rot_beta.apply(sheets)

    # NOTE Can add additional rotation in transverse direction, by specifying a
    # transverse angle gamma.

    # {f,n,s} in Piersanti et al. cross-fiber is sheet normal n
    # {F,T,S} in Bayer et al. cross-fiber is sheet normal S
    grid.cell_data["fiber"] = fibers
    grid.cell_data["cross-fiber"] = cross_fibers
    grid.cell_data["sheet"] = sheets

    grid.save("d-rbm-fibers.vtu")

    return grid.copy()
