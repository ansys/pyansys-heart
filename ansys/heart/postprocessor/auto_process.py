"""Script used to post process simulations automatically."""
import copy
import glob
import json
import os
import pathlib
from typing import List

from ansys.heart import LOG as LOGGER
from ansys.heart.postprocessor.Klotz_curve import EDPVR
from ansys.heart.postprocessor.SystemModelPost import SystemModelPost
from ansys.heart.postprocessor.aha17_strain import AhaStrainCalculator
from ansys.heart.postprocessor.dpf_utils import D3plotReader
from ansys.heart.postprocessor.exporter import LVContourExporter
from ansys.heart.preprocessor.mesh.objects import Cavity, SurfaceMesh
from ansys.heart.simulator.settings.settings import SimulationSettings
import numpy as np
import pyvista as pv


def zerop_post(directory, model):
    """
    Post-process zeropressure folder.

    Parameters
    ----------
    directory: where d3plot files locate
    model: Heart model

    Returns
    -------
    Dictionary of stress free simulation report.
    """
    folder = "post"
    os.makedirs(os.path.join(directory, folder), exist_ok=True)

    # read from d3plot
    data = D3plotReader(glob.glob(os.path.join(directory, "iter*.d3plot"))[-1])
    # read settings file
    setting = SimulationSettings()
    setting.load(os.path.join(directory, "simulation_settings.yml"))

    expected_time = setting.stress_free.analysis.end_time.to("millisecond").m

    if data.time[-1] != expected_time:
        Warning("Stress free computation is not converged, skip post process.")
        return None
    stress_free_coord = data.get_initial_coordinates()
    displacements = data.get_displacement()

    if len(model.cap_centroids) == 0:
        nodes = model.mesh.nodes
    else:
        # a center node for each cap has been created, add them into create the cavity
        nodes = np.vstack((model.mesh.nodes, np.zeros((len(model.cap_centroids), 3))))
        for cap_center in model.cap_centroids:
            nodes[cap_center.node_id] = cap_center.xyz

    # convergence information
    dst = np.linalg.norm(stress_free_coord + displacements[-1] - nodes, axis=1)
    error_mean = np.mean(dst)
    error_max = np.max(dst)

    # geometry information
    temp_mesh = copy.deepcopy(model.mesh)
    for key in temp_mesh.point_data.keys():
        temp_mesh.point_data.remove(key)
    for key in temp_mesh.cell_data.keys():
        temp_mesh.cell_data.remove(key)
    temp_mesh.save(os.path.join(directory, folder, "True_ED.vtk"))

    # Note: vtk files contain center cap node, but not cap mesh
    temp_mesh.points = stress_free_coord
    temp_mesh.save(os.path.join(directory, folder, "zerop.vtk"))
    temp_mesh.points = stress_free_coord + displacements[-1]
    temp_mesh.save(os.path.join(directory, folder, "Simu_ED.vtk"))

    def _post_cavity(name: str):
        """Extract cavity volume."""
        try:
            faces = (
                np.loadtxt(os.path.join(directory, name + ".segment"), delimiter=",", dtype=int) - 1
            )
        except FileExistsError:
            LOGGER.warning(f"Cannot find {name}.segment")

        volumes = []
        for i, dsp in enumerate(displacements):
            cavity_surface = SurfaceMesh(name=name, triangles=faces, nodes=stress_free_coord + dsp)
            cavity_surface.save(os.path.join(directory, folder, f"{name}_{i}.stl"))
            volumes.append(Cavity(cavity_surface).volume)

        return volumes

    # cavity information
    lv_volumes = _post_cavity("left_ventricle")

    for cavity in model.cavities:
        if cavity.name.lower() == "left ventricle":
            true_lv_ed_volume = cavity.volume

    volume_error = [(lv_volumes[-1] - true_lv_ed_volume) / true_lv_ed_volume]

    # Klotz curve information
    # unit is mL and mmHg
    lv_pr_mmhg = (
        setting.mechanics.boundary_conditions.end_diastolic_cavity_pressure["left_ventricle"]
        .to("mmHg")
        .m
    )

    klotz = EDPVR(true_lv_ed_volume / 1000, lv_pr_mmhg)
    sim_vol_ml = [v / 1000 for v in lv_volumes]
    sim_pr = lv_pr_mmhg * data.time / data.time[-1]

    fig = klotz.plot_EDPVR(simulation_data=[sim_vol_ml, sim_pr])
    fig.savefig(os.path.join(directory, folder, "klotz.png"))

    dct = {
        "Simulation output time (ms)": data.time.tolist(),
        "Left ventricle EOD pressure (mmHg)": lv_pr_mmhg,
        "True left ventricle volume (mm3)": true_lv_ed_volume,
        "Simulation Left ventricle volume (mm3)": lv_volumes,
        "Convergence": {
            "max_error (mm)": error_max,
            "mean_error (mm)": error_mean,
            "relative volume error (100%)": volume_error,
        },
    }

    # right ventricle exist
    if len(model.cavities) > 1:
        rv_volumes = _post_cavity("right_ventricle")
        for cavity in model.cavities:
            if cavity.name.lower() == "right ventricle":
                true_rv_ed_volume = cavity.volume

        volume_error.append((rv_volumes[-1] - true_rv_ed_volume) / true_rv_ed_volume)

        rv_pr_mmhg = (
            setting.mechanics.boundary_conditions.end_diastolic_cavity_pressure["right_ventricle"]
            .to("mmHg")
            .m
        )
        dct["Right ventricle EOD pressure (mmHg)"] = rv_pr_mmhg
        dct["True right ventricle volume"] = true_rv_ed_volume
        dct["Simulation Right ventricle volume"] = rv_volumes

    with open(os.path.join(directory, folder, "Post_report.json"), "w") as f:
        json.dump(dct, f)

    return dct


def mech_post(directory: pathlib.Path, model):
    """
    Post-process Main mechanical simulation folder.

    Parameters
    ----------
    directory: where d3plot files locate
    model: HeartModel

    """
    folder = "post"
    os.makedirs(os.path.join(directory, folder), exist_ok=True)

    system = SystemModelPost(directory)
    fig = system.plot_pv_loop()
    fig.savefig(os.path.join(directory, folder, "pv.png"))

    fig = system.plot_pressure_flow_volume(system.lv_system)
    fig.savefig(os.path.join(directory, folder, "lv.png"))

    if len(model.cavities) > 1:
        fig = system.plot_pressure_flow_volume(system.rv_system)
        fig.savefig(os.path.join(directory, folder, "rv.png"))

    #
    out_dir = directory / "post" / "pv"
    os.makedirs(out_dir, exist_ok=True)
    time = D3plotReader(directory / "d3plot").time / 1000  # to second
    for it, tt in enumerate(time):
        # assume heart beat once per 1s
        ef = system.get_ejection_fraction(t_start=time[-1] - 1, t_end=time[-1])
        fig = system.plot_pv_loop(t_start=0, t_end=tt, ef=ef)
        fig.savefig(out_dir / "pv_{0:d}.png".format(it))
    # build video with command
    # ffmpeg -f image2 -i pv_%d.png output.mp4

    exporter = LVContourExporter(os.path.join(directory, "d3plot"), model)

    model.compute_left_ventricle_anatomy_axis(first_cut_short_axis=0.2)
    exporter.export_contour_to_vtk("l4cv", model.l4cv_axis)
    exporter.export_contour_to_vtk("l2cv", model.l2cv_axis)
    normal = model.short_axis["normal"]
    p_start = model.short_axis["center"]
    for ap in model.left_ventricle.apex_points:  # use next()?
        if ap.name == "apex epicardium":
            p_end = ap.xyz

    for icut in range(2):
        p_cut = p_start + (p_end - p_start) * icut / 2
        cutter = {"center": p_cut, "normal": normal}
        exporter.export_contour_to_vtk(f"shor_{icut}", cutter)

    exporter.export_lvls_to_vtk("lvls")

    #
    out_dir = directory / "post" / "lrc_strain"
    os.makedirs(out_dir, exist_ok=True)
    aha_strain = AhaStrainCalculator(model, d3plot_file=directory / "d3plot")
    aha_strain.compute_aha_strain(out_dir, with_vtk=True)

    return


def read_uvc(
    directory,
):
    """Read UVC from d3plot files."""
    data = D3plotReader(os.path.join(directory, "apico-basal.d3plot"))
    grid = data.model.metadata.meshed_region.grid

    t = data.model.results.temperature.on_last_time_freq.eval()[0].data
    grid["apico-basal"] = t[::3]
    data = D3plotReader(os.path.join(directory, "transmural.d3plot"))
    t = data.model.results.temperature.on_last_time_freq.eval()[0].data
    grid["transmural"] = t[::3]

    data = D3plotReader(os.path.join(directory, "rotational.d3plot"))
    t = data.model.results.temperature.on_last_time_freq.eval()[0].data
    grid["rotational"] = t[::3]

    grid.set_active_scalars("transmural")
    grid.save(os.path.join(directory, "uvc.vtk"))
    return grid


def orthogonalization(grid) -> pv.UnstructuredGrid:
    """Gram–Schmidt orthogonalization."""
    norm = np.linalg.norm(grid["grad_trans"], axis=1)
    bad_cells = np.argwhere(norm == 0).ravel()

    print(
        f"{len(bad_cells)} cells has null gradient in trans direction."
        f" They should only at valve regions and can be checked from vtk file."
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
    grid = data.model.metadata.meshed_region.grid

    for name in field_list:
        data = D3plotReader(os.path.join(directory, name + ".d3plot"))
        t = data.model.results.temperature.on_last_time_freq.eval()[0].data
        grid[name] = t[::3]
        # grid.set_active_scalars(name)

    # note vtk gradient method shows warning/error for some cells
    grid2 = grid.point_data_to_cell_data()
    for name in field_list:
        derivative = grid2.compute_derivative(scalars=name, preference="cell")
        res = derivative["gradient"]
        grid2["grad_" + name] = res

    grid2.save("gradient.vtk")

    return grid2


def compute_la_fiber_cs(directory) -> pv.UnstructuredGrid:
    """Compute left atrium fibers coordinate system."""

    def bundle_selection(grid) -> pv.UnstructuredGrid:
        """Left atrium bundle selection."""
        # grid = pv.read(os.path.join(directory, "res.vtk"))

        # bundle selection
        tau_mv = 0.65
        tau_lpv = 0.65
        tau_rpv = 0.1

        grid["k"] = np.zeros((grid.n_cells, 3))
        grid["bundle"] = np.zeros(grid.n_cells, dtype=int)

        mask_mv = grid["r"] >= tau_mv
        grid["k"][mask_mv] = grid["grad_r"][mask_mv]
        grid["bundle"][mask_mv] = 1

        mask = np.invert(mask_mv) & (grid["v"] > tau_lpv)
        grid["k"][mask] = grid["grad_v"][mask]
        grid["bundle"][mask] = 2

        mask = np.invert(mask_mv) & (grid["v"] < tau_rpv)
        grid["k"][mask] = grid["grad_v"][mask]
        grid["bundle"][mask] = 3

        mask = grid["bundle"] == 0
        grid["k"][mask] = grid["grad_ab"][mask]

        # grid.save(os.path.join(directory, "res2.vtk"))

        return grid

    grid = get_gradient(directory, field_list=["trans", "ab", "v", "r"])
    # TODO sometimes, pv object broken when pass directly
    grid = pv.read("gradient.vtk")

    grid = bundle_selection(grid)

    grid.save("la_fiber.vtk")
    grid = pv.read("la_fiber.vtk")

    grid = orthogonalization(grid)
    grid.save("la_fiber.vtk")

    return grid


def compute_ra_fiber_cs(directory) -> pv.UnstructuredGrid:
    """Compute right atrium fibers coordinate system."""

    def bundle_selection(grid) -> pv.UnstructuredGrid:
        """Left atrium bundle selection."""
        # Ideal RA geometry
        tao_tv = 0.9
        tao_raw = 0.55
        tao_ct_minus = -0.18
        tao_ct_plus = -0.1
        tao_icv = 0.9
        tao_scv = 0.1
        tao_ib = 0.35
        tao_ras = 0.135

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
        ras = 7
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
                                tag[i] = ras
                            else:
                                k[i] = w_grad[i]
                                tag[i] = ras
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
                            tag[i] = ras
                        else:
                            k[i] = ab_grad[i]
                            tag[i] = raw_ist_raa

        grid["k"] = k
        grid["bundle"] = tag.astype(int)

        return grid

    grid = get_gradient(directory, field_list=["trans", "ab", "v", "r", "w"])

    grid = pv.read("gradient.vtk")

    grid = bundle_selection(grid)

    grid.save("ra_fiber.vtk")
    grid = pv.read("ra_fiber.vtk")

    grid = orthogonalization(grid)
    grid.save("ra_fiber.vtk")

    return grid
