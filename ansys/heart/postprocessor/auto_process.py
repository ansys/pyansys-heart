"""Script used to post process simulations automatically."""
import copy
import glob
import json
import os

from ansys.heart.postprocessor.D3plotPost import LVContourExporter
from ansys.heart.postprocessor.Klotz_curve import EDPVR
from ansys.heart.postprocessor.SystemModelPost import SystemModelPost
from ansys.heart.postprocessor.dpf_d3plot import D3plotReader
from ansys.heart.preprocessor.mesh.objects import Cavity, SurfaceMesh
import numpy as np


def zerop_post(directory, simulator):
    """
    Post-process zeropressure folder.

    Parameters
    ----------
    directory: where d3plot files locate
    simulator: Simulator

    Returns
    -------
    Dictionary of stress free simulation report.
    """
    folder = "post"
    os.makedirs(os.path.join(directory, folder), exist_ok=True)

    # read from d3plot
    data = D3plotReader(glob.glob(os.path.join(directory, "iter*.d3plot"))[-1])

    expected_time = simulator.settings.stress_free.analysis.end_time.to("millisecond").m
    if data.time[-1] != expected_time:
        Warning("Stress free computation is not converged, skip post process.")
        return None
    stress_free_coord = data.get_initial_coordinates()
    displacements = data.get_displacement()

    if len(simulator.model.cap_centroids) == 0:
        nodes = simulator.model.mesh.nodes
    else:
        # a center node for each cap has been created, add them into create the cavity
        nodes = np.vstack(
            (simulator.model.mesh.nodes, np.zeros((len(simulator.model.cap_centroids), 3)))
        )
        for cap_center in simulator.model.cap_centroids:
            nodes[cap_center.node_id] = cap_center.xyz

    # convergence information
    dst = np.linalg.norm(stress_free_coord + displacements[-1] - nodes, axis=1)
    error_mean = np.mean(dst)
    error_max = np.max(dst)

    # geometry information
    temp_mesh = copy.deepcopy(simulator.model.mesh)
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
            print(f"Cannot find {name}.segment")

        volumes = []
        for i, dsp in enumerate(displacements):
            cavity_surface = SurfaceMesh(name=name, triangles=faces, nodes=stress_free_coord + dsp)
            cavity_surface.save(os.path.join(directory, folder, f"{name}_{i}.stl"))
            volumes.append(Cavity(cavity_surface).volume)

        return volumes

    # cavity information
    lv_volumes = _post_cavity("left_ventricle")

    for cavity in simulator.model.cavities:
        if cavity.name.lower() == "left ventricle":
            true_lv_ed_volume = cavity.volume

    volume_error = [(lv_volumes[-1] - true_lv_ed_volume) / true_lv_ed_volume]

    # Klotz curve information
    # unit is mL and mmHg
    lv_pr_mmhg = (
        simulator.settings.mechanics.boundary_conditions.end_diastolic_cavity_pressure[
            "left_ventricle"
        ]
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
    if len(simulator.model.cavities) > 1:
        rv_volumes = _post_cavity("right_ventricle")
        for cavity in simulator.model.cavities:
            if cavity.name.lower() == "right ventricle":
                true_rv_ed_volume = cavity.volume

        volume_error.append((rv_volumes[-1] - true_rv_ed_volume) / true_rv_ed_volume)

        rv_pr_mmhg = (
            simulator.settings.mechanics.boundary_conditions.end_diastolic_cavity_pressure[
                "right_ventricle"
            ]
            .to("mmHg")
            .m
        )
        dct["Right ventricle EOD pressure (mmHg)"] = rv_pr_mmhg
        dct["True right ventricle volume"] = true_rv_ed_volume
        dct["Simulation Right ventricle volume"] = rv_volumes

    with open(os.path.join(directory, folder, "Post_report.json"), "w") as f:
        json.dump(dct, f)

    return dct


def mech_post(directory, model):
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

    return
