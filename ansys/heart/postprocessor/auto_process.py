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

"""Script used to post process simulations automatically."""

import copy
import glob
import json
import os
import pathlib

heart_version = os.getenv("ANSYS_HEART_MODEL_VERSION")
from ansys.heart.core import LOG as LOGGER
from ansys.heart.postprocessor.Klotz_curve import EDPVR
from ansys.heart.postprocessor.SystemModelPost import SystemModelPost
from ansys.heart.postprocessor.aha17_strain import AhaStrainCalculator
from ansys.heart.postprocessor.dpf_utils import D3plotReader
from ansys.heart.postprocessor.exporter import LVContourExporter
from ansys.heart.preprocessor.mesh.objects import Cavity
from ansys.heart.simulator.settings.settings import SimulationSettings
import matplotlib.pyplot as plt
import numpy as np

if not heart_version:
    heart_version = "v0.1"
if heart_version == "v0.2":
    from ansys.heart.preprocessor.models.v0_2.models import HeartModel
elif heart_version == "v0.1":
    from ansys.heart.preprocessor.models.v0_2.models import HeartModel


def zerop_post(directory: str, model: HeartModel) -> tuple[dict, np.ndarray, np.ndarray]:
    """Post-process zeropressure folder.

    Parameters
    ----------
    directory : str
        simulation folder path
    model : HeartModel
        model

    Returns
    -------
    tuple[dict, np.ndarray, np.ndarray]
        dictionary with convergence information
        stress free configuration
        computed end-of-diastolic configuration
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
        LOGGER.error("Stress free computation is not converged, skip post process.")
        exit()

    stress_free_coord = data.get_initial_coordinates()
    displacements = data.get_displacement()
    guess_ed_coord = stress_free_coord + displacements[-1]

    if len(model.cap_centroids) == 0 or heart_version == "v0.2":
        nodes = model.mesh.nodes
    else:  # TODO remove after migrating to v0.2
        # a center node for each cap has been created, add them into create the cavity
        nodes = np.vstack((model.mesh.nodes, np.zeros((len(model.cap_centroids), 3))))
        for cap_center in model.cap_centroids:
            nodes[cap_center.node_id] = cap_center.xyz

    # convergence information
    dst = np.linalg.norm(guess_ed_coord - nodes, axis=1)
    error_mean = np.mean(dst)
    error_max = np.max(dst)

    # geometry information
    temp_mesh = copy.deepcopy(model.mesh)
    for key in temp_mesh.point_data.keys():
        temp_mesh.point_data.remove(key)
    for key in temp_mesh.cell_data.keys():
        temp_mesh.cell_data.remove(key)
    temp_mesh.save(os.path.join(directory, folder, "True_ED.vtk"))

    temp_mesh.points = stress_free_coord
    temp_mesh.save(os.path.join(directory, folder, "zerop.vtk"))
    temp_mesh.points = stress_free_coord + displacements[-1]
    temp_mesh.save(os.path.join(directory, folder, "Simu_ED.vtk"))

    def compute_cavity_volume(cavity: Cavity) -> list:
        """Extract cavity volume."""
        volumes = []
        for i, dsp in enumerate(displacements):
            new_cavity = copy.deepcopy(cavity)
            new_cavity.surface.points = stress_free_coord + dsp
            new_cavity.surface.save(os.path.join(directory, folder, f"{cavity.name}_{i}.stl"))
            volumes.append(new_cavity.volume)

        return volumes

    # cavity information
    lv_volumes = compute_cavity_volume(model.left_ventricle.cavity)
    true_lv_ed_volume = model.left_ventricle.cavity.volume

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

        true_rv_ed_volume = model.right_ventricle.cavity.volume
        rv_volumes = compute_cavity_volume(model.right_ventricle.cavity)

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

    return dct, stress_free_coord, guess_ed_coord


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
    fig = system.plot_pv_loop(show_ed=True)
    fig.savefig(os.path.join(directory, folder, "pv.png"))

    fig = system.plot_pressure_flow_volume(system.lv_system)
    fig.savefig(os.path.join(directory, folder, "lv.png"))

    if len(model.cavities) > 1:
        fig = system.plot_pressure_flow_volume(system.rv_system)
        fig.savefig(os.path.join(directory, folder, "rv.png"))

    #
    out_dir = os.path.join(directory, "post", "pv")
    os.makedirs(out_dir, exist_ok=True)
    time = D3plotReader(os.path.join(directory, "d3plot")).time / 1000  # to second
    for it, tt in enumerate(time):
        # assume heart beat once per 0.8s # TODO
        ef = system.get_ejection_fraction(t_start=time[-1] - 0.8, t_end=time[-1])
        fig = system.plot_pv_loop(t_start=0, t_end=tt, ef=ef, show_ed=False)
        fig.savefig(os.path.join(out_dir, "pv_{0:d}.png".format(it)))
        plt.close()
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
    out_dir = os.path.join(directory, "post", "lrc_strain")
    os.makedirs(out_dir, exist_ok=True)
    aha_strain = AhaStrainCalculator(model, d3plot_file=os.path.join(directory, "d3plot"))
    aha_strain.compute_aha_strain(out_dir, with_vtk=True)

    return
