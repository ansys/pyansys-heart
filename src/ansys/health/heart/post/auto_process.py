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

"""Script used to postprocess simulations automatically."""

import copy
import glob
import json
import os
from pathlib import Path

import natsort
import numpy as np

from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.models import HeartModel
from ansys.health.heart.objects import Cavity
from ansys.health.heart.post.dpf_utils import D3plotReader, D3plotToVTKExporter
from ansys.health.heart.post.klotz_curve import EDPVR
from ansys.health.heart.post.pvloop import generate_pvloop
from ansys.health.heart.post.strain_calculator import AhaStrainCalculator
from ansys.health.heart.settings.settings import SimulationSettings


def _process_zerop_iteration_file(directory: str, model: HeartModel, filename: str):
    """Process a single iteration file from the zero-pressure simulation."""
    folder = "post"
    os.makedirs(os.path.join(directory, folder), exist_ok=True)

    data = D3plotReader(filename)

    report_filename = (
        Path(filename).parent
        / Path(folder)
        / Path(filename).name.replace(".d3plot", ".report.json")
    )

    # get load from settings file
    setting = SimulationSettings()
    setting.load(os.path.join(directory, "simulation_settings.yml"))

    stress_free_coord = data.get_initial_coordinates()
    displacements = [data.get_displacement_at(time=t) for t in data.time]
    guess_ed_coord = stress_free_coord + displacements[-1]

    nodes = model.mesh.points

    # convergence information
    dst = np.linalg.norm(guess_ed_coord - nodes, axis=1)
    error_mean = np.mean(dst)
    error_max = np.max(dst)

    # geometry information
    temp_mesh = copy.deepcopy(model.mesh)
    temp_mesh.clear_data()
    temp_mesh.save(os.path.join(directory, folder, "True_ED.vtk"))

    temp_mesh.points = stress_free_coord
    temp_mesh.save(os.path.join(directory, folder, "zerop.vtk"))
    temp_mesh.points = stress_free_coord + displacements[-1]
    temp_mesh.save(os.path.join(directory, folder, "Simu_ED.vtk"))

    report = {
        "Simulation output time (ms)": data.time.tolist(),
        "Convergence": {
            "max_error (mm)": error_max,
            "mean_error (mm)": error_mean,
        },
    }

    # Retrieve imposed cavity pressures
    imposed_cavity_pressures = {
        f"{k}": v.to("mmHg").m
        for k, v in setting.mechanics.boundary_conditions.end_diastolic_cavity_pressure.items()
    }
    # Maps boundary condition names to model cavity names
    bc_map = {
        "Left ventricle cavity": "left_ventricle",
        "Right ventricle cavity": "right_ventricle",
        "Left atrium cavity": "left_atrial",
        "Right atrium cavity": "left_atrial",
    }

    volume_info = {}

    # extract cavity information
    for cavity in model.cavities:
        cavity: Cavity
        cavity.surface = model.mesh.get_surface(cavity.surface.id)
        cavity.surface.compute_normals(inplace=True)
        true_ed_volume = cavity.volume
        inflated_volumes = []
        for i, dsp in enumerate(displacements):
            #! This needs to be refactored.
            new_cavity = copy.deepcopy(cavity)
            new_cavity.surface.points = (
                stress_free_coord[new_cavity.surface.global_node_ids_triangles]
                + dsp[new_cavity.surface.global_node_ids_triangles]
            )
            new_cavity.surface.save(os.path.join(directory, folder, f"{cavity.name}_{i}.vtk"))
            inflated_volumes.append(new_cavity.volume)

        volume_info[cavity.name] = {
            "imposed cavity pressure (mmHg)": imposed_cavity_pressures[bc_map[cavity.name]],
            "true end diastolic volume (mm3)": true_ed_volume,
            "simulated volumes (mm3)": inflated_volumes,
            "volume error (%)": (true_ed_volume - inflated_volumes[-1]) / true_ed_volume * 100,
        }

    # save cavity volumes in JSON
    report["Cavity volumes"] = volume_info

    with open(report_filename, "w") as f:
        json.dump(report, f, indent=4)

    return report, stress_free_coord, guess_ed_coord


def _plot_klotz_curve(directory: str, report: dict) -> None:
    """Plot the Klotz curve from the report."""
    lv_pr_mmhg = report["Cavity volumes"]["Left ventricle cavity"]["imposed cavity pressure (mmHg)"]
    lv_volumes = report["Cavity volumes"]["Left ventricle cavity"]["simulated volumes (mm3)"]
    true_lv_ed_volume = report["Cavity volumes"]["Left ventricle cavity"][
        "true end diastolic volume (mm3)"
    ]

    # Klotz curve information
    time_array = np.array(report["Simulation output time (ms)"])
    klotz = EDPVR(true_lv_ed_volume / 1000, lv_pr_mmhg)
    sim_vol_ml = [v / 1000 for v in lv_volumes]
    sim_pr = lv_pr_mmhg * time_array / time_array[-1]
    fig = klotz.plot_EDPVR(simulation_data=[sim_vol_ml, sim_pr])
    fig.savefig(os.path.join(directory, "post", "klotz.png"))

    return


def zerop_post(directory: str, model: HeartModel) -> tuple[dict, np.ndarray, np.ndarray]:
    """Postprocess the zero-pressure folder.

    Parameters
    ----------
    directory : str
        Path to the simulation folder.
    model : HeartModel
        Model to postprocess.

    Returns
    -------
    tuple[dict, np.ndarray, np.ndarray]
        Dictionary with convergence information,
        stress free configuration, and
        computed end-of-diastolic configuration.
    """
    # Iterate over all iteration files in the directory and generate a report.
    iter_files = natsort.natsorted(glob.glob(os.path.join(directory, "iter*.d3plot")))
    for file in iter_files:
        report, stress_free_coord, guess_ed_coord = _process_zerop_iteration_file(
            directory, model, file
        )

    # Generate the Klotz curve plot.
    _plot_klotz_curve(directory, report)

    # For backward compatibility, return the last report.
    with open(os.path.join(directory, "post", "Post_report.json"), "w") as f:
        json.dump(report, f, indent=4)

    return report, stress_free_coord, guess_ed_coord


def mech_post(directory: str, model: HeartModel) -> None:
    """Postprocess the mechanical simulation folder.

    Parameters
    ----------
    directory : str
        Path to the d3plot folder.
    model : HeartModel
        Heart model.
    """
    last_cycle_duration = 800
    folder = "post"
    os.makedirs(os.path.join(directory, folder), exist_ok=True)

    # create PV loop of last cycle
    out_dir = os.path.join(directory, "post", "pv")
    os.makedirs(out_dir, exist_ok=True)
    f = os.path.join(directory, "binout0000")
    if os.path.exists(f):
        generate_pvloop(f, out_dir=out_dir, t_to_keep=last_cycle_duration)
    else:
        f = os.path.join(directory, "binout")
        if os.path.exists(f):
            generate_pvloop(f, out_dir=out_dir, t_to_keep=last_cycle_duration)
        else:
            LOGGER.warning("Neither 'binout0000' nor 'binout' exists in the directory.")

    # write vtk files of last cycle
    out_dir = os.path.join(directory, "post", "vtks")
    os.makedirs(out_dir, exist_ok=True)
    exporter = D3plotToVTKExporter(os.path.join(directory, "d3plot"), t_to_keep=last_cycle_duration)
    for i, t in enumerate(exporter.save_time):
        # NOTE: the returned pv_object seems corrupted, I suspect it's a bug of pyvista
        exporter.convert_to_pvgrid_at_t(time=t, fname=os.path.join(out_dir, f"heart_{i}.vtu"))

    # compute strain of last cycle
    out_dir = os.path.join(directory, "post", "lrc_strain")
    os.makedirs(out_dir, exist_ok=True)
    aha_strain = AhaStrainCalculator(model, d3plot_file=os.path.join(directory, "d3plot"))
    aha_strain.compute_aha_strain(out_dir, write_vtk=True, t_to_keep=last_cycle_duration)

    return
