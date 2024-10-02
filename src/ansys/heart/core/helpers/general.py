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

"""Module containing some general methods."""

import distutils.spawn
import os
import shutil
import subprocess

from ansys.heart.core import LOG as LOGGER
from ansys.heart.core.models import HeartModel

# --------------------------------------------------------------------------
NCPU = 4
LSDYNAPATH = r"D:\wsl\lsdyna_mpp\mppdyna_d_sse2_linux86_64_intelmmpi"
OS = "wsl"

# LSDYNAPATH = r"C:\Program Files\ANSYS Inc\v222\ansys\bin\winx64\lsdyna_dp.exe"
# OS = "Win"


def _check_lsdyna_exe():
    """Check if LSDYNA executable exists."""
    if OS == "Win":
        result = str(distutils.spawn.find_executable(LSDYNAPATH))
    elif OS == "wsl":
        p = subprocess.Popen(["wsl", "which", LSDYNAPATH], stdout=subprocess.PIPE)
        result = p.stdout.read().decode()

    if result.strip() != LSDYNAPATH:
        LOGGER.error("Cannot find LSDYNA executable.")
        exit()

    return


def _clean_directory(directory: str):
    """Clean the directory by removing it and re-creating it."""
    if os.path.isdir(directory):
        shutil.rmtree(directory)
        os.makedirs(directory)
    else:
        os.makedirs(directory)

    return


def _run_lsdyna(sim_file: str, lsdynapath: str = LSDYNAPATH, ncpu: int = NCPU, options="", OS=OS):
    """Run lsdyna in wsl."""
    # extract_binout = False
    # os.chdir(pathlib.Path(sim_file).parent)

    # check_lsdyna_exe()

    if OS == "wsl":
        sim_file_wsl = (
            subprocess.run(["wsl", "wslpath", os.path.basename(sim_file)], capture_output=1)
            .stdout.decode()
            .strip()
        )
        lsdynapath_wsl = (
            subprocess.run(["wsl", "wslpath", lsdynapath.replace("\\", "/")], capture_output=1)
            .stdout.decode()
            .strip()
        )

        with open("run_lsdyna.sh", "w", newline="\n") as f:
            f.write("#!/usr/bin/env sh\n")
            f.write("echo start lsdyna in wsl...\n")
            f.write(f"mpirun -np {ncpu} {lsdynapath_wsl} i={sim_file_wsl} {options}\n")
            # if extract_binout:
            #     f.write("l2a iter3.binout0000")

        command = ["powershell", "-Command", "wsl", "-e", "bash", "-lic", "./run_lsdyna.sh"]

    elif OS == "Win":
        command = [lsdynapath, f"ncpu={ncpu}", f"i={sim_file}", f"{options}"]

    LOGGER.info("Running ls-dyna with command:")
    run_command_display = " ".join([str(s) for s in command])
    LOGGER.info(run_command_display)

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    # stdout, stderr = process.communicate()
    with open("simulationlog.log", "w") as f:
        for line in process.stdout:
            f.write(line.decode())

    if process.returncode != 0:
        Exception("Simulation error.")

    LOGGER.info("Finished")

    return


def model_summary(model: HeartModel, attributes: list = None) -> dict:
    """Generate a dictionary with model information.

    Parameters
    ----------
    model : HeartModel
        HeartModel for which to generate the summary dictionary
    attributes : list
        List of attributes to try to add to the dict.

    Returns
    -------
    dict
        Dictionary with model information.
    """
    sum_dict = {}
    sum_dict["GENERAL"] = {}

    try:
        sum_dict["GENERAL"]["total_num_tets"] = model.mesh.tetrahedrons.shape[0]
        sum_dict["GENERAL"]["total_num_nodes"] = model.mesh.nodes.shape[0]
    except TypeError:
        LOGGER.info("Failed to format General model information.")

    sum_dict["PARTS"] = {}
    sum_dict["CAVITIES"] = {}
    for ii, part in enumerate(model.parts):
        sum_dict["PARTS"][part.name] = {}
        sum_dict["PARTS"][part.name]["num_tets"] = len(part.element_ids)

        sum_dict["PARTS"][part.name]["SURFACES"] = {}
        sum_dict["PARTS"][part.name]["CAPS"] = {}

        for surface in part.surfaces:
            sum_dict["PARTS"][part.name]["SURFACES"][surface.name] = {}
            sum_dict["PARTS"][part.name]["SURFACES"][surface.name]["num_faces"] = (
                surface.triangles.shape[0]
            )

            if attributes:
                for attribute in attributes:
                    try:
                        sum_dict["PARTS"][part.name]["SURFACES"][surface.name][attribute] = getattr(
                            surface.clean(), attribute
                        )
                    except AttributeError:
                        pass

        for cap in part.caps:
            sum_dict["PARTS"][part.name]["CAPS"][cap.name] = {}
            sum_dict["PARTS"][part.name]["CAPS"][cap.name]["num_nodes"] = len(
                cap.global_node_ids_edge
            )

            if attributes:
                for attribute in attributes:
                    try:
                        sum_dict["PARTS"][part.name]["CAPS"][cap.name][attribute] = getattr(
                            cap, attribute
                        )
                    except AttributeError:
                        pass

    for cavity in model.cavities:
        sum_dict["CAVITIES"][cavity.name] = {}
        sum_dict["CAVITIES"][cavity.name]["volume"] = cavity.surface.volume

        if attributes:
            for attribute in attributes:
                try:
                    sum_dict["CAVITIES"][cavity.name][attribute] = getattr(cavity, attribute)
                except AttributeError:
                    pass

    return sum_dict
