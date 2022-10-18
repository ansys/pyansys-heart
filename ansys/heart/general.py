"""Module containing some general methods."""
import distutils.spawn
import os
import shutil
import subprocess

from ansys.heart.custom_logging import LOGGER

# --------------------------------------------------------------------------
NCPU = 4
# LSDYNAPATH = "/mnt/d/wsl/lsdyna_mpp/mppdyna_d_sse2_linux86_64_intelmmpi"
# OS = "wsl"

LSDYNAPATH = r"C:\Program Files\ANSYS Inc\v222\ansys\bin\winx64\lsdyna_sp.exe"
OS = "Win"


def check_lsdyna_exe():
    """Check if LSDYNA executable exists."""
    if OS == "Win":
        result = str(distutils.spawn.find_executable(LSDYNAPATH))
    elif OS == "wsl":
        p = subprocess.Popen(["wsl", "which", LSDYNAPATH], stdout=subprocess.PIPE)
        result = p.stdout.read().decode().strip()

    if result.strip() != LSDYNAPATH:
        LOGGER.error("Cannot find LSDYNA executable.")

    return


def clean_directory(directory: str):
    """Clean the directory by removing it and re-creating it."""
    if os.path.isdir(directory):
        shutil.rmtree(directory)
        os.makedirs(directory)
    else:
        os.makedirs(directory)

    return


def run_lsdyna(sim_file: str, lsdynapath: str = LSDYNAPATH, ncpu: int = NCPU, options=""):
    """Run lsdyna in wsl."""
    # extract_binout = False
    # os.chdir(pathlib.Path(sim_file).parent)

    check_lsdyna_exe()

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
