"""Module containing some general methods."""
import os
import shutil
import subprocess

NCPU = 4
LSDYNAPATH = "mppdyna_d_sse2_linux86_64_intelmmpi"

from ansys.heart.custom_logging import LOGGER


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
        f.write("echo start lsdyna...\n")
        f.write(f"mpirun -np {ncpu} {lsdynapath_wsl} i={sim_file_wsl} {options}\n")
        # if extract_binout:
        #     f.write("l2a iter3.binout0000")

    command = ["powershell", "-Command", "wsl", "-e", "bash", "-lic", "./run_lsdyna.sh"]

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
