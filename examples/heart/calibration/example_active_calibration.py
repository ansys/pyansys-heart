"""
todo
"""
import os
import subprocess

from ansys.heart.calibration.active_calibration import ActiveCalibration

if __name__ == "__main__":

    lsdyna_path = r"D:\my_path_to_ls_dyna\lsdyna_executable.exe"
    ncpu = 2
    dynatype = "smp"  # depend on your lsdyna_path

    # make sure pickle model contains right fiber information
    # i.e. fibergeneration has been run
    mdoel_path = r"path_to_model.pickle"

    # Create a folder that contains a LSOPT project for passive calibration
    target_folder = r"my_targe_folder"
    ActiveCalibration.create_calibration_project(
        target_folder, lsdyna_path, ncpu, dynatype, mdoel_path
    )

    # Run calibration
    os.chdir(target_folder)
    subprocess.run(["my_path_to_ls_opt", "ActiveCalibration.lsopt"])
