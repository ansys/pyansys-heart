import os.path
import pathlib

from ansys.heart.calibration.passive_calibration import create_calibration_folder

if __name__ == "__main__":
    """
    Run this script to create a folder for passive calibration

    Note:
    1. Suppose all k files related to zerop are here.
    2. Fiber orientation has been updated (solid_elements.k)

    How to:
    - Test run
       "D:\LSOPT\lsopt_2022R1_x64_win\LSOPT_2022R1\lsopt.exe" -b PassiveCalibration.lsopt
    - Run whole calibration without -b
    - You can clean the LSOPT project with --clean
    """
    path_to_case = pathlib.Path(__file__).absolute().parents[3]
    path_to_case = os.path.join(path_to_case, "tests", "heart", "calibration", "passive")

    create_calibration_folder(path_to_case)
