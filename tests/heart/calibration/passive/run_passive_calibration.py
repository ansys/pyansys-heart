import os
import pathlib
import sys

import numpy as np

sys.path.append(r"D:\pyheart-lib")
from ansys.heart.calibration.passive_calibration import PassiveCalibration

if __name__ == "__main__":
    """
    Run.bat is the main input of LS-OPT, it actives pyheart-lib venv and then run this script.

    Patient specific parameters are defined in material_calibration.k and will be called by main.k
    """
    path_to_case = pathlib.Path(__file__).absolute().parents[0]

    os.chdir(path_to_case)

    case = PassiveCalibration(path_to_case)
    message = case.run_lsdyna()
    print(message)
    case.load_results()
    error = case.compute_error()
    np.savetxt("result", [error])
    fig = case.plot()
    fig.savefig("vs.png")
