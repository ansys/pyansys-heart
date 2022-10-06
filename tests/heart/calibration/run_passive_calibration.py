import os
import pathlib
import sys

import numpy as np

sys.path.append(r"D:\pyheart-lib")
from ansys.heart.calibration.passive_calibration import PassiveCalibration

if __name__ == "__main__":
    "test"
    path_to_case = os.path.join(
        pathlib.Path(__file__).parents[3], "tests", "heart", "calibration", "case"
    )
    os.chdir(path_to_case)

    case = PassiveCalibration(path_to_case)
    case.run_lsdyna()
    case.load_results()
    error = case.compute_error()
    np.savetxt("result", [error])
    fig = case.plot()
    fig.savefig("vs.png")
