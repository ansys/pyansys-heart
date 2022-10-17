import os
import pathlib

from ansys.heart.calibration.passive_calibration import PassiveCalibration
from ansys.heart.general import run_lsdyna

if __name__ == "__main__":
    """
    Run.bat is the main input of LS-OPT, it actives pyheart-lib venv and then run this script.
    Patient specific parameters are defined in material_calibration.k and will be called by main.k
    Do not forget compute fiber at first.

    How to run a passive calibration:
    1. Copy the folder.
    2. adapt venv directory (where is pyheart-lib) in run.bat
    3. Test LSOPT project with a baseline run:
       "D:\LSOPT\lsopt_2022R1_x64_win\LSOPT_2022R1\lsopt.exe" -b PassiveCalibration.lsopt
    4. Run the previous command without -b

    """
    path_to_case = pathlib.Path(__file__).absolute().parents[0]

    os.chdir(path_to_case)

    case = PassiveCalibration(path_to_case)
    ncpu = 4
    lsdynapath = "mppdyna_d_sse2_linux86_64_intelmmpi"
    sim_file = os.path.join(path_to_case, "main.k")
    options = "case"

    run_lsdyna(sim_file, lsdynapath, ncpu=ncpu, options="case")
    case.load_results()
    error = case.compute_error()
    with open("result", "a") as f:
        f.write("{0:10.5e}".format(error))
    fig = case.plot()
    fig.savefig("vs.png")
