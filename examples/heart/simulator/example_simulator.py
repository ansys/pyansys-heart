"""Example to pre-process data from Strocchi2020 and Cristobal2021"""
import os
import pathlib

import ansys.heart.preprocessor.models as models
from ansys.heart.simulator.simulator import Simulator
from ansys.heart.simulator.support import run_preprocessor
import ansys.heart.writer.dynawriter as writers

if __name__ == "__main__":

    """Full Heart example

    1. Extracts simulation mesh
    2. Writes files for mechanics, zero-pressure, fiber generation, and purkinje

    Please change paths
    """

    path_to_case = os.path.join(
        pathlib.Path(__file__).parents[3], "downloads\\Strocchi2020\\01\\01.case"
    )
    workdir = os.path.join(pathlib.Path(path_to_case).parent, "BiVentricle")
    path_to_model = os.path.join(workdir, "heart_model.pickle")

    model = models.HeartModel.load_model(path_to_model)
    sim = Simulator(model)

    file_to_run = os.path.join(workdir, "purkinjegeneration", "main_left_ventricle.k")
    dynapath = "D:/Fortran/intelMPI/mppdyna_25AUG22"

    sim.run_lsdyna(sim_file=file_to_run, lsdynapath=dynapath, NCPU=1)
