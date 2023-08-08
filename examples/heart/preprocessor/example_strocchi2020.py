"""Example to process one of Strocchi 2020 cases."""
import os
import pathlib

from ansys.heart.misc.downloader import download_case, unpack_case
import ansys.heart.preprocessor.models_new as models
from ansys.heart.simulator.support import (
    get_input_geom_and_part_defintions_from_public_database,
    preprocess_model,
)

run_extraction = True
if run_extraction:
    # download case from remote repository
    case_num = 1  # patient number 1
    database = "Strocchi2020"

    download_folder: pathlib.Path = os.path.join(pathlib.Path(__file__).parents[3], "downloads")
    case_path: pathlib.Path = download_case(
        database=database, case_number=case_num, download_folder=download_folder, overwrite=False
    )
    mesh_path = os.path.join(
        pathlib.Path(case_path).parents[0], "%02d" % (case_num,), "%02d.case" % (case_num,)
    )

    if not os.path.isfile(mesh_path):
        unpack_case(case_path)

    input_geom, part_definitions = get_input_geom_and_part_defintions_from_public_database(
        mesh_path, model_type="BiVentricle", database="Strocchi2020"
    )

    info = models.ModelInfo(
        input_geom,
        scalar="surface-id",
        part_definitions=part_definitions,
        work_directory=(
            "D:\\development\\pyheart-lib\\pyheart-lib\\downloads"
            + "\\Strocchi2020\\01\\test_new_model"
        ),
    )

    model = preprocess_model(info, "BiVentricle", clean_workdir=False, use_wrapper=True)


import os
import pathlib

import ansys.heart.preprocessor.models_new as models
from ansys.heart.simulator.simulator import MechanicsSimulator

pickle_path = (
    "D:\\development\\pyheart-lib\\pyheart-lib\\downloads"
    + "\\Strocchi2020\\01\\test_new_model\\heart_model.pickle"
)
model = models.HeartModel.load_model(pickle_path)

lsdyna_path = (
    r"D:\development\dyna-versions\ls-dyna_smp_d_Dev_97584-g1b99fd817b_winx64_ifort190.exe"
)
simulator = MechanicsSimulator(
    model=model,
    lsdynapath=lsdyna_path,
    dynatype="smp",
    num_cpus=4,
    simulation_directory=os.path.join(model.info.workdir, "simulation-mechanics"),
)

# load default settings.
simulator.settings.load_defaults()
# compute the fiber orientation
simulator.compute_fibers()
# visualize computed fibers
simulator.model.plot_fibers(n_seed_points=2000)
# compute the stress free configuration
# simulator.compute_stress_free_configuration()
# do the main simulation
simulator.simulate()
print("done")

# from ansys.heart.simulator.simulator import MechanicsSimulator


# to do


# remove cells that are not associated with parts.
