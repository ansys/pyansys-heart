"""

"""
#################################################################
MPI = True
LSDYNA = "mppdyna"
NCPU = 8
#################################################################
import os
import shutil
import subprocess
from compute_volume import update_system_json


def run_lsdyna(sim_file, option=""):
    """
    Parameters
    ----------
    sim_file: "main.k" in most time
    option: like 'case' for zerop

    Returns
    -------

    """
    if MPI:
        commands = ["mpirun", "-np", str(NCPU), LSDYNA, "i=" + sim_file, option]
    else:
        commands = [LSDYNA, "i=" + sim_file, "ncpu=" + str(NCPU)]

    p = subprocess.run(commands, stdout=subprocess.PIPE)

    # info = p.stdout.decode()
    # if " E r r o r   t e r m i n a t i o n" in info:
    #     logging.error("LSYDNA simulation failed at this iteration")
    #     exit()
    # elif "N o r m a l    t e r m i n a t i o n" in info:
    #     logging.info("LSYDNA simulation terminates")


def get_zerop_guess_file():
    """
    Find the result file after zerop
    Returns
    -------

    """
    # todo: check if converged
    guess_files = []
    for file in os.listdir("."):
        if file[-5:] == "guess":
            guess_files.append(file)

    return guess_files[-1]


def main():
    """
    Run Fiber generation + Zerop + Closed loop
    Returns
    -------
    """

    # Run fiber generation
    os.chdir("fibergeneration")
    run_lsdyna("main.k")
    os.chdir("..")

    # todo: works only with BV, extend to 4C case
    shutil.copy2(
        os.path.join("fibergeneration", "element_solid_ortho.k"),
        os.path.join("zeropressure", "solid_elements.k"),
    )
    shutil.copy2(
        os.path.join("fibergeneration", "element_solid_ortho.k"),
        os.path.join("mechanics", "solid_elements.k"),
    )

    # run zerop simulation
    os.chdir("zeropressure")
    run_lsdyna("main.k", "case")

    # get the result file and copy it to main directory
    guess_file = get_zerop_guess_file()
    os.chdir("..")
    shutil.copy2(
        os.path.join("mechanics", "nodes.k"), os.path.join("mechanics", "nodes_eod.k"),
    )
    shutil.copy2(
        os.path.join("zeropressure", guess_file), os.path.join("mechanics", "nodes.k"),
    )

    os.chdir("mechanics")

    # change unstressed volume in Json file
    update_system_json("nodes.k")

    # run closed loop simulation
    run_lsdyna("main.k")
    return


if __name__ == "__main__":
    os.chdir("")
    main()
