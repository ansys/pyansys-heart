"""

"""
#################################################################
MPI = True
LSDYNA = "mppdyna"
NCPU = 6
#################################################################
import os
import shutil
import subprocess
from compute_volume import update_system_json
from pericardium_offset import add_pericardium_offset

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
    Run Zerop + Closed loop
    Returns
    -------
    """

    # run zerop simulation
    os.chdir("lsdyna_files_zeropressure")
    run_lsdyna("main.k", "case")

    # get the result file and copy it to main directory
    guess_file = get_zerop_guess_file()
    os.chdir("..")
    shutil.copy2(
        os.path.join("lsdyna_files", "nodes.k"),
        os.path.join("lsdyna_files", "nodes_eod.k"),
    )
    shutil.copy2(
        os.path.join("lsdyna_files_zeropressure", guess_file),
        os.path.join("lsdyna_files", "nodes.k"),
    )

    os.chdir("lsdyna_files")

    # change unstressed volume in Josn file
    update_system_json("nodes.k")

    # add pericardium springs offset
    add_pericardium_offset()

    # run closed loop simulation
    run_lsdyna("main.k")
    return


if __name__ == "__main__":
    os.chdir("")
    main()
