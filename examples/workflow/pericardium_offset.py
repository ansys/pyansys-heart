import os
import shutil

import numpy as np


def read_nodes_k(f):
    with open(f) as ff:
        lines = ff.readlines()
    data = []
    for line in lines:
        if line[0] != "*" and line[0] != "#":
            data.append(line)
    x0 = np.genfromtxt(data, delimiter=[8, 16, 16, 16])
    coord = x0[x0[:, 0].argsort()][:, 1:]  # sort by node ID
    return coord


def add_pericardium_offset():

    x0 = read_nodes_k("nodes.k")
    x1 = read_nodes_k("nodes_eod.k")

    shutil.copy2("pericardium.k", "pericardium_old.k")
    with open("pericardium_old.k") as f:
        peri_data = f.readlines()

    per_spring = np.genfromtxt(
        peri_data, skip_header=21, comments="*", delimiter=[8, 8, 8, 8, 8, 16]
    )
    per_nodes = per_spring[::3, 2].astype(int)

    per_x0 = x0[per_nodes - 1]
    per_x1 = x1[per_nodes - 1]
    per_offset = (per_x1 - per_x0).ravel()
    # np.savetxt(r"D:\pyheart-lib\examples\heart\workdir\x0.Csv",per_x0)
    # np.savetxt(r"D:\pyheart-lib\examples\heart\workdir\x1.Csv",per_x1)

    table = peri_data[21 : 21 + 3 * len(per_nodes)]
    for il, line in enumerate(table):
        table[il] = line[:56] + "{0:8d}".format(1) + "{0:16f}\n".format(per_offset[il])

    with open("pericardium.k", "w") as f:
        f.writelines(peri_data[0:21])
        f.writelines(table)
        f.write("*END\n")


if __name__ == "__main__":
    os.chdir(r"D:\wsl\user_load\test_offset")
    add_pericardium_offset()
