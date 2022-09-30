"""Methods for reading LS-DYNA binout file."""

import numpy as np


def read_nodout(fname: str) -> list:
    """
    Read nodout ASCII file.

    Parameters
    ----------
    fname: LSDYNA ASCII nodout file

    Returns
    -------
    list of node coordinates
    """
    with open(fname) as ff:
        data = ff.readlines()

    nb = int(data[-1].split()[0])
    x = []
    time = []
    for il, line in enumerate(data):
        if "x-coor" in line:
            x.append(il)
        elif "time" in line:
            t = float(line.split(" ")[-2])
            time.append(t)

    time = np.array(time[1::2])

    result = []
    for i in range(len(x)):
        result.append(np.loadtxt(data, skiprows=x[i] + 1, max_rows=nb)[:, 10:13])
    return time, result


if __name__ == "__main__":
    fname = r"nodout"
    time, coords = read_nodout(fname)
    from compute_volume import get_cavity_volume2

    lv_cavity = np.loadtxt("left_ventricle.segment", delimiter=",", dtype=int)
    rv_cavity = np.loadtxt("right_ventricle.segment", delimiter=",", dtype=int)
    for coord in coords:
        lv_volume = get_cavity_volume2(coord, lv_cavity)
        rv_volume = get_cavity_volume2(coord, rv_cavity)
        print(lv_volume)
