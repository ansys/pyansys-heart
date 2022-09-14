import numpy as np


def read_binout(fname):
    with open(fname) as ff:
        data = ff.readlines()

    nb = int(data[-1].split()[0])
    x = []
    for il, line in enumerate(data):
        if "x-coor" in line:
            x.append(il)

    if len(x) != 3:
        print("Cannot read BINOUT correctly")
        exit()

    coord0 = np.loadtxt(data, skiprows=x[0] + 1, max_rows=nb)[:, 10:13]
    coord1 = np.loadtxt(data, skiprows=x[2] + 1, max_rows=nb)[:, 10:13]
    return coord0, coord1


def get_error(coord0, coord1):
    dst = np.linalg.norm(coord1 - coord0, axis=1)
    mean = np.mean(dst)
    max = np.max(dst)
    print(max)
    print(mean)
    return mean, max


if __name__ == "__main__":
    fname = r"nodout"
    xm, x_end1 = read_binout(fname)
    from compute_volume import get_cavity_volume2

    lv_cavity = np.loadtxt("cavity_left_ventricle.segment", delimiter=",", dtype=int)
    rv_cavity = np.loadtxt("cavity_right_ventricle.segment", delimiter=",", dtype=int)
    lv0 = get_cavity_volume2(xm, lv_cavity)
    rv0 = get_cavity_volume2(xm, rv_cavity)
    lv_ed = get_cavity_volume2(x_end1, lv_cavity)
    rv_ed = get_cavity_volume2(x_end1, rv_cavity)
    print(lv0)
    print(rv0)
    print(lv_ed)
    print(rv_ed)
