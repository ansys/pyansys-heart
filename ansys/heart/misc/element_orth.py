"""parse LSDYNA *ELEMENT_SOLID_ORTHO keywords.
"""
import numpy as np


def read_orth_element_kfile(fn):
    def get_number_of_elements(file):
        lines = open(file).readlines()
        n = 0
        for line in lines:
            if line[0] == "*":
                n += 1
        return int((len(lines) - n) / 4)

    def generate_specific_rows(file, row_indices):
        with open(file) as f:
            lines = f.readlines()
        return [lines[i] for i in row_indices]

    nele = get_number_of_elements(fn)

    # skip first 1 row and read every 4 row
    skip_row = 1  # because the first line is *ELEMENT_SOLID_ORTHO
    every_row = 4

    # element ID and part ID
    indx = np.linspace(0, nele - 1, num=nele, dtype=int) * every_row + skip_row
    data = generate_specific_rows(fn, indx)
    ids = np.loadtxt(data, dtype=int)[:, :]

    # element connectivity
    indx = np.linspace(0, nele - 1, num=nele, dtype=int) * every_row + skip_row + 1
    data = generate_specific_rows(fn, indx)
    connect = np.loadtxt(data, dtype=int)[:, :]

    # fiber
    indx = np.linspace(0, nele - 1, num=nele, dtype=int) * every_row + skip_row + 2
    data = generate_specific_rows(fn, indx)
    fib = np.loadtxt(data)
    # sheet
    indx = np.linspace(0, nele - 1, num=nele, dtype=int) * every_row + skip_row + 3
    data = generate_specific_rows(fn, indx)
    sheet = np.loadtxt(data)

    # sort by ids
    index = np.argsort(ids[:, 0])
    elem_ids = ids[index, 0]
    part_ids = ids[index, 1]
    connect = connect[index]
    fib = fib[index]
    sheet = sheet[index]

    return elem_ids, part_ids, connect, fib, sheet


def write_orth_element_kfile(fname, elem_orth):

    with open(fname, "w") as f:
        f.write("*KEYWORDS\n")
        f.write("*ELEMENT_SOLID_ORTHO\n")

        for id, pid, elems, a, d in elem_orth:
            f.write("{0:8d}{1:8d}\n".format(id, pid))
            f.write((8 * "{:8d}" + "\n").format(*elems))
            f.write((3 * "{:16f}" + "\n").format(*a))
            f.write((3 * "{:16f}" + "\n").format(*d))

        f.write("*END\n")


def modify_ids_orth_elements():
    """
    Part ID is different from FiberGeneration module to simulation modules
    This script is to change them
    Returns
    -------

    """
    elem_ids, part_ids, connect, fib, sheet = read_orth_element_kfile("solid_elements.k")
    # Septum is a part of LV
    # Change part ID 3 to ID 1
    part_ids = np.where(part_ids == 3, 1, part_ids)
    write_orth_element_kfile("solid_elements.k", zip(elem_ids, part_ids, connect, fib, sheet))


if __name__ == "__main__":
    pass
