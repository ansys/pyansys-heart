"""Sample script to compare .k files generated in two different directories.

Notes
-----
Loops over the list of subfolders and compares the .k files. Prints out
the path to the .k file if the contents are not the same.
"""

import glob
import os
import pathlib


def normalize_line_endings(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n")


def read_file(file: pathlib.Path) -> str:
    with open(file, encoding="utf-8") as ref:
        return normalize_line_endings(ref.read())


if __name__ == "__main__":
    basedir = "D:\\development\\pyheart-lib\\pyheart-lib\\downloads\\Strocchi2020\\01"
    # paths to the main directories for which to compare the .k files.
    dir1 = os.path.join(basedir, "BiVentriclePyVista")
    dir_ref = os.path.join(basedir, "BiVentricle_reference")

    subfolders = [
        "electrophysiology",
        "fibergeneration",
        "mechanics",
        "purkinjegeneration",
        "zeropressuremechanics",
    ]
    for subfolder in subfolders:
        dir_ref_sub = os.path.join(dir_ref, subfolder)
        print("Checking {0}...".format(subfolder))
        list_of_k_files = glob.glob(os.path.join(dir_ref_sub, "*.k"))
        for file in list_of_k_files:
            file_ref = file
            file1 = os.path.join(dir1, subfolder, os.path.basename(file_ref))
            # compare the contents of the two files
            contents = read_file(file1)
            contents_ref = read_file(file_ref)
            if contents != contents_ref:
                # print("*****************")
                print("Found difference in: {0}.".format(file_ref))
                # print("*****************")
