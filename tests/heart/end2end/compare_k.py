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


def compare(folder1, folder2):
    list_of_k_files = glob.glob(os.path.join(folder1, "*.k"))
    for file in list_of_k_files:
        file_ref = file
        file1 = os.path.join(folder2, os.path.basename(file_ref))
        # compare the contents of the two files
        contents = read_file(file1)
        contents_ref = read_file(file_ref)
        if contents != contents_ref:
            print("Found difference in: {0}.".format(file_ref))


if __name__ == "__main__":
    folder1 = r"D:\ansysdev\Rodero2021\01\meca_fh\writer_zerop0"
    folder2 = r"D:\ansysdev\Rodero2021\01\meca_fh\writer_zerop"

    print("*****************")
    compare(folder1, folder2)
