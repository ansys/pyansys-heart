"""Sample script to compare .k files generated in two different directories.

Notes
-----
Loops over the list of subfolders and compares the .k files. Prints out
the path to the .k file if the contents are not the same.
"""

import argparse
import glob
import os
import pathlib


def normalize_line_endings(text: str) -> str:
    # ignore empty line
    text.replace("\n\n", "\n")

    return text.replace("\r\n", "\n").replace("\r", "\n")


def read_file(file: pathlib.Path) -> str:
    with open(file, encoding="utf-8") as ref:
        return normalize_line_endings(ref.read())


def compare(folder1: str, folder2: str, ignore: list = []) -> bool:
    """Compare k-files in two folders."""
    list_of_k_files = glob.glob(os.path.join(folder1, "*.k"))

    if len(list_of_k_files) == 0:
        print("No k-files found in {0}".format(folder1))
        return

    ignored = []
    not_found = []
    succeed = []
    failed = []

    for file in list_of_k_files:
        file_ref = pathlib.Path(file)
        file1 = pathlib.Path(os.path.join(folder2, os.path.basename(file_ref)))

        if file1.name in ignore:
            # print("%s ignored" % file1.name)
            ignored += [file1.name]
            continue

        if not os.path.isfile(file1):
            # print("Failed to find match for {0}".format(file_ref))
            not_found += [file1.name]
            continue

        # compare the contents of the two files
        contents = read_file(file1)
        contents_ref = read_file(file_ref)
        if contents != contents_ref:
            # print("Difference in: {0}.".format(file_ref))
            failed += [file1.name]
        else:
            # print("No difference found in: {0}".format(file_ref.name))
            succeed += [file1.name]

    summary_str = (
        "Total compared:\t{0}\nignored:{1}\nnot found:{2}\nfailed:\t{3}\nsucceed:{4}".format(
            len(list_of_k_files),
            "\t{0:d}\t{1}".format(len(ignored), ignored),
            "\t{0:d}\t{1}".format(len(not_found), not_found),
            "\t{0:d}\t{1}".format(len(failed), failed),
            "\t{0:d}\t{1}".format(len(succeed), succeed),
        )
    )
    print("*** Summary ****")
    print(summary_str)
    print("****************")

    # if any are not found or do not match return False
    if len(not_found + failed) != 0:
        return False
    # otherwise return True
    else:
        return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="k-file compare tool")

    parser.add_argument(
        "--folder1", help="Sets first folder to compare against.", type=str, required=True
    )

    parser.add_argument(
        "--folder2", help="Sets second folder to compare against.", type=str, required=True
    )

    parser.add_argument(
        "--ignore", nargs="+", help="k-file names to ignore.", required=False, default=[]
    )

    args = parser.parse_args()

    print("****************")
    compare(args.folder1, args.folder2, args.ignore)
