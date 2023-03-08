import glob as glob
import os
import pathlib

ROOT_FOLDER = os.path.join(pathlib.Path(__file__).parent)
from ansys.heart.misc.downloader import download_case, unpack_case

"""

Note
----
Note for VS Code/conda users: for the moment it seems that for proper pytest discovery in
VS Code's native testing framework you need to install pyfluent into the
base virtual environment.

"""


def pytest_sessionstart(session):
    print("Starting pytest session")
    workdir = get_workdir()
    pass


def get_assets_folder():
    return os.path.join(ROOT_FOLDER, "assets")


def download_asset(database: str = "Strocchi2020", casenumber: int = 1) -> pathlib.Path:
    """Download and unpack the requested asset if it is not yet available."""
    download_dir = os.path.join(get_assets_folder(), "cases")

    path_to_case1 = os.path.join(
        download_dir, database, f"{casenumber:02d}", f"{casenumber:02d}.case"
    )
    path_to_case2 = os.path.join(
        download_dir, database, f"{casenumber:02d}", f"{casenumber:02d}.vtk"
    )
    if os.path.isfile(path_to_case1):
        path_to_case = path_to_case1
    if os.path.isfile(path_to_case2):
        path_to_case = path_to_case2
    else:
        path_to_case = os.path.join(
            download_dir, database, f"{casenumber:02d}", f"{casenumber:02d}.vtk"
        )
        if database == "Strocchi2020":
            path_to_case = path_to_case.replace(".vtk", ".case")

    if not os.path.isfile(path_to_case):
        print("Downloading asset.")
        path_to_zip = download_case(database, casenumber, download_dir)
        unpack_case(path_to_zip)
        if database == "Strocchi2020":
            path_to_case = os.path.join(
                os.path.dirname(path_to_zip),
                path_to_zip.replace(".tar.gz", ""),
                path_to_zip.replace(".tar.gz", ".case"),
            )
        elif database == "Cristobal2021":
            path_to_case = path_to_zip.replace(".tar.gz", ".vtk")
    else:
        print("Asset already exists. Skip downloading.")
    return path_to_case


def get_workdir():
    return os.path.join(ROOT_FOLDER, "workdir_tests")


def create_directory(directory: str):
    """Creates directory"""
    print("Creating directory for tests")
    if not os.path.isdir(directory):
        os.makedirs(directory)
    return


def clean_workdir(directory: str):
    print("Cleaning working directory for tests")
    filelist = glob.glob(os.path.join(directory, "*"))
    print("Files to remove:")
    print(filelist)
    for file in filelist:
        print("Removing: %s" % file)
        os.remove(file)
    return


def normalize_line_endings(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n")


def read_file(file: pathlib.Path) -> str:
    with open(file, encoding="utf-8") as ref:
        return normalize_line_endings(ref.read())


def compare_string_with_file(output: str, reference_file: str) -> None:
    """compare the string in output, with the contents of reference_file
    normalize all line endinges to \\n
    """
    output = normalize_line_endings(output)
    ref_contents = read_file(reference_file)

    assert output == ref_contents


def clean_directory(directory: str):
    """Cleans the directory by removing it and re-creating it"""
    import shutil

    if os.path.isdir(directory):
        shutil.rmtree(directory)
        os.mkdir(directory)
    else:
        os.mkdir(directory)
    return


def remove_keys_from_dict(dictionary: dict, exclude_keys=[]):
    """Removes specific keys from the dictionary"""
    new_d = {k: dictionary[k] for k in set(list(dictionary.keys())) - set(exclude_keys)}
    return new_d
