import glob as glob
import os
import pathlib

ROOT_FOLDER = os.path.join(pathlib.Path(__file__).parent)

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
