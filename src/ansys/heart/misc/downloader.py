# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Auto downloads cases.

Auto downloads cases from the remote repositories of Strocchi et al 2020,
and Rodero et al 2021."""

import hashlib

# from importlib.resources import files
from importlib.resources import path as resource_path
import json
import os
from pathlib import Path, PurePath
import typing

from tqdm import tqdm
import validators

from ansys.heart.core import LOG as LOGGER

try:
    import wget  # type: ignore
except ImportError:
    LOGGER.error("wget not installed but required. Please install by: pip install wget")
    exit()

URLS = {
    "Strocchi2020": {"url": "https://zenodo.org/record/3890034", "num_cases": 24},
    "Rodero2021": {"url": "https://zenodo.org/record/4590294", "num_cases": 20},
}
VALID_DATABASES = list(URLS.keys())
DOWNLOAD_DIR = PurePath.joinpath(Path(__file__).parents[3], "downloads")

PATH_TO_HASHTABLE = resource_path(
    "ansys.heart.misc", "remote_repo_hash_table_sha256.json"
).__enter__()


def _format_download_urls():
    """Format the URLS for all cases."""
    download_urls = {}
    for database_name in URLS.keys():
        download_urls[database_name] = {}
        url = URLS[database_name]["url"]
        num_cases = URLS[database_name]["num_cases"]
        for case_number in range(1, num_cases + 1):
            download_urls[database_name][case_number] = "{:}/files/{:02d}.tar.gz?download=1".format(
                url, case_number
            )
    return download_urls


ALL_DOWNLOAD_URLS = _format_download_urls()


def download_case_from_zenodo(
    database: str,
    case_number: int,
    download_folder: Path,
    overwrite: bool = True,
    validate_hash: bool = True,
) -> Path:
    """Download a case from the remote repository

    Parameters
    ----------
    database : str
        name of the database. Either Strocchi2020 or Rodero2021
    case_number : int
        case number to download
    download_folder : Path
        path to the folder in which to download the case

    Returns
    -------
    bool
        flag indicating whether download and unpacking was successful

    Examples
    --------
    Download case 1 from the public repository (Strocchi2020) of pathological hearts.
    >>> path_to_case = download_case(
            database="Strocchi2020", case_number=1, download_folder="my/download/folder"
        )

    Download case 1 from the public repository (Rodero2021) of 'healthy' hearts.
    >>> path_to_case = download_case(
        database="Rodero2021", case_number=1, download_folder="my/download/folder"
        )
    """

    if database not in VALID_DATABASES:
        raise ValueError("Database not valid, please specify valid database: %s" % VALID_DATABASES)

    url = URLS[database]["url"]
    if case_number > URLS[database]["num_cases"]:
        raise ValueError(
            "Database {0} only has {1} cases".format(database, URLS[database]["num_cases"])
        )

    if database == "Rodero2021":
        save_dir = os.path.join(download_folder, database, "{:>02d}".format(case_number))
    elif database == "Strocchi2020":
        save_dir = os.path.join(download_folder, database)

    save_path = os.path.join(save_dir, "{:02d}.tar.gz".format(case_number))

    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)

    if not overwrite and os.path.isfile(save_path):
        LOGGER.warning(f"File {save_path} already exists. Skipping...")
        return save_path

    download_url = "{:}/files/{:02d}.tar.gz?download=1".format(url, case_number)

    if download_url not in list(ALL_DOWNLOAD_URLS[database].values()):
        LOGGER.error(f"{download_url} not valid.")
        return None

    # validate URL
    if not validators.url(download_url):
        LOGGER.error(f"{download_url} not valid")
        return None

    try:
        wget.download(download_url, save_path)
    except:
        LOGGER.error(f"Failed to download from {download_url}")
        return None

    if validate_hash:
        is_valid_file = _validate_hash_sha256(
            file_path=save_path,
            database=database,
            casenumber=case_number,
            path_hash_table=PATH_TO_HASHTABLE,
        )
    else:
        LOGGER.warning("Not validating hash. Proceed at own risk")
        is_valid_file = True
    if not is_valid_file:
        LOGGER.error("File data integrity can not be validated.")
        os.remove(save_path)

    return save_path


def _validate_hash_sha256(
    file_path: Path, database: str, casenumber: int, path_hash_table: Path
) -> bool:
    """Check the file's hash function against the expected sha256 hash function."""
    if os.path.isfile(path_hash_table):
        fid = open(path_hash_table, "r")
        sha256_table: dict = json.load(fid)
        fid.close()
        # convert strings to ints
        for key in sha256_table.keys():
            sha256_table[key] = {int(k): v for k, v in sha256_table[key].items()}
    else:
        raise FileExistsError("File does not exist")

    try:
        sha256_table[database][casenumber]
    except KeyError:
        raise KeyError(
            "{0} : {1} is not yet present in the hash table dictionary".format(database, casenumber)
        )

    sha256 = hashlib.sha256(open(file_path, "rb").read()).hexdigest()
    if sha256 == sha256_table[database][casenumber]:
        return True
    else:
        return False


def unpack_case(tar_path: Path):
    """Untar the downloaded tar-ball."""
    import tarfile

    try:
        tar_ball = tarfile.open(tar_path)
        tar_dir = os.path.dirname(tar_path)
        tar_ball.extractall(path=tar_dir)
        return True
    except:
        LOGGER.error("Unpacking failed...")
        return False


def download_all_cases(download_dir: str = None):
    """Download all supported cases.

    Parameters
    ----------
    download_dir : str
        Base directory where to download the cases to.

    Notes
    -----
    Note that downloading all cases may - depending on bandwidth - take substantial
    time.

    """
    if download_dir == None:
        download_dir == DOWNLOAD_DIR

    if not os.path.isdir(download_dir):
        raise FileExistsError(f"{download_dir} does not exist.")

    tar_files = []
    for database_name, subdict in URLS.items():
        num_cases = subdict["num_cases"]
        download_dir = PurePath.joinpath(download_dir)
        for ii in range(1, num_cases + 1):
            LOGGER.info("Downloading {0} : {1}".format(database_name, ii))
            path_to_tar_file = download_case_from_zenodo(database_name, ii, download_dir)
            tar_files = tar_files + path_to_tar_file
    return tar_files


def unpack_cases(list_of_tar_files: typing.List):
    """Un-tar the downloaded cases."""
    for file in tqdm(list_of_tar_files):
        unpack_case(file)
    return


if __name__ == "__main__":
    LOGGER.info("Protected")
