import hashlib
import json
import os

import pytest
import validators

from ansys.heart.misc.downloader import (
    _format_download_urls,
    download_case,
    unpack_case,
    validate_hash_sha256,
)
from tests.heart.conftest import get_workdir


@pytest.mark.parametrize(
    "database_name",
    ["Strocchi2020", "Rodero2021"],
)
def test_download_urls(database_name):
    """Test if URL still valid and exists."""
    all_download_urls = _format_download_urls()
    for case_num in all_download_urls[database_name]:
        url = all_download_urls[database_name][case_num]
        assert validators.url(url), "No valid URL for Case {0} of database {1}".format(
            case_num, database_name
        )


@pytest.mark.parametrize(
    "database_name",
    ["Rodero2021", "Strocchi2020"],
)
def test_download_and_unpack(database_name, tmp_path):
    """Test unpacking cases from different repositories."""
    import pathlib

    download_folder = os.path.join(tmp_path, "downloads")
    save_path = download_case(database_name, 1, download_folder, overwrite=True)

    # assert
    if database_name in ["Rodero2021", "Cristobal2021"]:
        expected_save_path = os.path.join(download_folder, database_name, "01", "01.tar.gz")
    elif database_name == "Strocchi2020":
        expected_save_path = os.path.join(download_folder, database_name, "01.tar.gz")

    assert os.path.isfile(save_path)

    assert save_path == expected_save_path

    unpack_case(save_path)

    expected_paths = [
        pathlib.Path(download_folder, database_name, "01", "01.vtk"),
        pathlib.Path(download_folder, database_name, "01", "01.case"),
    ]

    found = False
    for expected_path in expected_paths:
        if os.path.isfile(expected_path):
            found = True
            break
    assert found, ".vtk or .case file not in expected path."

    return


def test_validate_hash_function_001():
    """Test hash validator."""

    # create dummy data
    path_file1 = os.path.join(get_workdir(), "file1.txt")
    path_file2 = os.path.join(get_workdir(), "file2.txt")
    path_file3 = os.path.join(get_workdir(), "file3.txt")

    fid = open(path_file1, "w")
    fid.writelines("abc")
    fid.close()

    fid = open(path_file2, "w")
    fid.writelines("abc")
    fid.close()
    path_hash_table = "hash_table.json"

    fid = open(path_file3, "w")
    fid.writelines("abcd")
    fid.close()

    path_hash_table = "hash_table.json"

    # write dummy hash table
    hash_table = {"Rodero2021": {1: hashlib.sha256(open(path_file1, "rb").read()).hexdigest()}}
    with open(path_hash_table, "w") as outfile:
        json.dump(hash_table, outfile, indent=4, ensure_ascii=True)

    # hashes of the same file should be the same
    assert validate_hash_sha256(
        path_file2, "Rodero2021", casenumber=1, path_hash_table=path_hash_table
    ), "Expecting matching hash function"

    # hashes of two different files should be different
    assert not validate_hash_sha256(
        path_file3, "Rodero2021", casenumber=1, path_hash_table=path_hash_table
    ), "Expecting non-matching hash function"

    # cleanup
    os.remove(path_file1)
    os.remove(path_file2)
    os.remove(path_file3)
    os.remove(path_hash_table)

    return
