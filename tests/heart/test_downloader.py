import hashlib
import json
import os

from ansys.heart.misc.downloader import _format_download_urls, validate_hash_sha256
import pytest
import validators

from tests.heart.conftest import get_workdir


@pytest.mark.parametrize(
    "database_name",
    ["Strocchi2020", "Cristobal2021"],
)
def test_download_urls(database_name):
    """Test if URL still valid and exists."""
    all_download_urls = _format_download_urls()
    for case_num in all_download_urls[database_name]:
        url = all_download_urls[database_name][case_num]
        assert validators.url(url), "No valid URL for Case {0} of database {1}".format(
            case_num, database_name
        )


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
    hash_table = {"Cristobal2021": {1: hashlib.sha256(open(path_file1, "rb").read()).hexdigest()}}
    with open(path_hash_table, "w") as outfile:
        json.dump(hash_table, outfile, indent=4, ensure_ascii=True)

    # hashes of the same file should be the same
    assert validate_hash_sha256(
        path_file2, "Cristobal2021", casenumber=1, path_hash_table=path_hash_table
    ), "Expecting matching hash function"

    # hashes of two different files should be different
    assert not validate_hash_sha256(
        path_file3, "Cristobal2021", casenumber=1, path_hash_table=path_hash_table
    ), "Expecting non-matching hash function"

    # cleanup
    os.remove(path_file1)
    os.remove(path_file2)
    os.remove(path_file3)
    os.remove(path_hash_table)

    return
