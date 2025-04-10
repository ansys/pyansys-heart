# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
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

import hashlib
import os
import pathlib
import tarfile
import tempfile
import unittest.mock as mock

import pytest
import validators

from ansys.health.heart.utils.download import (
    _SHA256_TABLE,
    _format_download_urls,
    _infer_extraction_path_from_tar,
    _validate_hash_sha256,
    download_case_from_zenodo,
)


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
def test_download_case(database_name):
    """Test unpacking cases from different repositories."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        download_folder = tempdir

        with mock.patch("httpx.stream") as mock_download:
            save_path = download_case_from_zenodo(
                database_name, 1, download_folder, overwrite=True, validate_hash=False
            )

            if database_name == "Rodero2021":
                expected_save_path = os.path.join(download_folder, database_name, "01", "01.tar.gz")
            elif database_name == "Strocchi2020":
                expected_save_path = os.path.join(download_folder, database_name, "01.tar.gz")

            assert save_path == expected_save_path

            mock_download.assert_called_once()

    return


@pytest.mark.parametrize("tar_subpaths", (["01/01.case"], ["01.vtk"]))
def test_infer_extraction_path_from_tar(tar_subpaths):
    """Test unpacking a downloaded case."""
    # two configurations:
    # Strocchi2020 --> 01.tar.gz > 01/01.case
    # Rodero2021 --> 01.tar.gz > > 01.vtk
    with mock.patch("tarfile.open", return_value=mock.MagicMock(tarfile.TarFile)) as mock_taropen:
        mock_taropen.return_value.getnames.return_value = tar_subpaths
        with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
            path = _infer_extraction_path_from_tar(os.path.join(tempdir, "mytar.tar.gz"))

            expected_path = str(pathlib.Path(tempdir, tar_subpaths[0]))
            assert path == expected_path


def test_validate_hash_function_001():
    """Test hash validator."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        # create dummy data
        path_file1 = os.path.join(tempdir, "file1.txt")
        path_file2 = os.path.join(tempdir, "file2.txt")

        fid = open(path_file1, "w")
        fid.writelines("abc")
        fid.close()

        fid = open(path_file2, "w")
        fid.writelines("abcd\nblbla")
        fid.close()

        expected_hash = hashlib.sha256(open(path_file1, "rb").read()).hexdigest()

        # override original hash value with expected value.
        _SHA256_TABLE["Rodero2021"][1] = expected_hash

        assert _validate_hash_sha256(path_file1, "Rodero2021", casenumber=1), (
            "Expecting matching hash function"
        )

        # hashes of two different files should be different
        assert not _validate_hash_sha256(path_file2, "Rodero2021", casenumber=1), (
            "Expecting non-matching hash function"
        )

    return
