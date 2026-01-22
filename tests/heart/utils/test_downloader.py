# Copyright (C) 2023 - 2026 ANSYS, Inc. and/or its affiliates.
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
from pathlib import Path
import tempfile
from unittest import mock

import pytest

from ansys.health.heart.utils.download import (
    _SHA256_TABLE,
    _ZENODO_RECORDS,
    _get_file_download_url,
    _infer_extraction_path_from_tar,
    _validate_hash_sha256,
    download_case_from_zenodo,
)


@pytest.mark.parametrize(
    "database_name",
    ["Strocchi2020", "Rodero2021"],
)
def test_get_file_download_url_from_api(database_name):
    """Test fetching download URLs from Zenodo API."""
    record_id = _ZENODO_RECORDS[database_name]["record_id"]
    filename = "01.tar.gz"

    # Mock the httpx.get call to return a valid API response
    mock_response = mock.MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {
        "files": [
            {
                "key": filename,
                "size": 12345678,
                "links": {"self": f"https://zenodo.org/api/files/bucket-id/{filename}"},
            }
        ]
    }

    with mock.patch("httpx.get", return_value=mock_response) as mock_get:
        url, size = _get_file_download_url(record_id, filename)

        # Verify the API was called correctly
        mock_get.assert_called_once()
        call_args = mock_get.call_args
        assert f"/api/records/{record_id}" in call_args[0][0]

        # Verify the returned values
        assert url == f"https://zenodo.org/api/files/bucket-id/{filename}"
        assert size == 12345678


@pytest.mark.parametrize(
    "database_name",
    ["Rodero2021", "Strocchi2020"],
)
def test_download_case(database_name):
    """Test downloading cases from different repositories using Zenodo API."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        download_folder = Path(tempdir)

        # Mock the API call to get file metadata
        mock_api_response = mock.MagicMock()
        mock_api_response.status_code = 200
        mock_api_response.json.return_value = {
            "files": [
                {
                    "key": "01.tar.gz",
                    "size": 1024,
                    "links": {"self": "https://zenodo.org/api/files/bucket-id/01.tar.gz"},
                }
            ]
        }

        # Mock the download stream
        mock_stream_response = mock.MagicMock()
        mock_stream_response.headers = {"Content-Length": "1024"}
        mock_stream_response.num_bytes_downloaded = 1024
        mock_stream_response.iter_bytes = mock.MagicMock(return_value=[b"test data"])
        mock_stream_response.__enter__ = mock.MagicMock(return_value=mock_stream_response)
        mock_stream_response.__exit__ = mock.MagicMock(return_value=None)
        mock_stream_response.raise_for_status = mock.MagicMock()

        with mock.patch("httpx.get", return_value=mock_api_response):
            with mock.patch("httpx.stream", return_value=mock_stream_response):
                save_path = download_case_from_zenodo(
                    database_name, 1, download_folder, overwrite=True, validate_hash=False
                )

                if database_name == "Rodero2021":
                    expected_save_path = download_folder / database_name / "01" / "01.tar.gz"
                elif database_name == "Strocchi2020":
                    expected_save_path = download_folder / database_name / "01.tar.gz"

                assert save_path == expected_save_path
                assert save_path.exists()

    return


@pytest.mark.parametrize("tar_subpaths", (["01/01.case"], ["01.vtk"]))
def test_infer_extraction_path_from_tar(tar_subpaths):
    """Test unpacking a downloaded case."""
    # two configurations:
    # Strocchi2020 --> 01.tar.gz > 01/01.case
    # Rodero2021 --> 01.tar.gz > > 01.vtk

    # Mock tarfile.open to support context manager and return a mock tarball
    mock_tarball = mock.MagicMock()
    mock_tarball.getnames.return_value = tar_subpaths
    mock_tarfile_open = mock.MagicMock(return_value=mock_tarball)
    mock_tarball.__enter__.return_value = mock_tarball
    mock_tarball.__exit__.return_value = None

    with mock.patch("tarfile.open", mock_tarfile_open):
        with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
            tar_path = Path(tempdir) / "mytar.tar.gz"
            path = _infer_extraction_path_from_tar(tar_path)

            expected_path = str(Path(tempdir) / tar_subpaths[0])
            assert path == expected_path


def test_validate_hash_function_001():
    """Test hash validator."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        # create dummy data
        path_file1 = Path(tempdir) / "file1.txt"
        path_file2 = Path(tempdir) / "file2.txt"

        with open(path_file1, "w") as fid:
            fid.write("abc")

        with open(path_file2, "w") as fid:
            fid.write("abcd\nblbla")

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


def test_get_file_download_url_file_not_found():
    """Test that FileNotFoundError is raised when file not found in API response."""
    record_id = "3890034"
    filename = "nonexistent.tar.gz"

    mock_response = mock.MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {
        "files": [
            {
                "key": "01.tar.gz",
                "size": 12345,
                "links": {"self": "https://zenodo.org/api/files/bucket-id/01.tar.gz"},
            }
        ]
    }

    with mock.patch("httpx.get", return_value=mock_response):
        with pytest.raises(FileNotFoundError, match="File 'nonexistent.tar.gz' not found"):
            _get_file_download_url(record_id, filename)


def test_download_case_invalid_database():
    """Test that invalid database raises DatabaseNotSupportedError."""
    from ansys.health.heart.exceptions import DatabaseNotSupportedError

    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        with pytest.raises(DatabaseNotSupportedError):
            download_case_from_zenodo("InvalidDB", 1, Path(tempdir))


def test_download_case_invalid_case_number():
    """Test that invalid case number raises ValueError."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        with pytest.raises(ValueError, match="Case number .* is invalid"):
            download_case_from_zenodo("Strocchi2020", 999, Path(tempdir))


def test_download_case_http_error_handling():
    """Test that HTTP errors are handled gracefully."""
    import httpx

    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        # Mock successful API call
        mock_api_response = mock.MagicMock()
        mock_api_response.status_code = 200
        mock_api_response.json.return_value = {
            "files": [
                {
                    "key": "01.tar.gz",
                    "size": 1024,
                    "links": {"self": "https://zenodo.org/api/files/bucket-id/01.tar.gz"},
                }
            ]
        }

        # Mock failed download stream with proper httpx exception
        mock_stream_response = mock.MagicMock()
        mock_stream_response.raise_for_status = mock.MagicMock(
            side_effect=httpx.HTTPStatusError(
                "500 Server Error", request=mock.MagicMock(), response=mock.MagicMock()
            )
        )
        mock_stream_response.response = mock.MagicMock()
        mock_stream_response.response.status_code = 500
        mock_stream_response.__enter__ = mock.MagicMock(return_value=mock_stream_response)
        mock_stream_response.__exit__ = mock.MagicMock(return_value=None)

        with mock.patch("httpx.get", return_value=mock_api_response):
            with mock.patch("httpx.stream", return_value=mock_stream_response):
                result = download_case_from_zenodo(
                    "Strocchi2020", 1, Path(tempdir), validate_hash=False
                )
                assert result is None
