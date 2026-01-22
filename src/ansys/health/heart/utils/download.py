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

"""Module containing methods to download cases from public databases."""

import hashlib
from pathlib import Path
import tarfile

import httpx
import rich.progress

from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.exceptions import DatabaseNotSupportedError

_ZENODO_RECORDS = {
    "Strocchi2020": {"record_id": "3890034", "num_cases": 24},
    "Rodero2021": {"record_id": "4590294", "num_cases": 20},
}
_VALID_DATABASES = list(_ZENODO_RECORDS.keys())
_DOWNLOAD_DIR = Path.home() / "pyansys-heart" / "downloads"
_ZENODO_API_BASE = "https://zenodo.org/api"

# checksum table for Rodero 2021 et al and Strocchi 2020 et al.
_SHA256_TABLE = {
    "Rodero2021": {
        1: "a89f4c667d56e49e368b7bab8ba67e7eb42def0aa076a700ad82487e98a18aaf",
        2: "a9c2e24928770473630ad226c6f3759cd0fd4bc68241bf7f83508964b06756eb",
        3: "5c21bd551e365f1f3817bbd8413710364253de1ab027d46c7993318c7add28d6",
        4: "c7fc82b239d8e64e6a52b96c402a06ac64dc26098ccf364bcc39ce7ce9a7d44e",
        5: "1c24680bd2dedc1fed876dc4dd37052d0986235f0a8c0742a57d983c27db7b1f",
        6: "ed465189d3c35acd586b9ce3fd1b7bfbe9aa11202df14bfa253ab53a0cc27038",
        7: "d4a0b430694d21ecb0acfe48f9c5c7fd0f54e03e4f833a7ed0db5ad046e0c5a2",
        8: "7d0e737462107a5831214750c37b40dcc4add5a3e6ed9ec993cb642ae3514708",
        9: "4a227cf5e27680c1d4d32402e2a3b85f3d87c26881251e07467113629b8a508f",
        10: "f5c130dfa2d7ef2e0dbf5af3ed101ba494ecc0bf18b262fa87750eb5bdc234f2",
        11: "eb3dfabe4c30a91e1437c3ebea68a0a8eddb3dd9b04c55bf8acaf46074c6548e",
        12: "31bf4de87f7b49c9193dc151668290b05382e5a62ecaa1a090eba3d9fd6fe454",
        13: "62b32451464bb173aa6616169e16c0a1caeca2a395ccfaad75af9e492e529ce6",
        14: "58302a951a46f8d8d2e8e13a93d9b2717f36ec8eb1a68562de8e74fd839e9322",
        15: "59c8e095a3defadc7ec25cb84ad0de344beb55b262edc57c1c2327c714b6e977",
        16: "85ca99c11f0fb8ff38f9d59dfaa0cd1182827ac587edb66d73b331e9fc52058e",
        17: "a6aac90fbdddf534f0f009c0b27b2c505ee63e7364673b4366d48b0123d9ff32",
        18: "b50919a711dc3914cc0a4dd9cabc19679a9f36be9ac6eb26a23ca6799f501720",
        19: "e763e5683da4a48bd8944f794f6cbccd3949109c713e258ebfae53033ffda043",
        20: "348526bb6fcb91b2258a46b05dd49e87cc759ea7fdf05049e798bdbe091cc523",
    },
    "Strocchi2020": {
        1: "2cd411b11e2787416dcd76a4fe99699efe35492bb6a38bcd096432c882507aac",
        2: "114f1194c59dbbbe358af514b0bb6594d0bb90df2e338f8e2147f22fd5030af0",
        3: "c0d871e9f72e5ba044cfd201391952e7f13580b44cff1f83041ac1d7a761f182",
        4: "f2a84273a418a9af4d4915e0114b1eba455abd7cc194f0ecc643f553bd5d3ce4",
        5: "037635e41cb10a4182523139e2b8d94ea58fcf41c7ee475d50105de8ca502b38",
        6: "a06917b7615c4bbb8c69202db624e02901f1e8cbb0982d027578abe1739482eb",
        7: "8493d0f3c71421323c09f622343e56438f61c4c7dfebdc50f842bd332c4d2182",
        8: "918ed5400cdb58ccc903edd3983f822f02bf3fdd218d8f79884ca59a6661853d",
        9: "25767eddca70fd1fa1d06e6cf8ea8b4791ef7da35c56b026c5c586ec43b84cb1",
        10: "9639e508b78fce99068cd37ce5c5b47e2e814f9922ecbb2d6859ab70020dc00c",
        11: "9fb1ac1c5679f1445e048dd74b9b55b765a3e05b672fdafbfe7a140a61edb4f1",
        12: "28da3ead2c65860ee681980f26d3232da75d510ff07338744574b3b6b1991f83",
        13: "9420d52ef9740193e1d642517c3f0779c96a4a3cd132a9de3480d46e7a54eb7a",
        14: "2bec0aa0412c8c4b0a4148bafeeba9d3124b6c914d80e3827d47c2d8cdd7e271",
        15: "d6c3fba6430b36b2c88a26f6fe0a21aeb14567f222b18cad268a7c8246ac6d14",
        16: "06334b08286ab2c56dfc34408d0f4b34042aa96eb55be43e8faebaad909a7c5a",
        17: "46050e2ccccaf63c1f036b2a1b6f8f3b5bf2940cc6b5d4f7087b204bf0a66859",
        18: "83a5e5b1c28bb855ef5327ec9e15f7947df42d410f334a1d5b31286ce0caad75",
        19: "1c301f8ea2ee0ded0d810fbb7da4e8065d2608ddacfa1ff0bcf21dc642645c61",
        20: "d97e19ac311492238797a3b3b7bb34f5542fed0a1cae4f6c776f216eb96f566e",
        21: "fe9e06a7c1e192b39a1bf6e7ae4218d96283164ea9217a39d6f9206dbfc7c413",
        22: "8ad46acd187875425ac85b12222767c2172a29ac30220736a461c077521c5a17",
        23: "7accf3273939ace04dc89ffbf41bb408785d6530c87b56a31e0c3badedbbae8f",
        24: "6a4cb391e38c6a63c75c324892f66924e4d02a29cb1313bce61a8e057949139a",
    },
}


def _get_file_download_url(record_id: str, filename: str) -> tuple[str, int, str]:
    """Get download URL, size, and checksum for a specific file from Zenodo API.

    Parameters
    ----------
    record_id : str
        Zenodo record identifier.
    filename : str
        Name of file to download (e.g., "01.tar.gz").

    Returns
    -------
    tuple[str, int, str]
        Download URL, file size in bytes, and checksum string (format: "algorithm:hash").

    Raises
    ------
    httpx.HTTPStatusError
        If API request fails.
    FileNotFoundError
        If file not found in record.
    ValueError
        If file metadata is incomplete.

    Examples
    --------
    >>> url, size, checksum = _get_file_download_url("3890034", "01.tar.gz")
    >>> print(checksum)  # e.g., "md5:2942bfabb3d05332b66eb128e0842cff"
    """
    api_url = f"{_ZENODO_API_BASE}/records/{record_id}"

    try:
        response = httpx.get(api_url, timeout=10.0, follow_redirects=True)
        response.raise_for_status()
        metadata = response.json()
    except httpx.HTTPStatusError as e:
        error_msg = f"Failed to fetch Zenodo record {record_id}: {e.response.status_code}"
        LOGGER.error(error_msg)
        raise
    except httpx.RequestError as e:
        error_msg = f"Network error fetching Zenodo record {record_id}: {e}"
        LOGGER.error(error_msg)
        raise

    try:
        # Find the requested file in the files array
        file_info = next(
            (f for f in metadata.get("files", []) if f["key"] == filename),
            None,
        )

        if file_info is None:
            available_files = [f["key"] for f in metadata.get("files", [])]
            raise FileNotFoundError(
                f"File '{filename}' not found in record {record_id}. "
                f"Available files: {available_files}"
            )

        download_url = file_info["links"]["self"]
        file_size = file_info["size"]
        checksum = file_info["checksum"]  # Format: "md5:hash_value" or similar

        return download_url, file_size, checksum

    except KeyError as e:
        error_msg = f"Invalid metadata structure for record {record_id}: missing {e}"
        LOGGER.error(error_msg)
        raise ValueError(error_msg) from e


def download_case_from_zenodo(
    database: str,
    case_number: int,
    download_folder: Path,
    overwrite: bool = True,
    validate_hash: bool = True,
) -> Path | None:
    """Download a case from the remote repository using Zenodo REST API.

    Parameters
    ----------
    database : str
        name of the database. Options are ``'Strocchi2020'`` or ``'Rodero2021'``.
    case_number : int
        Case number to download.
    download_folder : Path
        Path to the folder to download the case to.
    overwrite : bool, default: True
        Whether to overwrite existing files.
    validate_hash : bool, default: True
        Whether to validate file integrity using checksums.

    Returns
    -------
    Path or None
        Path to the tarball that contains the VTK/CASE files, or None if download failed.

    Raises
    ------
    DatabaseNotSupportedError
        If database name is not recognized.
    ValueError
        If case number is invalid for the database.

    Examples
    --------
    Download case 1 from the public repository (``'Strocchi2020'``) of pathological hearts.

    >>> path_to_tar_file = download_case_from_zenodo(
    ...     database="Strocchi2020", case_number=1, download_folder="my/download/folder"
    ... )

    Download case 1 from the public repository (``'Rodero2021'``) of healthy hearts.

    >>> path_to_tar_file = download_case_from_zenodo(
    ...     database="Rodero2021", case_number=1, download_folder="my/download/folder"
    ... )
    """
    # Validate database
    if database not in _VALID_DATABASES:
        raise DatabaseNotSupportedError(database, f"Please choose one of {_VALID_DATABASES}")

    # Validate case number
    num_cases = _ZENODO_RECORDS[database]["num_cases"]
    if not 1 <= case_number <= num_cases:
        raise ValueError(
            f"Case number {case_number} is invalid for {database}. Valid range: 1-{num_cases}"
        )

    # Setup paths using pathlib
    download_folder = Path(download_folder)
    if database == "Rodero2021":
        save_dir = download_folder / database / f"{case_number:02d}"
    else:  # Strocchi2020
        save_dir = download_folder / database

    save_path = save_dir / f"{case_number:02d}.tar.gz"

    # Create directory
    save_dir.mkdir(parents=True, exist_ok=True)

    # Check if file exists
    if not overwrite and save_path.exists():
        LOGGER.warning(f"File {save_path} already exists. Skipping...")
        return save_path

    # Get download URL from Zenodo API
    try:
        record_id = _ZENODO_RECORDS[database]["record_id"]
        filename = f"{case_number:02d}.tar.gz"
        download_url, expected_size, api_checksum = _get_file_download_url(record_id, filename)
        LOGGER.info(f"Downloading {database} case {case_number} from Zenodo API...")
        LOGGER.debug(f"Download URL: {download_url}")
        LOGGER.debug(f"Expected size: {expected_size / 1024 / 1024:.2f} MB")
        LOGGER.debug(f"API checksum: {api_checksum}")

    except (FileNotFoundError, ValueError, httpx.HTTPStatusError) as e:
        LOGGER.error(f"Failed to get download URL: {e}")
        return None

    # Download file with progress bar
    try:
        with httpx.stream(
            "GET",
            download_url,
            follow_redirects=True,
            timeout=300.0,  # 5 minute timeout
        ) as response:
            response.raise_for_status()

            total = int(response.headers.get("Content-Length", expected_size))

            with rich.progress.Progress(
                "[progress.percentage]{task.percentage:>3.0f}%",
                rich.progress.BarColumn(bar_width=None),
                rich.progress.DownloadColumn(),
                rich.progress.TransferSpeedColumn(),
            ) as progress:
                download_task = progress.add_task("Download", total=total)
                with open(save_path, "wb") as fp:
                    for chunk in response.iter_bytes():
                        fp.write(chunk)
                        progress.update(
                            download_task,
                            completed=response.num_bytes_downloaded,
                        )

        LOGGER.info(f"Download completed: {save_path}")

    except httpx.HTTPStatusError as e:
        error_msg = f"HTTP error downloading from {download_url}: {e.response.status_code}"
        LOGGER.error(error_msg)
        if save_path.exists():
            save_path.unlink()
        return None

    except httpx.RequestError as e:
        error_msg = f"Network error downloading from {download_url}: {e}"
        LOGGER.error(error_msg)
        if save_path.exists():
            save_path.unlink()
        return None

    except OSError as e:
        error_msg = f"File system error writing to {save_path}: {e}"
        LOGGER.error(error_msg)
        if save_path.exists():
            save_path.unlink()
        return None

    # Validate checksum
    if validate_hash:
        try:
            LOGGER.debug(f"Validating checksum: {api_checksum}")
            is_valid_file = _validate_checksum(save_path, api_checksum)

            if not is_valid_file:
                LOGGER.error(
                    f"Checksum validation failed for {save_path}. "
                    f"File may be corrupted. Expected: {api_checksum}"
                )
                save_path.unlink()
                return None

            LOGGER.info("Checksum validation passed.")

        except Exception as e:
            LOGGER.error(f"Error validating checksum: {e}")
            save_path.unlink()
            return None
    else:
        LOGGER.warning(
            "Skipping checksum validation. Downloaded file integrity cannot be guaranteed."
        )

    return save_path


def _validate_checksum(file_path: Path, expected_checksum: str) -> bool:
    """Validate file checksum against expected value from Zenodo API.

    Parameters
    ----------
    file_path : Path
        Path to the file to validate.
    expected_checksum : str
        Expected checksum in format "algorithm:hash" (e.g., "md5:abc123...").

    Returns
    -------
    bool
        True if checksum matches, False otherwise.

    Raises
    ------
    ValueError
        If checksum format is invalid or algorithm is unsupported.

    Examples
    --------
    >>> is_valid = _validate_checksum(Path("file.tar.gz"), "md5:2942bfabb3d05332b66eb128e0842cff")
    """
    try:
        # Parse checksum format: "algorithm:hash"
        if ":" in expected_checksum:
            algorithm, expected_hash = expected_checksum.split(":", 1)
        else:
            # Default to MD5 if no algorithm specified (Zenodo default)
            algorithm = "md5"
            expected_hash = expected_checksum

        # Normalize algorithm name
        algorithm = algorithm.lower()

        # Get the appropriate hash function
        if algorithm == "md5":
            hash_func = hashlib.md5()
        elif algorithm in ("sha256", "sha-256"):
            hash_func = hashlib.sha256()
        elif algorithm in ("sha512", "sha-512"):
            hash_func = hashlib.sha512()
        elif algorithm == "sha1":
            hash_func = hashlib.sha1()
        else:
            raise ValueError(f"Unsupported checksum algorithm: {algorithm}")

        # Calculate file checksum
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                hash_func.update(chunk)

        computed_hash = hash_func.hexdigest()

        # Compare checksums (case-insensitive)
        return computed_hash.lower() == expected_hash.lower()

    except Exception as e:
        LOGGER.error(f"Error validating checksum: {e}")
        raise


def _validate_hash_sha256(file_path: Path, database: str, casenumber: int) -> bool:
    """Check the file's hash function against the expected sha256 hash function."""
    try:
        _SHA256_TABLE[database][casenumber]
    except KeyError:
        raise KeyError(
            "{0} : {1} is not yet present in the hash table dictionary".format(database, casenumber)
        )

    sha256 = hashlib.sha256(open(file_path, "rb").read()).hexdigest()
    if sha256 == _SHA256_TABLE[database][casenumber]:
        return True
    else:
        return False


def _infer_extraction_path_from_tar(tar_path: str | Path) -> str:
    """Infer the path to the relevant CASE or VTK file from the tarball path."""
    tar_path = Path(tar_path)
    with tarfile.open(tar_path, "r:gz") as tarball:
        names = tarball.getnames()

    # Order matters: check if .case file exists before .vtk file
    sub_path = next((name for name in names if name.endswith(".case")), None)
    if not sub_path:
        sub_path = next((name for name in names if name.endswith(".vtk")), None)

    if sub_path is None:
        LOGGER.error(f"No relevant files are found in {tar_path}.")
        return str(tar_path)

    path = (tar_path.parent / sub_path).resolve()
    return str(path)


def _get_members_to_unpack(tar_ball: tarfile.TarFile) -> list:
    """Get the members to unpack from the tarball.

    Notes
    -----
    This ignores the large VTK files for the Strocchi 2020 archives.
    """
    if len(tar_ball.getnames()) > 1:
        members_to_unpack = [
            member for member in tar_ball.getmembers() if not member.name.endswith(".vtk")
        ]
    else:
        members_to_unpack = tar_ball.getmembers()
    return members_to_unpack


def _is_safe_tar_member(member: tarfile.TarInfo, target_dir: str):
    """Get safe members, prevent absolute paths and path traversal."""
    member_path = Path(target_dir) / member.name
    abs_target_dir = Path(target_dir).resolve()
    abs_member_path = member_path.resolve()
    return abs_member_path.is_relative_to(abs_target_dir)


def unpack_case(tar_path: Path, reduce_size: bool = True) -> str | bool:
    r"""Unpack the downloaded tarball file.

    Parameters
    ----------
    tar_path : Path
        Path to TAR.GZ file.
    reduce_size : bool, default: True
        Whether to reduce the size of the unpacked files by removing the VTK file for the
        Strocchi database.

    Examples
    --------
    >>> from ansys.health.heart.utils.download import unpack_case
    >>> path = unpack_case("Rodero2021\\01.tar.gz")

    Returns
    -------
    str
        Path to the CASE or VTK file.
    """
    try:
        tar_path = Path(tar_path)
        with tarfile.open(tar_path, "r:gz") as tar_ball:
            tar_dir = tar_path.parent
            if reduce_size:
                members = _get_members_to_unpack(tar_ball)
            else:
                members = tar_ball.getmembers()

            # Validate members
            unsafe_members = [m for m in members if not _is_safe_tar_member(m, str(tar_dir))]
            if unsafe_members:
                names = [m.name for m in unsafe_members]
                raise ValueError(f"Unsafe tar members detected in '{tar_path}': {names}")
            # Ignore bandit warnings, because members are validated
            tar_ball.extractall(path=tar_dir, members=members)  # nosec: B202

            path = _infer_extraction_path_from_tar(tar_path)
            return path

    except Exception as exception:
        LOGGER.error(f"Unpacking failed. {exception}")
        return False


def download_all_cases(download_dir: str = None) -> list[str]:
    """Download all supported cases.

    Parameters
    ----------
    download_dir : str
        Base directory to download cases to.

    Examples
    --------
    >>> from ansys.health.heart.utils.download import download_all_cases
    >>> tar_files = download_all_cases("my-downloads")

    To unpack all cases, you can use the ``unpack_cases()`` method:

    >>> from ansys.health.heart.utils.download import unpack_cases
    >>> unpack_cases(tar_files)

    Notes
    -----
    Note that depending on bandwidth, downloading all cases might take a lot of
    time.

    """
    if download_dir is None:
        download_dir = _DOWNLOAD_DIR
    else:
        download_dir = Path(download_dir)

    if not download_dir.exists():
        raise FileExistsError(f"{download_dir} does not exist.")

    tar_files = []
    for database_name, metadata in _ZENODO_RECORDS.items():
        num_cases = metadata["num_cases"]
        for ii in range(1, num_cases + 1):
            LOGGER.info("Downloading {0} : {1}".format(database_name, ii))
            path_to_tar_file = download_case_from_zenodo(database_name, ii, download_dir)
            if path_to_tar_file:
                tar_files.append(str(path_to_tar_file))
    return tar_files


def unpack_cases(list_of_tar_files: list[Path | str]) -> None:
    """Unpack a list of TAR files.

    Parameters
    ----------
    list_of_tar_files : list[Path or str]
        List of TAR files to unpack.

    Examples
    --------
    >>> from ansys.health.heart.utils.download import unpack_cases
    >>> unpack_cases(["01.tar.gz", "02.tar.gz"])
    """
    for file in list_of_tar_files:
        LOGGER.info(f"Unpacking {file}...")
        unpack_case(Path(file))
    return


if __name__ == "__main__":
    download_case_from_zenodo("Rodero2021", 1, ".")
