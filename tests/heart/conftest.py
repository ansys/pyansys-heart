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

import logging as deflogging
import os
import pathlib
import sys

import pytest

import ansys.health.heart.models as models
from ansys.health.heart.utils.download import download_case_from_zenodo, unpack_case

ROOT_FOLDER = os.path.join(pathlib.Path(__file__).parent)

"""

Notes
-----
Note for VS Code/conda users: for the moment it seems that for proper pytest discovery in
VS Code's native testing framework you need to install pyfluent into the
base virtual environment.

"""


def get_assets_folder():
    return os.path.join(ROOT_FOLDER, "assets")


def download_asset(
    database: str = "Strocchi2020", casenumber: int = 1, clean_folder: bool = False
) -> pathlib.Path:
    """Download and unpack the requested asset if it is not yet available."""
    download_dir = os.path.join(get_assets_folder(), "cases")

    if database not in ["Strocchi2020", "Rodero2021"]:
        raise ValueError("Only Strocchi2020 supported for tests.")

    # find case name recursively.
    if database == "Strocchi2020":
        case_path = os.path.join(
            download_dir,
            database,
            "{:02d}".format(casenumber),
            "{:02d}.case".format(casenumber),
        )
    elif database == "Rodero2021":
        case_path = os.path.join(
            download_dir,
            database,
            "{:02d}".format(casenumber),
            "{:02d}.vtk".format(casenumber),
        )

    if os.path.isfile(case_path):
        print("File already exists...")
        return case_path

    print("Downloading asset.")
    path_to_zip = download_case_from_zenodo(database, casenumber, download_dir)
    unpack_case(path_to_zip)

    # remove .vtk file to reduce size (relevant for Github cache)
    if database == "Strocchi2020":
        path_to_vtk = case_path.replace(".case", "-350um.vtk")
    elif database == "Rodero2021":
        path_to_vtk = case_path

    if clean_folder:
        if os.path.isfile(path_to_vtk):
            print(f"Removing .vtk file {path_to_vtk}")
            os.remove(path_to_vtk)

        if os.path.isfile(path_to_zip):
            os.remove(path_to_zip)

    if os.path.isfile(case_path):
        return case_path
    else:
        raise FileExistsError("File not found.")


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


def remove_keys_from_dict(dictionary: dict, exclude_keys=[]):
    """Removes specific keys from the dictionary."""
    new_d = {k: dictionary[k] for k in set(list(dictionary.keys())) - set(exclude_keys)}
    return new_d


@pytest.fixture
def fake_record():
    def inner_fake_record(
        logger,
        msg="This is a message",
        instance_name="172.1.1.1:52000",
        handler_index=0,
        name_logger=None,
        level=deflogging.DEBUG,
        filename="fn",
        lno=0,
        args=(),
        exc_info=None,
        extra={},
    ):
        sinfo = None
        if not name_logger:
            name_logger = logger.name

        if "instance_name" not in extra.keys():
            extra["instance_name"] = instance_name

        record = logger.makeRecord(
            name_logger,
            level,
            filename,
            lno,
            msg,
            args=args,
            exc_info=exc_info,
            extra=extra,
            sinfo=sinfo,
        )
        handler = logger.handlers[handler_index]
        return handler.format(record)

    return inner_fake_record


def is_debugging():
    return "debugpy" in sys.modules


# enable_stop_on_exceptions if the debugger is running during a test
if is_debugging():

    @pytest.hookimpl(tryfirst=True)
    def pytest_exception_interact(call):
        raise call.excinfo.value

    @pytest.hookimpl(tryfirst=True)
    def pytest_internalerror(excinfo):
        raise excinfo.value


def get_fourchamber() -> models.FourChamber:
    vtu_file = os.path.join(
        get_assets_folder(),
        "reference_models",
        "strocchi2020",
        "01",
        "FourChamber",
        "heart_model.vtu",
    )

    json_file = os.path.join(
        get_assets_folder(),
        "reference_models",
        "strocchi2020",
        "01",
        "FourChamber",
        "heart_model.partinfo.json",
    )

    model: models.FourChamber = models.FourChamber(working_directory=".")

    model.load_model_from_mesh(vtu_file, json_file)

    return model


def get_fullheart() -> models.FullHeart:
    vtu_file = os.path.join(
        get_assets_folder(),
        "reference_models",
        "strocchi2020",
        "01",
        "FullHeart",
        "heart_model.vtu",
    )

    json_file = os.path.join(
        get_assets_folder(),
        "reference_models",
        "strocchi2020",
        "01",
        "FullHeart",
        "heart_model.partinfo.json",
    )

    model: models.FullHeart = models.FullHeart(working_directory=".")

    model.load_model_from_mesh(vtu_file, json_file)

    return model
