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

import glob
import os
import pathlib
import tempfile

import pytest

from ansys.health.heart.writer.laplace import LaplaceWriter
from tests.heart.conftest import get_assets_folder, get_fourchamber
from tests.heart.end2end.compare_k import read_file


@pytest.fixture
def fourchamber():
    return get_fourchamber()


@pytest.mark.k_file_writer
def test_uvc(fourchamber):
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:
        to_test_folder = os.path.join(workdir, "uvc")
        fourchamber._workdir = workdir
        writer = LaplaceWriter(fourchamber, "uvc")
        writer.update()
        writer.export(to_test_folder)

        ref_folder = os.path.join(
            get_assets_folder(),
            "reference_models",
            "strocchi2020",
            "01",
            "FourChamber",
            "uvc",
        )
        compare_outputs(to_test_folder, ref_folder)


@pytest.mark.k_file_writer
def test_la_fiber(fourchamber):
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:
        to_test_folder = os.path.join(workdir, "la_fiber")
        writer = LaplaceWriter(fourchamber, "la_fiber")
        writer.update()
        writer.export(to_test_folder)

        ref_folder = os.path.join(
            get_assets_folder(),
            "reference_models",
            "strocchi2020",
            "01",
            "FourChamber",
            "la_fiber",
        )
        compare_outputs(to_test_folder, ref_folder)


@pytest.mark.k_file_writer
def test_ra_top(fourchamber):
    start = [-46.4863, 133.29, 405.961]
    end = [-33.7271, 134.605, 332.155]
    mid = [-43.5525, 148.992, 367.271]

    writer = LaplaceWriter(fourchamber, "ra_fiber", raa=[-33, 82, 417], top=[start, mid, end])
    ids = writer._find_top_nodeset_by_geodesic(writer.target)

    assert len(ids) == 61
    assert ids[0] == 1224
    assert ids[-1] == 23575


@pytest.mark.k_file_writer
def test_ra_fiber(fourchamber):
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:
        to_test_folder = os.path.join(workdir, "ra_fiber")
        writer = LaplaceWriter(fourchamber, "ra_fiber", raa=[-33, 82, 417])
        writer.update()
        writer.export(to_test_folder)

        ref_folder = os.path.join(
            get_assets_folder(),
            "reference_models",
            "strocchi2020",
            "01",
            "FourChamber",
            "ra_fiber",
        )
        compare_outputs(to_test_folder, ref_folder)


@pytest.mark.k_file_writer
def test_drbm(fourchamber):
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as workdir:
        to_test_folder = os.path.join(workdir, "drbm")
        fourchamber._workdir = workdir
        writer = LaplaceWriter(fourchamber, "D-RBM")
        writer.update()
        writer.export(to_test_folder)

        ref_folder = os.path.join(
            get_assets_folder(),
            "reference_models",
            "strocchi2020",
            "01",
            "FourChamber",
            "drbm",
        )
        compare_outputs(to_test_folder, ref_folder)


def compare_outputs(to_test_folder, ref_folder):
    ref_files = glob.glob(os.path.join(ref_folder, "*.k"))
    # compare each of the reference files to the files that were generated.
    for ref_file in ref_files:
        file_to_compare = os.path.join(to_test_folder, pathlib.Path(ref_file).name)
        assert read_file(ref_file) == read_file(file_to_compare), (
            f"File {pathlib.Path(ref_file).name} does not match."
        )
