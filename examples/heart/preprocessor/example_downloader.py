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

"""Extract a full heart mesh from the public database of 24 pathological hearts
by Strocchi et al (2020)."""

import os
import pathlib

import pyvista as pv

from ansys.heart.core.downloader import download_case_from_zenodo, unpack_case

PROJECT_DIRECTORY = pathlib.Path(__file__).absolute().parents[3]

if __name__ == "__main__":
    # download case from remote repository
    case_num = 1  # patient number 1
    database = "Strocchi2020"
    download_folder: pathlib.Path = os.path.join(PROJECT_DIRECTORY, "downloads")
    case_path: pathlib.Path = download_case_from_zenodo(
        database=database, case_number=case_num, download_folder=download_folder, overwrite=False
    )
    unpack_case(case_path)
    mesh_path = os.path.join(
        pathlib.Path(case_path).parents[0], "%02d" % (case_num,), "%02d.case" % (case_num,)
    )
    mesh = pv.read(mesh_path)

    plotter = pv.Plotter()
    plotter.add_mesh(mesh, color="white")
    plotter.show()
    print("done")
