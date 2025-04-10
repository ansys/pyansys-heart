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

"""

Download a PyAnsys - Heart compatible case from Zenodo.
-------------------------------------------------------
This example shows you how to download a Strocchi 2020 or Rodero 2021 case from the Zenodo
database.
"""

###############################################################################
# .. note::
#    You can also manually download the .case or .vtk files from the Strocchi2020
#    and Rodero2021 databases first. See:
#
#    - https://zenodo.org/records/3890034
#    - https://zenodo.org/records/4590294
#
#    Alternatively you can make use of the download
#    module instead. See the example below.


# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory and generated model.

from pathlib import Path

from ansys.health.heart.utils.download import download_case_from_zenodo, unpack_case

# Download the tar file of Rodero2021 from the Zenodo database.
download_dir = Path.home() / "pyansys-heart" / "downloads"
tar_file = download_case_from_zenodo("Rodero2021", 1, download_dir, overwrite=True)

# Unpack the tar file and get the path to the input .vtk/.case file.
path = unpack_case(tar_file)

print(path)
