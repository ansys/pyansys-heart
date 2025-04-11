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

Download a PyAnsys Heart-compatible case from Zenodo
----------------------------------------------------
This example shows how to download a Strocchi 2020 or Rodero 2021 case from the Zenodo
database.
"""

# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths, including that of the working
# directory and generated model.

import glob
from pathlib import Path

from ansys.health.heart.utils.download import download_case_from_zenodo, unpack_case

# Download the TAR file for Rodero 2021 from the Zenodo database.
download_dir = Path.home() / "pyansys-heart" / "downloads"
tar_file = download_case_from_zenodo("Rodero2021", 1, download_dir)

# Unpack the TAR file and get the path to the input CASE or VTK file.
path = unpack_case(tar_file)

print(path)

# list all files
glob.glob("downloads" + "/**/*.*", recursive=True)

###############################################################################
# .. note::
#    You can also manually download the CASE or VTK files from the Strocchi 2020
#    and Rodero 2021 databases. For more information, see:
#
#    - `A Publicly Available Virtual Cohort of Four-chamber Heart Meshes for
#      Cardiac Electro-mechanics Simulations <https://zenodo.org/records/3890034>`_
#    - `Virtual cohort of adult healthy four-chamber heart meshes from CT images <https://zenodo.org/records/4590294>`_
#
#    Alternatively, you can simply click one of the following buttons to download
#    a CASE file in an IPYNB, PY, or ZIP format.
