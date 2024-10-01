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

import hashlib
import json
import os

from ansys.heart.misc.downloader import download_case_from_zenodo

if __name__ == "__main__":
    """Creates a reference hash table used to check integrity of data files.

    Notes
    -----
    Use this script to download all .tar.gz files from both the Strocchi2020 and
    Rodero2021 public datasets and generate a hash table which is used to
    verify file integrity.
    """
    databases = {"Strocchi2020": 24, "Rodero2021": 20}
    base_folder = "D:\\development\\PyAnsys-Heart\\PyAnsys-Heart\\downloads"
    path_to_hash_table = os.path.join(base_folder, "remote_repo_hash_table_sha256.json")
    if os.path.isfile(path_to_hash_table):
        # read
        fid = open(path_to_hash_table, "r")
        sha256_table: dict = json.load(fid)
        fid.close()
        # convert strings to ints
        for key in sha256_table.keys():
            sha256_table[key] = {int(k): v for k, v in sha256_table[key].items()}
    else:
        sha256_table = {}

    for database in databases:
        num_cases = databases[database]
        download_folder = os.path.join(base_folder)
        try:
            sha256_table[database]
        except KeyError:
            sha256_table[database] = {}

        for ii in range(1, num_cases + 1):
            try:
                sha256_table[database][ii]
                key_exists = True
                print("{0}:{1} Entry already exists".format(database, ii))
                continue
            except KeyError:
                key_exists = False

            print("\nDownloading... %d" % ii)
            path_to_case = download_case_from_zenodo(
                database, ii, download_folder, validate_hash=False
            )
            sha256 = hashlib.sha256(open(path_to_case, "rb").read()).hexdigest()
            sha256_table[database][ii] = sha256

            with open(path_to_hash_table, "w") as outfile:
                json.dump(sha256_table, outfile, indent=4, ensure_ascii=True)
