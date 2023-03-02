"""Extract a full heart mesh from the public database of 24 pathological hearts
by Strocchi et al (2020)."""

import os
import pathlib

from ansys.heart.misc.downloader import download_case, unpack_case
import ansys.heart.preprocessor.models as models

if __name__ == "__main__":
    # download case from remote repository
    case_num = 1  # patient number 1
    database = "Strocchi2020"
    download_folder: pathlib.Path = os.path.join(pathlib.Path(__file__).parents[3], "download")
    case_path: pathlib.Path = download_case(
        database=database, case_number=case_num, download_folder=download_folder, overwrite=False
    )
    unpack_case(case_path)

    print("done")
