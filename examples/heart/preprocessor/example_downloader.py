"""Extract a full heart mesh from the public database of 24 pathological hearts
by Strocchi et al (2020)."""

import os
import pathlib

from ansys.heart.misc.downloader import download_case, unpack_case
import pyvista as pv

PROJECT_DIRECTORY = pathlib.Path(__file__).absolute().parents[3]

if __name__ == "__main__":
    # download case from remote repository
    case_num = 1  # patient number 1
    database = "Strocchi2020"
    download_folder: pathlib.Path = os.path.join(PROJECT_DIRECTORY, "downloads")
    case_path: pathlib.Path = download_case(
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
