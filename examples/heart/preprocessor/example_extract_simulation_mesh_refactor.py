"""Example to pre-process data from Strocchi2020 and Cristobal2021"""
import os
import pathlib


from ansys.heart.preprocessor.models import LeftVentricle, BiVentricle, FourChamber, FullHeart
from ansys.heart.preprocessor.models import ModelInfo


if __name__ == "__main__":
    work_dir = pathlib.Path.joinpath(pathlib.Path(__file__).parents[1], "workdir")
    info = ModelInfo(database="Strocchi2020", work_directory=work_dir, path_to_case="test.case")

    model = LeftVentricle(info)

    # get parts:

    pass
