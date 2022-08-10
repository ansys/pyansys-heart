"""Example to pre-process data from Strocchi2020 and Cristobal2021"""
import os
import pathlib


from ansys.heart.preprocessor.models import LeftVentricle, BiVentricle, FourChamber, FullHeart
from ansys.heart.preprocessor.models import ModelInfo


if __name__ == "__main__":

    # path_to_case = "D:\\development\\pyheart-lib\\pyheart-lib\\downloads\\Strocchi2020\\01\\01.case"
    # workdir = os.path.join(pathlib.Path(path_to_case).parent, "BiVentricleRefactored")
    # info = ModelInfo(database="Strocchi2020", work_directory=workdir, path_to_case=path_to_case)

    # info.clean_workdir(remove_all=True)
    # info.create_workdir()
    # model = BiVentricle(info)
    # model.read_input_mesh()
    # model.extract_simulation_mesh()

    # map data from original mesh to new mesh

    path_to_case = "D:\\development\\pyheart-lib\\pyheart-lib\\downloads\\Strocchi2020\\01\\01.case"
    workdir = os.path.join(pathlib.Path(path_to_case).parent, "FullHeartRefactored")
    info = ModelInfo(database="Strocchi2020", work_directory=workdir, path_to_case=path_to_case)

    info.clean_workdir(remove_all=True)
    info.create_workdir()
    model = FullHeart(info)
    model.read_input_mesh()
    model.extract_simulation_mesh()
    # model.dump_model()

    pass
