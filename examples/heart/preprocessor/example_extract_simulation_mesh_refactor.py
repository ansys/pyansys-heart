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
    # model.extract_simulation_mesh()
    # model.dump_model()
    # model.print_info()

    path_to_case = "D:\\development\\pyheart-lib\\pyheart-lib\\downloads\\Strocchi2020\\01\\01.case"
    workdir = os.path.join(pathlib.Path(path_to_case).parent, "FullHeartRefactored")
    info = ModelInfo(database="Strocchi2020", work_directory=workdir, path_to_case=path_to_case)

    info.clean_workdir(remove_all=True)
    info.create_workdir()
    info.dump_info()
    model = FullHeart(info)
    model.extract_simulation_mesh()
    model.dump_model()
    model.print_info()
    model.info.clean_workdir([".stl", ".vtk", ".jou", ".log", ".trn"])

    model.load_model(os.path.join(model.info.workdir, "heart_model.pickle"))
    # model.dump_model()

    pass
