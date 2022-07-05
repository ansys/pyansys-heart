from .example_create_lsdyna_files import create_ls_dyna_files

if __name__ == "__main__":

    path_model_info = (
        "D:\\development\\pyheart-lib\\pyheart-lib\\downloads\\Strocchi2020_simplified\\workdir\\"
    )

    create_ls_dyna_files(path_model_info, writer_type="Mechanics")
    create_ls_dyna_files(path_model_info, writer_type="ZeroPressure")
    create_ls_dyna_files(path_model_info, writer_type="FiberGeneration")
