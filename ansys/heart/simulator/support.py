"""Provide methods to support workflows."""
import logging
import os
import pathlib as Path

LOGGER = logging.getLogger("pyheart_global.simulator")
import ansys.heart.preprocessor.models.v0_1.models as models


def run_preprocessor(
    model_type: models.HeartModel,
    database: str,
    path_original_mesh: Path,
    work_directory: Path,
    path_to_model: Path = None,
    mesh_size: float = 2.0,
    add_blood_pool: bool = False,
    clean_workdir: bool = True,
):
    """Run the preprocessor with the given input arguments.

    Parameters
    ----------
    model_type : models.HeartModel
        Type of model. Valid values include: LeftVentricle, BiVentricle, FourChamber, FullHeart
    database : str
        Name of the database. Either "Strocchi2020" or "Rodero2021"
    path_original_mesh : Path
        Path to the input mesh file
    work_directory : Path
        Working directory
    path_to_model : Path, optional
        Path to the model, by default None, writes as "heart_model.pickle"
    mesh_size : float, optional
        Size used for remeshing the volume, by default 2.0
    clean_workdir : bool, optional
        Flag indicating whether to clean the working directory, by default True
    """
    if not path_to_model:
        path_to_model = os.path.join(work_directory, "heart_model.pickle")

    LOGGER.info("##############################")
    LOGGER.info("## Launching preprocessor: ###")
    LOGGER.info("##############################")
    LOGGER.info("## Model type: {:<12} ##".format(model_type.__name__))
    LOGGER.info("## Database: {:<14} ##".format(database))
    LOGGER.info("## Remeshing: {:<13} ##".format(str(True)))
    LOGGER.info("##      Size: {:<13} ##".format(mesh_size))
    LOGGER.info("##############################")
    LOGGER.info("## path mesh: %s" % path_original_mesh)
    LOGGER.info("## working directory: %s" % work_directory)
    LOGGER.info("## store model in: %s" % path_to_model)
    LOGGER.info("##############################")

    if not os.path.isdir(work_directory):
        os.makedirs(work_directory)

    # instantiate model information
    info = models.ModelInfo(
        database=database,
        path_to_case=path_original_mesh,
        work_directory=work_directory,
        path_to_model=path_to_model,
        add_blood_pool=add_blood_pool,
    )

    info.mesh_size = mesh_size

    info.clean_workdir(remove_all=True)
    info.create_workdir()
    info.dump_info()

    if model_type == models.LeftVentricle:
        model = models.LeftVentricle(info)

    elif model_type == models.BiVentricle:
        model = models.BiVentricle(info)

    elif model_type == models.FourChamber:
        model = models.FourChamber(info)

    elif model_type == models.FullHeart:
        model = models.FullHeart(info)

    model.extract_simulation_mesh()
    model.dump_model(path_to_model)
    model.print_info()
    if clean_workdir:
        model.info.clean_workdir([".stl", ".vtk", ".jou", ".log"])

    return model
