"""Functional test to determine whether generated biventricle model has all the
expected features."""

import copy
import glob
import os
import pathlib
import shutil
import sys

import ansys.heart.preprocessor.models.v0_1.models as models
from ansys.heart.simulator.support import run_preprocessor
from ansys.heart.writer.dynawriter import (
    ElectrophysiologyDynaWriter,
    FiberGenerationDynaWriter,
    MechanicsDynaWriter,
    PurkinjeGenerationDynaWriter,
    ZeroPressureMechanicsDynaWriter,
)
import pytest

from tests.heart.common import (
    compare_cavity_volume,
    compare_generated_mesh,
    compare_part_names,
    compare_surface_names,
)
from tests.heart.conftest import download_asset, get_assets_folder, get_workdir
from tests.heart.end2end.compare_k import read_file

# marks all tests with the 'requires_fluent' tag after this line
pytestmark = pytest.mark.requires_fluent


# run this fixture first
@pytest.fixture(autouse=True, scope="module")
def extract_bi_ventricle():
    """Extract BiVentricular model which is similar to the reference model.

    Note
    ----
    Do this once as fixture.
    """

    assets_folder = get_assets_folder()
    # path_to_case = os.path.join(assets_folder, "cases", "Strocchi2020", "01", "01.case")
    # path_to_case = os.path.join(assets_folder, "cases", "Rodero2021", "01", "01.vtk")

    path_to_case = download_asset("Strocchi2020", casenumber=1, clean_folder=True)

    # # load model to test against
    # path_to_reference_stats = os.path.join(
    #     assets_folder,
    #     "reference_models",
    #     "strocchi2020",
    #     "01",
    #     "BiVentricle",
    #     "heart_model.pickle",
    # )
    path_to_reference_stats = path_to_reference_stats = os.path.join(
        assets_folder,
        "reference_models",
        "strocchi2020",
        "01",
        "BiVentricle",
        "stats_reference_model_BiVentricle.json",
    )
    assert os.path.isfile(path_to_case)
    assert os.path.isfile(path_to_reference_stats)

    global ref_stats
    import json

    with open(path_to_reference_stats, "r") as openfile:
        ref_stats = json.load(openfile)

    workdir = os.path.join(get_workdir(), "BiVentricle")
    path_to_model = os.path.join(workdir, "heart_model.pickle")

    # copy mesh file
    if not os.path.isdir(workdir):
        os.makedirs(workdir)
    shutil.copyfile(
        os.path.join(
            assets_folder,
            "reference_models",
            "strocchi2020",
            "01",
            "BiVentricle",
            "fluent_volume_mesh.msh.h5",
        ),
        os.path.join(workdir, "fluent_volume_mesh.msh.h5"),
    )

    global model
    model = run_preprocessor(
        model_type=models.BiVentricle,
        database="Strocchi2020",
        path_original_mesh=path_to_case,
        work_directory=workdir,
        path_to_model=path_to_model,
        mesh_size=2,
        skip_meshing=True,
    )

    yield

    # cleanup
    try:
        shutil.rmtree(workdir)
    except:
        print("Failed to cleanup.")


def test_part_names():
    compare_part_names(model, ref_stats)
    pass


def test_surface_names():
    compare_surface_names(model, ref_stats)
    pass


def test_cavities_volumes():
    compare_cavity_volume(model, ref_stats)
    pass


@pytest.mark.xfail(
    sys.platform == "linux", reason="Mesh generation slightly different for Linux version of Fluent"
)
def test_mesh():
    """Test the number of tetrahedrons and triangles in the volume mesh and surface meshes"""
    compare_generated_mesh(model, ref_stats)


# functional tests to determine whether any change was made to
# LS-DYNA input.
@pytest.mark.xfail(reason="This test is currently too strict - rewrite or disable")
@pytest.mark.parametrize(
    "writer_class",
    [
        ElectrophysiologyDynaWriter,
        MechanicsDynaWriter,
        ZeroPressureMechanicsDynaWriter,
        FiberGenerationDynaWriter,
        PurkinjeGenerationDynaWriter,
    ],
)
def test_writers(writer_class):
    """Test whether all writers yield the same .k files as the reference model.

    Notes
    -----
    This skips over most .k files that contain mesh related info.
    """
    writer = writer_class(copy.deepcopy(model))
    ref_folder = os.path.join(
        get_assets_folder(),
        "reference_models",
        "strocchi2020",
        "01",
        "BiVentricle",
        "k_files",
        writer_class.__name__,
    )
    to_test_folder = os.path.join(get_workdir(), "biventricle", writer_class.__name__)
    writer.update()
    writer.export(to_test_folder)

    ref_files = glob.glob(os.path.join(ref_folder, "*.k"))
    # compare each of the reference files to the files that were generated.
    for ref_file in ref_files:
        file_to_compare = os.path.join(to_test_folder, pathlib.Path(ref_file).name)
        assert read_file(ref_file) == read_file(
            file_to_compare
        ), f"File {pathlib.Path(ref_file).name} does not match."

    return
