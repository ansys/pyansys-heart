# from conftest import get_workdir, clean_directory
import os

from ansys.heart.preprocessor.input import _InputModel
from ansys.heart.preprocessor.mesh.fluenthdf5 import FluentMesh
from ansys.heart.preprocessor.mesh.mesher import mesh_from_good_quality_input_model
import pyvista as pv

from tests.heart.conftest import clean_directory, get_workdir

# import shutil


def test_meshing_from_clean():
    """Prep surface/boundary data for creating different cell zones."""
    sphere1 = pv.Sphere()
    sphere2 = pv.Sphere()
    sphere2.points += [0.2, 0.0, 0.0]
    sphere = sphere1 + sphere2
    union = sphere1.boolean_union(sphere2)
    intersection = sphere1.boolean_intersection(sphere2)
    sphere = union + intersection

    intersection.cell_data["CellSource"][intersection.cell_data["CellSource"] == 0] = 2
    intersection.cell_data["CellSource"][intersection.cell_data["CellSource"] == 1] = 3

    sphere = union + intersection

    model: _InputModel = _InputModel(
        sphere,
        part_definitions={
            "Part1": {"id": 1, "enclosed_by_boundaries": {"triangles_001": 0, "triangles_002": 3}},
            "Part2": {"id": 2, "enclosed_by_boundaries": {"triangles_003": 1, "triangles_004": 2}},
            "Part3": {"id": 3, "enclosed_by_boundaries": {"triangles_004": 2, "triangles_005": 3}},
        },
        scalar="CellSource",
    )

    write_dir = os.path.join(get_workdir(), "mesher")

    clean_directory(write_dir)

    # write_dir = r"D:\development\pyheart-lib\pyheart-lib\tests\heart\workdir_tests\mesher"
    if not os.path.isdir(write_dir):
        os.makedirs(write_dir)

    model.write_part_boundaries(write_dir)

    mesh_file = os.path.join(write_dir, "test_mesh.msh.h5")
    mesh_from_good_quality_input_model(
        model,
        write_dir,
        mesh_file,
    )

    # check if mesh is as expected.
    mesh = FluentMesh()
    mesh.load_mesh(mesh_file)

    assert len(mesh.cell_zones) == 3
    assert ["triangles_001", "triangles_002", "triangles_003", "triangles_004"] == [
        fz.name for fz in mesh.face_zones if "interior" not in fz.name
    ]

    os.remove(mesh_file)

    pass
