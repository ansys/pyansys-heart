# from conftest import get_workdir, clean_directory
import os

from ansys.heart.preprocessor.input import _InputModel
from ansys.heart.preprocessor.mesh.mesher import (
    mesh_from_manifold_input_model,
    mesh_from_non_manifold_input_model,
)
import numpy as np
import pyvista as pv

from tests.heart.conftest import clean_directory, get_workdir


def test_meshing_for_manifold():
    """Test meshing method for a clean manifold input model."""
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

    write_dir = os.path.join(get_workdir(), "mesher1")

    # write_dir = r"D:\development\pyheart-lib\pyheart-lib\tests\heart\workdir_tests\mesher"
    if not os.path.isdir(write_dir):
        os.makedirs(write_dir)
    else:
        clean_directory(write_dir)

    mesh_file = os.path.join(write_dir, "test_mesh.msh.h5")
    mesh = mesh_from_manifold_input_model(model, write_dir, mesh_file, mesh_size=0.02)

    assert len(mesh.cell_zones) == 3
    assert ["triangles_001", "triangles_002", "triangles_003", "triangles_004"] == [
        fz.name for fz in mesh.face_zones if "interior" not in fz.name
    ]

    os.remove(mesh_file)

    pass


def test_meshing_for_non_manifold():
    """Test meshing method for a dirty non-manifold input model."""

    # prepare mock input data: two non-connecting boxes.
    box1 = pv.Box(bounds=(-1, 0, -1, 1, -1, 1), quads=False)
    box2 = pv.Box(bounds=(0.01, 1, -0.9, 0.9, -0.9, 0.9), quads=False)
    box = box1 + box2
    # box.plot(show_edges=True)

    # split surface by normals so that each side is a separate surface.
    box1.cell_data.set_scalars(int(0), "surface-id")
    box1.compute_normals(feature_angle=89, inplace=True)

    box2.cell_data.set_scalars(int(0), "surface-id")
    box2.compute_normals(feature_angle=89, inplace=True)

    # split input surface by normal direction.
    for ii, normal in enumerate(np.unique(box1.cell_data["Normals"], axis=0)):
        mask = np.all(normal == box1.cell_data["Normals"], axis=1)
        box1.cell_data["surface-id"][mask] = int(ii)

    for ii, normal in enumerate(np.unique(box2.cell_data["Normals"], axis=0)):
        mask = np.all(normal == box2.cell_data["Normals"], axis=1)
        box2.cell_data["surface-id"][mask] = int(ii + 6)

    box = box1 + box2
    box.cell_data["surface-id"] = np.array(box.cell_data["surface-id"], dtype=int)

    # prepare input model.
    model: _InputModel = _InputModel(
        box,
        part_definitions={
            "Part1": {
                "id": 1,
                "enclosed_by_boundaries": {"s1": 0, "s2": 1, "s3": 2, "s4": 3, "s5": 4, "s6": 5},
            },
            "Part2": {
                "id": 2,
                "enclosed_by_boundaries": {
                    "s7": 6,
                    "s8": 7,
                    "s9": 8,
                    "s10": 9,
                    "s11": 10,
                    "s12": 11,
                },
            },
        },
        scalar="surface-id",
    )

    write_dir = os.path.join(get_workdir(), "mesher2")

    clean_directory(write_dir)

    if not os.path.isdir(write_dir):
        os.makedirs(write_dir)

    # call meshing method.
    mesh_file = os.path.join(write_dir, "test_mesh.msh.h5")
    fluent_mesh = mesh_from_non_manifold_input_model(model, write_dir, mesh_file, mesh_size=0.1)

    assert len(fluent_mesh.cell_zones) == 2
    assert sorted(["s1", "s2", "s3", "s4", "s5", "s6", "s8", "s9", "s10", "s11", "s12"]) == sorted(
        [fz.name for fz in fluent_mesh.face_zones if "interior" not in fz.name]
    )

    os.remove(mesh_file)

    pass
