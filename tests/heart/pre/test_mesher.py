# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# from conftest import get_workdir, clean_directory
import glob
import os
import shutil
import tempfile
from unittest import mock

import numpy as np
import pytest
import pyvista as pv

from ansys.health.heart.pre.input import _InputModel
import ansys.health.heart.pre.mesher as mesher

pytestmark = pytest.mark.requires_fluent


@pytest.fixture(scope="session", autouse=True)
def clean_up_temp_dirs():
    tmpdirs = glob.glob(os.path.join(tempfile.gettempdir(), ".pyansys-heart*"))
    yield
    os.chdir(os.path.dirname(__file__))
    import time

    # Sleep to give Fluent process time to shutdown - otherwise tempdirs cannot
    # be removed due to access violation.
    time.sleep(5)

    tmpdirs1 = glob.glob(os.path.join(tempfile.gettempdir(), ".pyansys-heart*"))
    tmp_dirs_remove = list(set(tmpdirs1) - set(tmpdirs))
    for tmp_dir in tmp_dirs_remove:
        try:
            shutil.rmtree(tmp_dir)
        except Exception:
            pass


def test_get_fluent_meshing_session(monkeypatch):
    """Test passing input arguments to the fluent meshing session launcher."""

    expected_keys = [
        "mode",
        "precision",
        "processor_count",
        "start_transcript",
        "ui_mode",
        "product_version",
        "start_container",
    ]

    with mock.patch("ansys.health.heart.pre.mesher.pyfluent.launch_fluent") as mock_launch:
        mesher._uses_container = False
        mesher._get_fluent_meshing_session(".")
        assert list(mock_launch.call_args.kwargs.keys()) == expected_keys

        # Set additional arguments
        mock_launch.reset_mock()
        mesher._extra_launch_kwargs = {"additional_arguments": "-ssh"}
        mesher._get_fluent_meshing_session(".")
        assert mock_launch.call_args.kwargs["additional_arguments"] == "-ssh"

        # Set number of processors
        monkeypatch.delenv(name="PYANSYS_HEART_NUM_CPU", raising=False)
        mock_launch.reset_mock()
        mesher._num_cpus = 4
        mesher._get_fluent_meshing_session(".")
        assert mock_launch.call_args.kwargs["processor_count"] == 4

        # Set number of processors through environment variable.
        monkeypatch.setenv(name="PYANSYS_HEART_NUM_CPU", value="100")
        mock_launch.reset_mock()
        mesher._get_fluent_meshing_session(".")
        assert mock_launch.call_args.kwargs["processor_count"] == 100

        # Use container.
        mock_launch.reset_mock()
        mesher._uses_container = True
        mesher._get_fluent_meshing_session(".")
        assert mock_launch.call_args.kwargs["processor_count"] == 1
        assert mock_launch.call_args.kwargs["start_container"] is True
        assert mock_launch.call_args.kwargs["container_dict"]


@pytest.mark.parametrize(
    "part_names,global_size,size_per_part,expected",
    (
        [["part1", "part2"], 1.0, None, {"part1": 1.0, "part2": 1.0}],
        [["part1", "part2"], 1.0, {"part1": 2.0}, {"part1": 2.0, "part2": 1.0}],
        [["part1", "part2"], 1.0, {"Part1": 2.0}, {"part1": 2.0, "part2": 1.0}],
        [["part1", "part2"], 1.0, {"Part3": 2.0}, {"part1": 1.0, "part2": 1.0}],
        [["Part 1", "part2"], 1.0, {"Part 1": 2.0}, {"part_1": 2.0, "part2": 1.0}],
    ),
)
def test_get_mesh_size_per_part(part_names, global_size, size_per_part, expected):
    """Test setting/overriding mesh size per part."""
    size_per_part1 = mesher._update_size_per_part(part_names, global_size, size_per_part)
    assert size_per_part1 == expected


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

    tmpdir = tempfile.TemporaryDirectory(prefix=".pyansys-heart")

    mesh_file = os.path.join(tmpdir.name, "test_mesh.msh.h5")
    mesh = mesher.mesh_from_manifold_input_model(model, tmpdir.name, mesh_file, mesh_size=0.02)

    assert len(mesh.volume_ids) == 3
    assert ["triangles_001", "triangles_002", "triangles_003", "triangles_004"] == sorted(
        [surface_name for surface_name in mesh.surface_names if "interior" not in surface_name]
    )

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

    tmpdir = tempfile.TemporaryDirectory(prefix=".pyansys-heart")

    # call meshing method.
    mesh_file = os.path.join(tmpdir.name, "test_mesh.msh.h5")

    # specify incomplete mesh-size-per-part. Part 2 should use 0.1 for mesh-size
    vtk_mesh = mesher.mesh_from_non_manifold_input_model(
        model,
        tmpdir.name,
        mesh_file,
        global_mesh_size=0.1,
        _global_wrap_size=0.1,
        mesh_size_per_part={"Part1": 0.15},
    )

    assert len(vtk_mesh.volume_ids) == 2
    assert sorted([volume_name for volume_name in vtk_mesh.volume_names]) == sorted(
        model.part_names
    )
    assert sorted(["s1", "s2", "s3", "s4", "s5", "s6", "s8", "s9", "s10", "s11", "s12"]) == sorted(
        [surface_name for surface_name in vtk_mesh.surface_names if "interior" not in surface_name]
    )
    # assert that more cells exist in Part 2, even though box size is smaller.
    assert vtk_mesh.get_volume(1).n_cells < vtk_mesh.get_volume(2).n_cells

    pass
