# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
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
import hashlib
import os
import shutil
import tempfile
from unittest.mock import MagicMock, Mock, patch

import numpy as np
import pytest
import pyvista as pv
from pyvista import examples as pyvista_examples

import ansys.heart.core.models as models
from ansys.heart.simulator.settings.settings import DynaSettings

# import after mocking.
import ansys.heart.simulator.simulator as simulators
from ansys.heart.simulator.simulator import run_lsdyna


def _get_md5(filename):
    return hashlib.md5(open(filename, "rb").read()).hexdigest()


@pytest.fixture
def base_simulator(mocker) -> simulators.BaseSimulator:
    mocker.patch.object(shutil, "which", return_value=1)
    model = Mock(spec=models.FourChamber).return_value
    model.left_atrium.endocardium = 1
    model.right_atrium.endocardium = 1
    polydata = pv.PolyData(np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]]), [3, 0, 1, 2])
    model.mesh = polydata

    setting = Mock(spec=DynaSettings)
    setting.lsdyna_path = ""
    simulation_directory = "."
    simulator = simulators.BaseSimulator(model, setting, simulation_directory)

    return simulator


@pytest.fixture
def mechanics_simulator(mocker) -> simulators.MechanicsSimulator:
    mocker.patch.object(shutil, "which", return_value=1)
    model = Mock(spec=models.FourChamber).return_value
    model.left_atrium.endocardium = 1
    model.right_atrium.endocardium = 1
    polydata = pv.PolyData(np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]]), [3, 0, 1, 2])
    model.mesh = polydata

    setting = Mock(spec=DynaSettings)
    setting.lsdyna_path = ""
    simulation_directory = "."
    simulator = simulators.MechanicsSimulator(model, setting, simulation_directory)
    return simulator


@pytest.fixture
def mock_laplace():
    with patch.object(
        simulators.BaseSimulator,
        "run_laplace_problem",
    ) as mock_laplace:
        # We only check if run_laplace_problem() is called by correct input
        mock_laplace.side_effect = Exception("ignore output")
        mock_laplace.return_value = "target"
        yield mock_laplace


#! TODO: test does not capture errors in pass through of kwargs
@pytest.mark.parametrize("appendage", [None, [0, 0, 0]])
def test_compute_left_atrial_fiber(base_simulator, mock_laplace, appendage):
    try:
        base_simulator.compute_left_atrial_fiber(appendage=appendage)
    except Exception as e:
        assert str(e) == "ignore output"

    if appendage is None:
        mock_laplace.assert_called_once_with(
            os.path.join(base_simulator.root_directory, "la_fiber"), "la_fiber", laa=None
        )
    else:
        mock_laplace.assert_called_once_with(
            os.path.join(base_simulator.root_directory, "la_fiber"),
            "la_fiber",
            laa=[0, 0, 0],
        )


@pytest.mark.parametrize("top", [None, [1, 0, 0]])
def test_compute_right_atrial_fiber(base_simulator, mock_laplace, top):
    try:
        base_simulator.compute_right_atrial_fiber([0, 0, 0], top=top)
    except Exception as e:
        assert str(e) == "ignore output"

    if top is None:
        mock_laplace.call_args = (
            (
                os.path.join(base_simulator.root_directory, "ra_fiber"),
                "ra_fiber",
            ),
            {"raa": np.array([0, 0, 0]), "top": None},
        )
    else:
        mock_laplace.call_args = (
            (
                os.path.join(base_simulator.root_directory, "ra_fiber"),
                "ra_fiber",
            ),
            {"raa": np.array([0, 0, 0]), "top": [1, 0, 0]},
        )


def test_compute_uvc(base_simulator, mock_laplace):
    try:
        base_simulator.compute_uhc()
    except Exception as e:
        assert str(e) == "ignore output"

    mock_laplace.assert_called_once_with(
        os.path.join(base_simulator.root_directory, "uvc"),
        "uvc",
    )


@pytest.mark.parametrize(
    "simulator_type",
    (
        simulators.BaseSimulator,
        simulators.EPSimulator,
        simulators.EPMechanicsSimulator,
        simulators.MechanicsSimulator,
    ),
)
def test_simulator_inits(mocker, simulator_type):
    """Test inits of all simulators."""
    # mock which
    mocker.patch.object(shutil, "which", return_value=1)
    # test init
    model = MagicMock(spec=models.FourChamber)
    model.workdir = os.getcwd()
    simulator = simulator_type(model=model, dyna_settings=None)

    assert simulator.dyna_settings.__str__() == DynaSettings().__str__()


def test_base_simulator_load_default_settings(mocker):
    """Test loading defaults."""
    from ansys.heart.simulator.settings.settings import SimulationSettings

    mocker.patch.object(shutil, "which", return_value=1)
    # test init
    model = Mock(spec=models.FourChamber)
    model.workdir = os.getcwd()
    model.mesh = pyvista_examples.load_hexbeam()
    model.mesh.cell_data["fiber"] = np.zeros((model.mesh.n_cells, 3), dtype=float)
    model.mesh.cell_data["sheet"] = np.zeros((model.mesh.n_cells, 3), dtype=float)
    fiber = np.eye(3)
    sheet = np.eye(3)

    simulator = simulators.BaseSimulator(model=model, dyna_settings=None)
    simulator.load_default_settings() == SimulationSettings()

    # mock methods
    mocker.patch("ansys.heart.simulator.simulator.BaseSimulator._write_fibers", return_value=".")
    mocker.patch("ansys.heart.simulator.simulator.BaseSimulator._run_dyna", return_value=True)
    mocker.patch(
        "ansys.heart.simulator.simulator._read_orth_element_kfile",
        return_value=(np.array([1, 2, 3]), [], [], fiber, sheet),
    )
    mocker.patch("ansys.heart.core.models.HeartModel.dump_model")

    simulator.compute_fibers()


@patch("subprocess.Popen")
@pytest.mark.parametrize("settings", [None, Mock(DynaSettings)])
def test_run_dyna(mock_subproc_popen, settings):
    """Test run_dyna with mock settings and patched Popen."""
    curr_dir = os.getcwd()
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        tmp_file = os.path.join(tempdir, "main.k")
        with open(tmp_file, "w") as f:
            f.write("*INCLUDE\n")

        run_lsdyna(tmp_file, settings, curr_dir)
        assert mock_subproc_popen.assert_called_once


@pytest.mark.parametrize(
    "folder_name,zerop_folder,auto_post,initial_stress",
    [
        ("main-mechanics1", None, True, False),
        ("main-mechanics1", None, False, False),
        ("main-mechanics1", "zero-pressure", False, True),
    ],
)
#  ids=["combination1, combination2"])
def test_mechanics_simulator_simulate(
    folder_name,
    zerop_folder,
    auto_post,
    initial_stress,
    mocker,
    mechanics_simulator: simulators.MechanicsSimulator,
):
    """Test the simulate method of the mechanics simulator."""
    # set up mocks
    mock_mech_post = mocker.patch("ansys.heart.simulator.simulator.mech_post")
    mock_run_dyna = mocker.patch("ansys.heart.simulator.simulator.MechanicsSimulator._run_dyna")
    mock_write_main = mocker.patch(
        "ansys.heart.simulator.simulator.MechanicsSimulator._write_main_simulation_files"
    )

    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        mechanics_simulator.root_directory = tempdir
        mechanics_simulator.initial_stress = initial_stress
        # create temp files
        if initial_stress:
            zerop_folder = os.path.join(tempdir, "zero-pressure")
            os.makedirs(zerop_folder)

            # generate set of unique lsda files.
            # TODO: may want to have a separate test for this
            # TODO: in order to check the exact condition the
            # TODO: lsda file is used.
            for ii in range(0, 11):
                filename = os.path.join(zerop_folder, f"iter{ii}.dynain.lsda")
                with open(filename, "w") as tmpfile:
                    tmpfile.write(f"lsda_file_{ii}")

            md5_ref = _get_md5(filename)

        mechanics_simulator.simulate(
            folder_name=folder_name, zerop_folder=zerop_folder, auto_post=auto_post
        )

        assert os.path.isdir(os.path.join(tempdir, folder_name))

        mock_run_dyna.assert_called_once()
        mock_write_main.assert_called_once()
        if auto_post:
            mock_mech_post.assert_called_once()
        else:
            mock_mech_post.assert_not_called()

        mock_mech_post.reset_mock()
        mock_run_dyna.reset_mock()
        mock_write_main.reset_mock()

        if initial_stress:
            # TODO: assert if correct file is copied by confirming hash. Will need
            # TODO: unique files for that.
            md5 = _get_md5(os.path.join(tempdir, folder_name, "dynain.lsda"))
            assert md5 == md5_ref
