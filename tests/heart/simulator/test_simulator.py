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
import os
from unittest.mock import Mock, patch

import numpy as np
import pytest
import pyvista as pv

from ansys.heart.core.models import FourChamber
from ansys.heart.simulator.settings.settings import DynaSettings
from ansys.heart.simulator.simulator import BaseSimulator


@pytest.fixture
def mock_which():
    # mock lsdyna path check
    with patch("ansys.heart.simulator.simulator.which") as mocked_which:
        # whatever return value
        mock_which.return_value = 1
        yield mock_which


@pytest.fixture
def simulator(mock_which) -> BaseSimulator:
    model = Mock(spec=FourChamber).return_value
    model.left_atrium.endocardium = 1
    model.right_atrium.endocardium = 1
    polydata = pv.PolyData(np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]]), [3, 0, 1, 2])
    model.mesh = polydata

    setting = Mock(spec=DynaSettings)
    setting.lsdyna_path = ""
    simulation_directory = "."
    simulator = BaseSimulator(model, setting, simulation_directory)
    # simulator.run_laplace_problem = MagicMock(return_value=polydata)

    return simulator


@pytest.fixture
def mock_laplace():
    with patch.object(
        BaseSimulator,
        "run_laplace_problem",
    ) as mock_laplace:
        # We only check if run_laplace_problem() is called by correct input
        mock_laplace.side_effect = Exception("ignore output")
        mock_laplace.return_value = "target"
        yield mock_laplace


@pytest.mark.parametrize("appendage", [None, [0, 0, 0]])
def test_compute_left_atrial_fiber(simulator, mock_laplace, appendage):

    try:
        simulator.compute_left_atrial_fiber(appendage=appendage)
    except Exception as e:
        assert str(e) == "ignore output"

    if appendage is None:
        mock_laplace.assert_called_once_with(
            os.path.join(simulator.root_directory, "la_fiber"), "la_fiber", laa=None
        )
    else:
        mock_laplace.assert_called_once_with(
            os.path.join(simulator.root_directory, "la_fiber"),
            "la_fiber",
            laa=[0, 0, 0],
        )


@pytest.mark.xfail(reason="assert_called_once_with get an error with these inputs")
@pytest.mark.parametrize("top", [None, [1, 0, 0]])
def test_compute_right_atrial_fiber(simulator, mock_laplace, top):

    try:
        simulator.compute_right_atrial_fiber([0, 0, 0], top=top)
    except Exception as e:
        assert str(e) == "ignore output"

    if top is None:
        mock_laplace.assert_called_once_with(
            os.path.join(simulator.root_directory, "ra_fiber"),
            "ra_fiber",
            raa=np.array([0, 0, 0]),
            top=None,
        )
    else:
        mock_laplace.assert_called_once_with(
            os.path.join(simulator.root_directory, "ra_fiber"),
            "ra_fiber",
            raa=np.array([0, 0, 0]),
            top=[1, 0, 0],
        )


def test_compute_uvc(simulator, mock_laplace):

    try:
        simulator.compute_uhc()
    except Exception as e:
        assert str(e) == "ignore output"

    mock_laplace.assert_called_once_with(
        os.path.join(simulator.root_directory, "uvc"),
        "uvc",
    )
