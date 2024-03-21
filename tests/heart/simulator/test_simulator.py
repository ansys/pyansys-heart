from unittest.mock import Mock, patch

from ansys.heart.preprocessor.models.v0_1.models import FourChamber
from ansys.heart.simulator.settings.settings import DynaSettings
from ansys.heart.simulator.simulator import BaseSimulator
import numpy as np
import pytest
import pyvista as pv


@pytest.fixture
def mock_which():
    # otherwise, need a legal dyna path
    with patch("ansys.heart.simulator.simulator.which") as mocked_which:
        # whatever return value
        mock_which.return_value = 1
        yield mock_which


@pytest.fixture
def mocked_BaseSimulator(mock_which) -> BaseSimulator:
    model = Mock(spec=FourChamber).return_value
    model.left_atrium.endocardium = 1
    model.right_atrium.endocardium = 1
    polydata = pv.PolyData(np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]]), [3, 0, 1, 2])
    model.mesh = polydata

    setting = Mock(spec=DynaSettings)
    setting.lsdyna_path = ""
    simulation_directory = "dummy"
    simulator = BaseSimulator(model, setting, simulation_directory)
    # simulator.run_laplace_problem = MagicMock(return_value=polydata)

    return simulator


@pytest.mark.parametrize("appendage", [None, [0, 0, 0]])
def test_compute_left_atrial_fiber(mocked_BaseSimulator, appendage):
    with patch(
        "ansys.heart.simulator.simulator.BaseSimulator.compute_left_atrial_fiber",
        side_effect=Mock(),
    ) as mock_la_fiber:
        mocked_BaseSimulator.compute_left_atrial_fiber(appendage=appendage)

        if appendage is None:
            mock_la_fiber.assert_called_once()
        else:
            mock_la_fiber.assert_called_once_with(appendage=[0, 0, 0])
