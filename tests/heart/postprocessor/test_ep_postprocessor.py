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

import glob
import os
import tempfile
import unittest.mock as mock

import numpy as np
import pytest
import pyvista as pv
pv.OFF_SCREEN = True

from ansys.heart.core.models import FullHeart
from ansys.heart.postprocessor.ep_postprocessor import EPpostprocessor


#! @kelhouari can we get more sensible mock data?
def _create_mock_ECG_data() -> tuple:
    """Create mock ECG data."""
    data = np.arange(10 * 12).reshape((12, 10)).T
    time = data[:, 0]
    ecg_data = data[:, 1:11]

    return time, ecg_data, data


def _create_mock_electrode_positions() -> np.ndarray:
    """Create mock electrode positions."""
    return pv.Sphere().points[0:12, :]


@pytest.fixture
def _mock_ep_postprocessor(mocker) -> EPpostprocessor:
    """Get mock EP postprocessor object."""
    mocker.patch("ansys.heart.postprocessor.ep_postprocessor.D3plotReader")
    mock_model = mock.Mock(FullHeart)
    return EPpostprocessor(".", mock_model)


@pytest.mark.parametrize("to_plot", [False, True], ids=["Plot=False", "Plot=True"])
def test_compute_12lead_ECG(_mock_ep_postprocessor: EPpostprocessor, to_plot, mocker):
    """Test 12 lead ECG computation."""

    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        # mock post path
        # mock plt.show (surpesses interactive figure)
        mock_post = mocker.patch(
            "ansys.heart.postprocessor.ep_postprocessor.EPpostprocessor.create_post_folder",
            return_value=tempdir,
        )
        mock_show = mocker.patch("ansys.heart.postprocessor.ep_postprocessor.plt.show")

        time, ecg_data, _ = _create_mock_ECG_data()
        ecg = _mock_ep_postprocessor.compute_12_lead_ECGs(ECGs=ecg_data, times=time, plot=to_plot)

        #! TODO: add assertion.
        if to_plot:
            assert os.path.isfile(os.path.join(tempdir, "12LeadECGs.png"))
            mock_post.assert_called_once()
            mock_show.assert_called_once()


def test_read_ECGs(_mock_ep_postprocessor: EPpostprocessor):
    """Read ECGs"""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        data_expected = _create_mock_ECG_data()[-1]
        ecg_data_file = os.path.join(tempdir, "ecg.data")
        np.savetxt(ecg_data_file, data_expected, delimiter=" ", header="h1\nh2\nh3\nh4")
        ecg, times = _mock_ep_postprocessor.read_ECGs(ecg_data_file)
        assert np.allclose(times, data_expected[:, 0])
        assert np.allclose(ecg, data_expected[:, 1:11])


#! TODO: implement sensible asserts.
#! TODO: reduce overlap with test_export_transmembrane_to_vtk
def test_compute_ECGs(_mock_ep_postprocessor: EPpostprocessor, mocker):
    """Test the ECG computation."""
    # mocks the following:
    # vm, times = self.get_transmembrane_potential()
    # self.reader.meshgrid
    # used dummy data for:
    # electrodes: (electrode positions)

    _mock_ep_postprocessor.reader = mock.Mock()
    _mock_ep_postprocessor.reader.meshgrid = pv.examples.load_tetbeam()
    vm = np.ones((10, _mock_ep_postprocessor.reader.meshgrid.n_points))
    times = np.arange(0, 10)

    mock_get_transmembrane = mocker.patch(
        "ansys.heart.postprocessor.ep_postprocessor.EPpostprocessor.get_transmembrane_potential",
        return_value=(vm, times),
    )

    electrodes = _create_mock_electrode_positions()

    ecg_computed, time_computed = _mock_ep_postprocessor.compute_ECGs(electrodes=electrodes)

    mock_get_transmembrane.assert_called_once()

    pass


def test_export_transmembrane_to_vtk(_mock_ep_postprocessor: EPpostprocessor, mocker):
    """Test exporting to VTK."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        _mock_ep_postprocessor.reader = mock.Mock()
        _mock_ep_postprocessor.reader.meshgrid = pv.examples.load_tetbeam()
        vm = np.ones((10, _mock_ep_postprocessor.reader.meshgrid.n_points))
        times = np.arange(0, 10)
        mock_get_transmembrane = mocker.patch(
            "ansys.heart.postprocessor.ep_postprocessor.EPpostprocessor.get_transmembrane_potential",  # noqa E501
            return_value=(vm, times),
        )
        mock_post = mocker.patch(
            "ansys.heart.postprocessor.ep_postprocessor.EPpostprocessor.create_post_folder",
            return_value=tempdir,
        )

        _mock_ep_postprocessor.export_transmembrane_to_vtk()

        mock_post.assert_called_once()
        mock_get_transmembrane.assert_called_once()

        assert len(glob.glob(os.path.join(tempdir, "*.vtk"))) == 10

        #! TODO: do we need asserts?
        _mock_ep_postprocessor.animate_transmembrane()
