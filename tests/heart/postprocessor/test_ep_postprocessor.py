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

import glob
import os
import tempfile
import unittest.mock as mock

import numpy as np
import pytest
import pyvista as pv

from ansys.heart.core.models import FullHeart
from ansys.heart.core.postprocessor.dpf_utils import EPpostprocessor


# TODO: test_compute_12lead_ECG requires better asserts and we may
# TODO: want to add better test data that allows testing internals
# TODO: of the method.
# TODO: see issue #973
def _create_mock_ECG_data() -> tuple:  # noqa N802
    """Create mock ECG data."""
    data = np.arange(10 * 12).reshape((12, 10)).T
    time = data[:, 0]
    ecg_data = data[:, 1:11]

    return time, ecg_data, data


def _create_mock_electrode_positions() -> np.ndarray:
    """Create mock electrode positions."""
    return pv.Sphere().points[0:12, :]


@pytest.fixture
def _mock_ep_postprocessor():
    """Get mock EP postprocessor object."""
    with mock.patch("ansys.heart.core.postprocessor.dpf_utils.D3plotReader"):
        mock_model = mock.Mock(FullHeart)

        yield EPpostprocessor(".", mock_model)


@pytest.mark.parametrize("to_plot", [False, True], ids=["Plot=False", "Plot=True"])
def test_compute_12lead_ECG(to_plot, _mock_ep_postprocessor: EPpostprocessor):  # noqa: N802
    """Test 12 lead ECG computation."""
    time, ecg_data, _ = _create_mock_ECG_data()

    with mock.patch(
        "ansys.heart.core.postprocessor.dpf_utils.EPpostprocessor.create_post_folder",
        return_value="",
    ) as mock_post:
        with mock.patch("ansys.heart.core.postprocessor.dpf_utils.plt.show") as mock_show:
            with mock.patch("ansys.heart.core.postprocessor.dpf_utils.plt.savefig") as mock_save:
                _mock_ep_postprocessor.compute_12_lead_ECGs(ECGs=ecg_data, times=time, plot=to_plot)
                if to_plot:
                    mock_post.assert_called_once()
                    mock_show.assert_called_once()
                    mock_save.assert_called_once()


def test_read_ECGs(_mock_ep_postprocessor: EPpostprocessor):  # noqa N802
    """Read ECGs."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        data_expected = _create_mock_ECG_data()[-1]
        ecg_data_file = os.path.join(tempdir, "ecg.data")
        np.savetxt(ecg_data_file, data_expected, delimiter=" ", header="h1\nh2\nh3\nh4")
        ecg, times = _mock_ep_postprocessor.read_ECGs(ecg_data_file)
        assert np.allclose(times, data_expected[:, 0])
        assert np.allclose(ecg, data_expected[:, 1:11])


# TODO: implement sensible asserts.
# TODO: reduce overlap with test_export_transmembrane_to_vtk
def test_compute_ECGs(_mock_ep_postprocessor: EPpostprocessor):  # noqa N802
    """Test the ECG computation."""
    _mock_ep_postprocessor.reader = mock.Mock()
    _mock_ep_postprocessor.reader.meshgrid = pv.examples.load_tetbeam()
    vm = np.ones((10, _mock_ep_postprocessor.reader.meshgrid.n_points))
    times = np.arange(0, 10)

    with mock.patch(
        "ansys.heart.core.postprocessor.dpf_utils.EPpostprocessor.get_transmembrane_potential",
        return_value=(vm, times),
    ) as mock_get_transmembrane:
        electrodes = _create_mock_electrode_positions()

        ecg_computed, time_computed = _mock_ep_postprocessor.compute_ECGs(electrodes=electrodes)

        mock_get_transmembrane.assert_called_once()

    pass


@mock.patch("pyvista.Plotter.show")
@mock.patch("pyvista.Plotter.update_scalars")
@mock.patch("pyvista.Plotter.update")
def test_export_transmembrane_to_vtk(
    mock_update, mock_update_scalars, mock_show, _mock_ep_postprocessor: EPpostprocessor
):
    """Test exporting to VTK."""
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        _mock_ep_postprocessor.reader = mock.Mock()
        _mock_ep_postprocessor.reader.meshgrid = pv.examples.load_tetbeam()
        vm = np.ones((10, _mock_ep_postprocessor.reader.meshgrid.n_points))
        times = np.arange(0, 10)

        with mock.patch(
            "ansys.heart.core.postprocessor.dpf_utils.EPpostprocessor.get_transmembrane_potential",  # noqa E501
            return_value=(vm, times),
        ) as mock_get_transmembrane:
            with mock.patch(
                "ansys.heart.core.postprocessor.dpf_utils.EPpostprocessor.create_post_folder",
                return_value=tempdir,
            ) as mock_post:
                _mock_ep_postprocessor.export_transmembrane_to_vtk()

                mock_post.assert_called_once()
                mock_get_transmembrane.assert_called_once()

                assert len(glob.glob(os.path.join(tempdir, "*.vtk"))) == 10

                _mock_ep_postprocessor.animate_transmembrane()

                assert mock_update.call_count == vm.shape[0]
                assert mock_update_scalars.call_count == vm.shape[0]
                mock_show.assert_called_once()
