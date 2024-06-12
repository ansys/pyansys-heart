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

"""
Methods for reading LS-DYNA binout file.

Note: depend on qd, can be replaced other Ansys's package like by pydpf, lsreader?
"""

import matplotlib.pyplot as plt
import numpy as np

try:
    from qd.cae.Binout import Binout
except ImportError:
    raise ImportError("qd not found.")


class Elout:
    """
    Class of element output.

    Read stress and history variables (28) for solid elements.
    """

    def __init__(self, fn):
        """
        Init Elout.

        Parameters
        ----------
        fn: binout file from LS-DYNA
        """
        self.binout = Binout(fn)
        self.time = self.binout.read("elout", "solid_hist", "time")
        self.ids = self.binout.read("elout", "solid_hist", "ids")[0]
        self.solid_nb = len(self.ids)

    def get_stress(self):
        """
        Get stress.

        Returns
        -------
        numpy array of time * element * 6 stress.
        stress-xx   stress-yy   stress-zz   stress-xy   stress-yz   stress-zx
        """
        data = self.binout.read("elout", "solid_hist", "data").reshape(len(self.time), -1, 9)
        if not np.all(data[:, :, 0] == 4) or data.shape[1] != self.solid_nb:
            raise Exception("Cannot read stress correctly.")
        else:
            return data[:, :, 1:7]

    def get_history_variable(self):
        """
        Get history variables.

        Returns
        -------
        numpy array of time * element * 27 history variables.
        """
        data = self.binout.read("elout", "solid_hist", "hist").reshape(len(self.time), -1, 28)
        if not np.all(data[:, :, 0] == 4) or data.shape[1] != self.solid_nb:
            raise Exception("Cannot read hv correctly.")
        else:
            return data[:, :, 1:]


class NodOut:
    """Nodout class."""

    def __init__(self, fn):
        """Init NodOut."""
        self.binout = Binout(fn)
        self.time = self.binout.read("nodout", "time")
        self.ids = self.binout.read("nodout", "ids")

        # Suppose all nodes are saved in order
        check_ok = np.array_equal(self.ids, np.linspace(1, len(self.ids + 1), num=len(self.ids)))
        if not check_ok:
            Exception("Nodout file is not complete, cannot continue...")

    def get_coordinates(self):
        """Get coordinates."""
        x_coord = self.binout.read("nodout", "x_coordinate")
        y_coord = self.binout.read("nodout", "y_coordinate")
        z_coord = self.binout.read("nodout", "z_coordinate")
        coords = np.stack((x_coord, y_coord, z_coord), axis=2)
        return coords


class IcvOut:
    """IcvOut class."""

    def __init__(self, fn):
        """Init IcvOut."""
        self.binout = Binout(fn)
        time = self.binout.read("icvout", "time")  # s
        pressure = self.binout.read("icvout", "ICV_Pressure")
        volume = self.binout.read("icvout", "ICV_Volume")
        flow = self.binout.read("icvout", "ICVI_flow_rate")

        # Fix small bug in LSDYNA Output: Volume is 0 at t0
        # V1 = V0 + dt * (1-gamma) * Q0
        # negative flow is inflow for the cavity
        # 0.4 is from 1-gamma
        volume[0] = volume[1] + (time[1] - time[0]) * flow[0] * 0.4

        self.time = time
        self.pressure = pressure.reshape(-1, pressure.ndim)
        self.flow = flow.reshape(-1, pressure.ndim)
        self.volume = volume.reshape(-1, pressure.ndim)


if __name__ == "__main__":
    a = NodOut("binout0000")
    from ansys.heart.preprocessor.models import HeartModel

    model = HeartModel.load_model("heart_model.pickle")

    apex_id = model.left_ventricle.apex_points[0].node_id
    mv_ids = model.left_ventricle.caps[0].node_ids
    coords = a.get_coordinates()
    apex = coords[:, apex_id, :]
    mv_center = np.mean(coords[:, mv_ids, :], axis=1)
    dst = np.linalg.norm(mv_center - apex, axis=1)
    plt.plot(a.time[2:], dst[2:])
    plt.show()
