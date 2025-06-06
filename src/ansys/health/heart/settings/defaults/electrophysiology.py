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

"""Module contains default values for EP (Electrophysiology) simulations."""

from pint import Quantity

heart = {
    "cycles": 1,
    "beat_time": Quantity(800, "ms"),
}

"""Generic analysis settings."""
analysis = {
    "end_time": heart["cycles"] * heart["beat_time"],
    "dtmin": Quantity(0.0, "ms"),
    "dtmax": Quantity(1.0, "ms"),
    "dt_d3plot": Quantity(10, "ms"),
    "solvertype": "Monodomain",
}


"""Material settings."""
material = {
    "myocardium": {
        "velocity_fiber": Quantity(0.7, "mm/ms"),  # mm/ms in case of eikonal model
        "velocity_sheet": Quantity(0.2, "mm/ms"),  # mm/ms in case of eikonal model
        "velocity_sheet_normal": Quantity(0.2, "mm/ms"),  # mm/ms in case of eikonal model
        "sigma_fiber": Quantity(0.5, "mS/mm"),  # mS/mm
        "sigma_sheet": Quantity(0.1, "mS/mm"),  # mS/mm
        "sigma_sheet_normal": Quantity(0.1, "mS/mm"),  # mS/mm
        "sigma_passive": Quantity(1.0, "mS/mm"),  # mS/mm
        "beta": Quantity(140, "1/mm"),
        "cm": Quantity(0.01, "uF/mm^2"),  # uF/mm^2
        "lambda": Quantity(0.2, "dimensionless"),
        "percent_endo": Quantity(0.17, "dimensionless"),
        "percent_mid": Quantity(0.41, "dimensionless"),
    },
    "beam": {
        "velocity": Quantity(1, "mm/ms"),  # mm/ms in case of eikonal model
        "sigma": Quantity(1, "mS/mm"),  # mS/mm
        "beta": Quantity(140, "1/mm"),
        "cm": Quantity(0.001, "uF/mm^2"),  # uF/mm^2
        "lambda": Quantity(0.2, "dimensionless"),
        "pmjrestype": Quantity(1),
        "pmjres": Quantity(0.001, "1/mS"),  # 1/mS
    },
}

"""Stimulation settings."""
stimulation = {
    "stimdefaults": {
        "node_ids": None,
        "t_start": Quantity(0.0, "ms"),
        "period": Quantity(800.0, "ms"),
        "duration": Quantity(2, "ms"),
        "amplitude": Quantity(50, "uF/mm^3"),
    },
}
