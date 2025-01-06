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

"""Module contains default values for fiber generation."""

from pint import Quantity

"""Ventricle fiber angles."""
angles = {
    "alpha_endo": Quantity(-60, "degree"),
    "alpha_epi": Quantity(60, "degree"),
    "beta_endo": Quantity(-65, "degree"),
    "beta_epi": Quantity(25, "degree"),
    "beta_endo_septum": Quantity(-65, "degree"),
    "beta_epi_septum": Quantity(25, "degree"),
}

# From https://doi.org/10.1016/j.cma.2020.113468, table 4, for idealized geometry
# paper has a typo, lpv should smaller than rpv
la_bundle = {"tau_mv": 0.65, "tau_lpv": 0.10, "tau_rpv": 0.65}

ra_bundle = {
    "tau_tv": 0.9,
    "tau_raw": 0.55,
    "tau_ct_minus": -0.18,
    "tau_ct_plus": -0.1,
    "tau_icv": 0.9,
    "tau_scv": 0.1,
    "tau_ib": 0.135,  # paper has a typo, ras should larger than ib
    "tau_ras": 0.35,
}
