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

"""Module contains default values for mechanics simulations."""

from pint import Quantity

heart = {
    "cycles": 1,
    "beat_time": Quantity(800, "ms"),
}

"""Generic analysis settings."""
analysis = {
    "end_time": heart["cycles"] * heart["beat_time"],
    "dtmin": Quantity(5.0, "ms"),
    "dtmax": Quantity(5.0, "ms"),
    "dt_d3plot": heart["beat_time"] / 40,
    "dt_icvout": Quantity(5.0, "ms"),
    "global_damping": Quantity(0.1, "1/ms"),
    "stiffness_damping": Quantity(-0.2, "ms"),
}

"""Material settings."""
"""
reference:
    for actype=1:
    https://doi.org/10.1152/ajpheart.01226.2004
    https://doi.org/10.1152/japplphysiol.00255.2014
    for actype=3
    fiber stretch dependency is removed
    prescribed active stress in
    https://doi.org/10.1371/journal.pone.0235145
    HB must be 800ms
"""

material = {
    "myocardium": {
        # # Sack et.al
        # "isotropic": {
        #     "rho": Quantity(0.001, "g/mm^3"),
        #     "kappa": Quantity(1.0, "MPa"),
        #     "k1": Quantity(1.05e-3, "MPa"),
        #     "k2": Quantity(7.542),
        # },
        # "anisotropic": {
        #     "k1f": Quantity(3.465e-3, "MPa"),
        #     "k2f": Quantity(14.472, "dimensionless"),
        #     "k1s": Quantity(0.481e-3, "MPa"),
        #     "k2s": Quantity(12.548, "dimensionless"),
        #     "k1fs": Quantity(0.283e-3, "MPa"),
        #     "k2fs": Quantity(3.088, "dimensionless"),
        # },
        # # ambit
        # "isotropic": {
        #     "rho": Quantity(0.001, "g/mm^3"),
        #     "kappa": Quantity(1.0, "MPa"),
        #     "k1": Quantity(0.059e-3, "MPa"),
        #     "k2": Quantity(8.023),
        # },
        # "anisotropic": {
        #     "k1f": Quantity(18.472e-3, "MPa"),
        #     "k2f": Quantity(16.026, "dimensionless"),
        #     "k1s": Quantity(2.481e-3, "MPa"),
        #     "k2s": Quantity(11.120, "dimensionless"),
        #     "k1fs": Quantity(0.216e-3, "MPa"),
        #     "k2fs": Quantity(11.436, "dimensionless"),
        # },
        # # Gultekin et al.
        # "isotropic": {
        #     "rho": Quantity(0.001, "g/mm^3"),
        #     "kappa": Quantity(1.0, "MPa"),
        #     "k1": Quantity(0.40e-3, "MPa"),
        #     "k2": Quantity(6.55),
        # },
        # "anisotropic": {
        #     "k1f": Quantity(3.05e-3, "MPa"),
        #     "k2f": Quantity(29.05, "dimensionless"),
        #     "k1s": Quantity(1.25e-3, "MPa"),
        #     "k2s": Quantity(36.65, "dimensionless"),
        #     "k1fs": Quantity(0.15e-3, "MPa"),
        #     "k2fs": Quantity(6.28, "dimensionless"),
        # },
        "isotropic": {
            "rho": Quantity(0.001, "g/mm^3"),
            "kappa": Quantity(1.0, "MPa"),
            "k1": Quantity(0.00236, "MPa"),
            "k2": Quantity(1.75),
        },
        "anisotropic": {
            "k1f": Quantity(0.00049, "MPa"),
            "k2f": Quantity(9.01, "dimensionless"),
        },
        # Note: Mechanical simulation uses actype=1, EP-Mechanical simulation uses actype=3
        # related parameters are hard coded
        # For more advanced control, use Material class
        "active": {
            "beat_time": heart["beat_time"],
            "taumax": Quantity(0.125, "MPa"),
            "ss": 0.0,
            "sn": 0.0,
        },
    },
    "passive": {
        "type": "MAT295",
        "rho": Quantity(0.001, "g/mm^3"),
        "itype": -1,
        "kappa": Quantity(10.0, "MPa"),
        "mu1": Quantity(0.1, "MPa"),
        "alpha1": 2,
    },
}

"""Boundary condition settings."""
boundary_conditions = {
    "robin": {
        # pericardium: https://doi.org/10.1016/j.jbiomech.2020.109645
        "ventricle": {
            "penalty_function": [0.25, 25],
            "stiffness": Quantity(0.05, "MPa/mm"),
            "damper": Quantity(0.005, "kPa*s/mm"),
        },
        # from ambit
        "atrial": {
            "stiffness": Quantity(0.075e-3, "MPa/mm"),
            "damper": Quantity(0.005, "kPa*s/mm"),
        },
    },
    "valve": {
        "stiffness": Quantity(0.002, "MPa/mm"),
        "scale_factor": {"normal": 0.5, "radial": 1.0},
    },
    "end_diastolic_cavity_pressure": {
        # TODO: align names with cavity/part names.
        # https://doi.org/10.3389/fphys.2018.00539
        "left_ventricle": Quantity(15, "mmHg"),
        "left_atrial": Quantity(15, "mmHg"),
        "right_ventricle": Quantity(8, "mmHg"),
        "right_atrial": Quantity(8, "mmHg"),
        # # https://doi.org/10.1016/j.jbiomech.2020.109645
        # "left_ventricle": Quantity(18.0, "mmHg"),
        # "right_ventricle": Quantity(9.54, "mmHg"),
    },
}

"""System model parameters."""
# 2wk model, with parameters deduced in a physiological range, from Eindhoven lecture notes
co = Quantity(5, "L/min")
tau = Quantity(2, "s")
pee = Quantity(100, "mmHg")

rp = 0.97 * pee / co
ca = tau / rp
ra = 0.03 * rp
rv = ra

system_model = {
    "name": "open-loop",
    "left_ventricle": {
        "constants": {
            # preload resistance
            "Rv": rv,
            # Z: after load diode resistance
            "Ra": ra,
            # R
            "Rp": rp,
            # C
            "Ca": ca,
            # constant preload, i.e. ED pressure
            "Pven": boundary_conditions["end_diastolic_cavity_pressure"]["left_ventricle"],
        },
        "initial_value": {"part": Quantity(70.0, "mmHg")},
    },
    "right_ventricle": {
        "constants": {
            # preload resistance
            "Rv": rv * 0.5,
            # Z: after load diode resistance
            "Ra": ra * 0.35,
            # R
            "Rp": rp * 0.125,
            # C
            "Ca": ca * 4.5,
            # constant preload, i.e. ED pressure
            "Pven": boundary_conditions["end_diastolic_cavity_pressure"]["right_ventricle"],
        },
        "initial_value": {"part": Quantity(15.0, "mmHg")},
    },
}

# 3wk model found in: https://doi.org/10.1016/j.jbiomech.2020.109645
system_model3 = {
    "name": "ConstantPreloadWindkesselAfterload",
    "left_ventricle": {
        "constants": {
            # Diode resistance
            "Rv": Quantity(0.05, "mmHg*s/mL"),
            # Z
            "Ra": Quantity(0.13, "mmHg*s/mL"),
            # R
            "Rp": Quantity(5.76, "mmHg*s/mL"),
            # C
            "Ca": Quantity(0.85, "mL/mmHg"),
            # constant preload, i.e. ED pressure
            "Pven": boundary_conditions["end_diastolic_cavity_pressure"]["left_ventricle"],
        },
        "initial_value": {"part": Quantity(70.0, "mmHg")},
    },
    "right_ventricle": {
        "constants": {
            # Diode resistance
            "Rv": Quantity(0.05, "mmHg*s/mL"),
            # Z
            "Ra": Quantity(0.13 * 0.35, "mmHg*s/mL"),
            # R
            "Rp": Quantity(5.76 * 0.125, "mmHg*s/mL"),
            # C
            "Ca": Quantity(0.85 * 4.5, "mL/mmHg"),
            # constant preload, i.e. ED pressure
            "Pven": boundary_conditions["end_diastolic_cavity_pressure"]["right_ventricle"],
        },
        "initial_value": {"part": Quantity(15.0, "mmHg")},
    },
}
