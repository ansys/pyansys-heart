"""Module contains default values for mechanics simulations."""
from pint import Quantity

heart = {
    "cycles": 3,
    "beat_time": Quantity(800, "ms"),
}

"""Generic analysis settings."""
analysis = {
    "end_time": heart["cycles"] * heart["beat_time"],
    "dtmin": Quantity(10.0, "ms"),
    "dtmax": Quantity(10.0, "ms"),
    "dt_d3plot": heart["beat_time"] / 20,
    "dt_icvout": Quantity(10.0, "ms"),
    "global_damping": Quantity(0.1, "1/ms"),
}

"""Material settings."""
"""
reference:
    for actype=1:
    https://doi.org/10.1152/ajpheart.01226.2004
    https://doi.org/10.1152/japplphysiol.00255.2014
    for actype=2
    TODO
    for actype=3
    TODO
    prescribed active stress
    https://doi.org/10.1371/journal.pone.0235145
"""

material = {
    "myocardium": {
        # # Sack et.al
        # "isotropic": {
        #     "rho": Quantity(0.001, "g/mm^3"),
        #     "nu": 0.49,
        #     "k1": Quantity(1.05e-3, "MPa"),
        #     "k2": Quantity(7.542),
        # },
        # "anisotropic": {
        #     "k1f": Quantity(3.465e-3, "MPa"),
        #     "k2f": Quantity(14.472, "dimensionless"),
        #     "k1f": Quantity(0.481e-3, "MPa"),
        #     "k2f": Quantity(12.548, "dimensionless"),
        #     "k1f": Quantity(0.283, "MPa"),
        #     "k2f": Quantity(3.088, "dimensionless"),
        # },
        "isotropic": {
            "rho": Quantity(0.001, "g/mm^3"),
            "nu": 0.499,
            "k1": Quantity(0.0023599999999999997, "MPa"),
            "k2": Quantity(1.75),
        },
        "anisotropic": {
            "k1f": Quantity(0.00049, "MPa"),
            "k2f": Quantity(9.01, "dimensionless"),
        },
        "active": {
            "actype": 1,
            "beat_time": heart["beat_time"],
            "taumax": Quantity(0.125, "MPa"),
            "ss": 0.0,
            "sn": 0.0,
        },
    },
    "atrium": {
        "type": "NeoHook",  # or 'MAT295'
        "rho": Quantity(0.001, "g/mm^3"),
        "itype": -1,
        "mu1": Quantity(0.0349, "MPa"),
        "alpha1": 2,
    },
    # null cap: center node is fixed
    # rigid cap: more difficult to converge
    "cap": {
        "type": "null",
        "rho": Quantity(0.001, "g/mm^3"),
    },
}

"""Boundary condition settings."""
# pericardium: https://doi.org/10.1016/j.jbiomech.2020.109645
boundary_conditions = {
    "pericardium": {"penalty_function": [0.25, 25], "spring_stiffness": Quantity(0.05, "MPa/mm")},
    "valve": {
        "stiffness": Quantity(0.002, "MPa/mm"),
        "scale_factor": {"normal": 0.5, "radial": 1.0},
    },
    "end_diastolic_cavity_pressure": {
        ## https://doi.org/10.3389/fphys.2018.00539
        "left_ventricle": Quantity(15, "mmHg"),  # the same for left atrium if exist
        "right_ventricle": Quantity(8, "mmHg"),  # the same for right atrium if exist
        # # https://doi.org/10.1016/j.jbiomech.2020.109645
        # "left_ventricle": Quantity(18.0, "mmHg"),
        # "right_ventricle": Quantity(9.54, "mmHg"),
    },
}

"""System model parameters."""
# 2wk model, with parameters deduced in a physiological range
co = Quantity(5, "L/min")
tau = Quantity(2, "s")
pee = Quantity(100, "mmHg")

rp = 0.97 * pee / co
ca = tau / rp
ra = 0.03 * rp
rv = 0.2 * ra

system_model = {
    "name": "ConstantPreloadWindkesselAfterload",
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
