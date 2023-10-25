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
    "dt_d3plot": Quantity(50.0, "ms"),
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
    prescribed active stress
    TODO
"""

material = {
    "myocardium": {
        "isotropic": {
            "rho": Quantity(0.001, "g/mm^3"),
            "nu": 0.499,
            "k1": Quantity(0.0023599999999999997, "MPa"),
            "k2": Quantity(1.75),
        },
        "anisotropic": {"k1f": Quantity(0.00049, "MPa"), "k2f": Quantity(9.01, "dimensionless")},
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
    "cap": {
        "type": "null",
        "rho": Quantity(0.001, "g/mm^3"),
        "itype": -1,
        "mu1": Quantity(2.0, "MPa"),
        "alpha1": 2,
        "thickness": Quantity(5.0, "mm"),
    },
}

"""Boundary condition settings."""
boundary_conditions = {
    "pericardium": {"penalty_function": [0.65, 25], "spring_stiffness": Quantity(0.05, "MPa/mm")},
    "valve": {
        "biventricle": Quantity(0.002, "MPa/mm"),
        "fourchamber": Quantity(0.02, "MPa/mm"),
        "scale_factor": {"normal": 0.5, "radial": 1.0},
    },
    "end_diastolic_cavity_pressure": {
        "left_ventricle": Quantity(0.002, "MPa"),
        "right_ventricle": Quantity(0.0005333299999999999, "MPa"),
    },
}

"""System model parameters."""
system_model = {
    "name": "ConstantPreloadWindkesselAfterload",
    "left_ventricle": {
        "constants": {
            "Rv": Quantity(6.65e-06, "MPa/(mm^3/ms)"),
            "Ra": Quantity(1.729e-05, "MPa/(mm^3/ms)"),
            "Rp": Quantity(0.000766, "MPa/(mm^3/ms)"),
            "Ca": Quantity(6390000.0, "mm^3/MPa"),
            "Pven": Quantity(0.002, "MPa"),
        },
        "initial_value": {"part": Quantity(0.00931, "MPa")},
    },
    "right_ventricle": {
        "constants": {
            "Rv": Quantity(6.65e-06, "MPa/(mm^3/ms)"),
            "Ra": Quantity(6.0514999999999996e-06, "MPa/(mm^3/ms)"),
            "Rp": Quantity(9.575e-05, "MPa/(mm^3/ms)"),
            "Ca": Quantity(28755000.0, "mm^3/MPa"),
            "Pven": Quantity(0.0005333299999999999, "MPa"),
        },
        "initial_value": {"part": Quantity(0.0019950000000000002, "MPa")},
    },
}
