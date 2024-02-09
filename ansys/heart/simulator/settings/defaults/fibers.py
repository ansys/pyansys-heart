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

# Quarteroni idealized geometry
# paper has a type, lpv should < rpv
la_bundle = {"tau_mv": 0.65, "tau_lpv": 0.10, "tau_rpv": 0.65}

# la_bundle = {"tau_mv": 0.85, "tau_lpv": 0.4, "tau_rpv": 0.9}

ra_bundle = {
    "tau_tv ": 0.9,
    "tau_raw ": 0.55,
    "tau_ct_minus ": -0.18,
    "tau_ct_plus ": -0.1,
    "tau_icv ": 0.9,
    "tau_scv ": 0.1,
    "tau_ib ": 0.135,  # paper has a typo, ras should > ras
    "tau_ras ": 0.35,  # paper has a typo, ras should > ras
}
