"""Module contains default values for fiber generation."""
from pint import Quantity

"""Fiber angles."""
angles = {
    "alpha_endo": Quantity(-60, "degree"),
    "alpha_epi": Quantity(60, "degree"),
    "beta_endo": Quantity(-65, "degree"),
    "beta_epi": Quantity(25, "degree"),
    "beta_endo_septum": Quantity(-65, "degree"),
    "beta_epi_septum": Quantity(25, "degree"),
}
