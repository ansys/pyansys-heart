"""Module contains default values for fiber generation."""
from pint import Quantity

"""Fiber angles."""
angles = {
    "alpha_endo": Quantity(-60, "dimensionless"),
    "alpha_epi": Quantity(60, "dimensionless"),
    "beta_endo": Quantity(-65, "dimensionless"),
    "beta_epi": Quantity(25, "dimensionless"),
    "beta_endo_septum": Quantity(-65, "dimensionless"),
    "beta_epi_septum": Quantity(25, "dimensionless"),
}
