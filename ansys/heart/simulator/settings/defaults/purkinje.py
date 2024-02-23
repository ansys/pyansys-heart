"""Module contains default values for Purkinje generation."""
from pint import Quantity

"""Construction parameters."""
build = {
    "edgelen": Quantity(0.5, "mm"),
    "ngen": Quantity(200, "dimensionless"),
    "nbrinit": Quantity(3, "dimensionless"),
    "nsplit": Quantity(2, "dimensionless"),
    "pmjtype": Quantity(1, "dimensionless"),
    "pmjradius": Quantity(0.7, "dimensionless"),
}
