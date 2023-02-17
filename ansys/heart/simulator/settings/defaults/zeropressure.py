"""Module contains default values for zero-pressure/stress-free-configuration simulations."""
from pint import Quantity

"""Generic analysis settings."""
analysis = {
    "end_time": Quantity(1000.0, "ms"),
    "dtmin": Quantity(10.0, "ms"),
    "dtmax": Quantity(100.0, "ms"),
    "dt_d3plot": Quantity(100.0, "ms"),
    "dt_nodout": Quantity(200, "ms"),
}
