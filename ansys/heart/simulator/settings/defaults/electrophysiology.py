"""Module contains default values for electrophysiology simulations."""
from pint import Quantity

"""Generic analysis settings."""
analysis = {
    "end_time": Quantity(3000.0, "ms"),
    "dtmin": Quantity(10.0, "ms"),
    "dtmax": Quantity(10.0, "ms"),
    "dt_d3plot": Quantity(50.0, "ms"),
}

boundary_conditions = {
    "stimulation_position": [0.0, 0.0, 0.0]
}