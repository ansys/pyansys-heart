"""Module contains default values for electrophysiology simulations."""
from pint import Quantity

heart = {
    "cycles": 1,
    "beat_time": Quantity(800, "ms"),
}

"""Generic analysis settings."""
analysis = {
    "end_time": heart["cycles"] * heart["beat_time"],
    "dtmin": Quantity(0.0, "ms"),
    "dtmax": Quantity(1.0, "ms"),
    "dt_d3plot": Quantity(10, "ms"),
    "solvertype": "Monodomain",
}


"""Material settings."""
material = {
    "myocardium": {
        "velocity_fiber": Quantity(1, "mm/ms"),  # mm/ms in case of eikonal model
        "velocity_sheet": Quantity(1, "mm/ms"),  # mm/ms in case of eikonal model
        "velocity_sheet_normal": Quantity(1, "mm/ms"),  # mm/ms in case of eikonal model
        "sigma_fiber": Quantity(0.5, "mS/mm"),  # mS/mm
        "sigma_sheet": Quantity(0.1, "mS/mm"),  # mS/mm
        "sigma_sheet_normal": Quantity(0.1, "mS/mm"),  # mS/mm
        "sigma_passive": Quantity(1.0, "mS/mm"),  # mS/mm
        "beta": Quantity(140, "1/mm"),
        "cm": Quantity(0.01, "uF/mm^2"),  # uF/mm^2
        "lambda": Quantity(0.2, "dimensionless"),
    },
    "beam": {
        "velocity": Quantity(1, "mm/ms"),  # mm/ms in case of eikonal model
        "sigma": Quantity(1, "mS/mm"),  # mS/mm
        "beta": Quantity(140, "1/mm"),
        "cm": Quantity(0.01, "uF/mm^2"),  # uF/mm^2
        "lambda": Quantity(0.2, "dimensionless"),
        "pmjrestype": Quantity(1),
        "pmjres": Quantity(0.001, "1/mS"),  # 1/mS
    },
}
