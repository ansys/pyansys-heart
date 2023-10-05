"""
Use dynalib to create some commonly used material cards and their default values.

Note
----
E.g.:
Mat295
Mat077
MatNull

"""
import logging

from ansys.dyna.keywords import keywords

LOGGER = logging.getLogger("pyheart_global.writer")

# from importlib.resources import files
from importlib.resources import path as resource_path

# import custom keywords in separate namespace
from ansys.heart.writer import custom_dynalib_keywords as custom_keywords
import numpy as np
import pandas as pd


class MaterialCap(keywords.MatNull):
    """Material of the closing cap/valves.

    Parameters
    ----------
    keywords : keywords.MatNull
        Inherits from Null type material
    """

    def __init__(self, mid: int = 1):
        super().__init__(mid=mid, ro=1.04e-6)


class MaterialNeoHook(custom_keywords.Mat077H):
    """Material for the atrium.

    Parameters
    ----------
    Mat077H : Parent class
        Parent class from which this material is derived
    """

    def __init__(
        self,
        mid: int,
        rho: float,
        c10: float,
        poisson_ratio: float = 0.499,
    ):
        super().__init__(mid=mid, ro=rho, pr=poisson_ratio, n=0, c10=c10)
        return


class MaterialHGOMyocardium(keywords.Mat295):
    """HGO Material model - derived from Mat295."""

    def __init__(self, mid: int = 1, iso_user=None, anisotropy_user=None, active_user=None):
        # Default parameters
        base = {"aopt": 2.0, "itype": -3, "beta": 0.0}
        # Update parameters with user's input
        base.update(iso_user)

        # initialize isotropic module
        super().__init__(
            mid=mid,
        )
        for k, v in base.items():
            setattr(self, k, v)
        if anisotropy_user is not None:
            # Default parameters
            self.atype = -1
            self.intype = 0
            self.nf = 1
            self.ftype = 1

            common_sheet_fiber = {
                "a": 0.0,
                "b": 1.0,
                "fcid": 0,
                "ftype": self.ftype,
            }

            # change type if key exist
            if "k1s" in anisotropy_user.keys():
                self.nf = 2
                if "k1fs" in anisotropy_user.keys():
                    self.intype = 1

            fiber = {
                "theta": 0,
                "k1": anisotropy_user["k1f"],
                "k2": anisotropy_user["k2f"],
            }
            fiber.update(common_sheet_fiber)

            if self.nf == 1:
                self.anisotropic_settings = pd.DataFrame([fiber])
            elif self.nf == 2:
                sheet = {
                    "theta": 90,
                    "k1": anisotropy_user["k1s"],
                    "k2": anisotropy_user["k2s"],
                }
                sheet.update(common_sheet_fiber)
                self.anisotropic_settings = pd.DataFrame([fiber, sheet])

            if self.intype == 1:
                self.coupling_k1 = anisotropy_user["k1fs"]
                self.coupling_k2 = anisotropy_user["k2fs"]

            # set dummy a/d vectors
            # these values are indeed replaced by *ELEMENT_SOLID_ORTHO or
            # by dynain.lsda file if they exist
            self.a1 = 1.0
            self.a2 = 0.0
            self.a3 = 0.0

            self.d1 = 0.0
            self.d2 = 1.0
            self.d3 = 0.0

        if active_user is not None:
            if active_user["actype"] == 1:
                active = {
                    "acdir": 1,  # active along first/fiber direction
                    "acthr": active_user["ca2ionm"] / 2,
                    "sf": 1.0,
                    "sn": 0.0,
                    "n": 2,
                    "stf": 0.0,
                    "b": 4.75,
                    "l0": 1.58,
                    "l": 1.85,
                    "mr": 1048.9,  # ms*um^-1
                }
            elif active_user["actype"] == 2:
                # Default parameters
                active = {
                    # "actype": 2,
                    "acdir": 1,
                    "acthr": 0.1,
                    "sf": 1.0,
                    "ss": 0.02,
                    "sn": 0.02,
                    # "ca2ionm": 4.35,
                    "n": 2,
                    "stf": 0.0,
                    "b": 4.75,
                    "l0": 1.58,
                    "l": 1.78,
                    "eta": 1.45,
                }
            # Update parameters with user's input
            active.update(active_user)
            # transfer into keywords
            for k, v in active.items():
                setattr(self, k, v)


def active_curve(
    curve_type: str = "Strocchi2020",
    endtime: float = 15,
    timestep: float = 1e-2,
):
    """Compute various (normalized) curves used for the active module.

    Parameters
    ----------
    curve_name : str
        Type of curve to compute
    """
    # time array
    # T = np.arange( 0, endtime, timestep )
    # NOTE: needs cleaning up
    if curve_type == "Strocchi2020":
        # parameters used in Strocchi:

        # NOTE: in milliseconds
        LOGGER.warning("End-time set to 1000 ms: Strocchi uses 800 ms")
        t_end = 1000  # note in Strocchi this seems to be 800 ms actually
        t = np.arange(0, t_end, timestep * 1e3)
        # Tpeak = 125
        tau_r = 130
        tau_d = 100
        t_dur = 550

        # this is normalized in y: that is Tpeak is not used
        active_stress = (np.tanh(t / tau_r)) ** 2 * (np.tanh((t_dur - t) / tau_d)) ** 2
        active_stress[t > t_dur] = 0

        # repeat dataset nCycles times:
        t = t / 1000  # to seconds
        t_end = t_end / 1000

        # number of cycles to return
        nCycles = int(np.ceil(endtime / t_end))

        time_array = t  # time array
        calcium_array = active_stress  # active stress array
        for ii in range(1, nCycles):
            time_array = np.append(time_array, t + ii * t_end)
            calcium_array = np.append(calcium_array, active_stress)

    # used for generating multi beats with model actype 1
    elif curve_type == "constant":
        nb_beats = 10
        period = 1.0  # in second

        # define shape pattern
        value = np.array([0, 1, 1])
        time = np.array([0, 0.001 * period, 0.9 * period])

        # repeat for every period
        time_array = time
        for i in range(1, nb_beats):
            time_array = np.append(time_array, time + period * i)
        calcium_array = np.tile(value, nb_beats)

        # append last point
        time_array = np.append(time_array, period * nb_beats)
        calcium_array = np.append(calcium_array, 0.0)

    elif curve_type == "TrueCalcium":
        file_path = resource_path("ansys.heart.writer", "calcium_from_EP.txt").__enter__()
        a = np.loadtxt(file_path)
        time_array = a[:, 0] / 1000
        calcium_array = a[:, 1]

    # import matplotlib.pyplot as plt
    # plt.plot(time_array, calcium_array)
    # plt.show()

    return time_array, calcium_array


if __name__ == "__main__":
    active_curve(curve_type="TrueCalcium")

    kw = MaterialCap()
    # test
    dct_iso = {"rho": 1, "k1": 1, "k2": 1}
    dct_aniso = {"k1f": 1, "k2f": 2}
    dct_aniso.update({"k1s": 1, "k2s": 2})
    # dct_aniso.update({"k1fs": 1, "k2fs": 2})
    dct_active = {"actype": 1, "acid": 15, "taumax": 125, "ca2ionm": 4.35}
    kw = MaterialHGOMyocardium(
        mid=1, iso_user=dct_iso, anisotropy_user=dct_aniso, active_user=dct_active
    )
    print(kw)
