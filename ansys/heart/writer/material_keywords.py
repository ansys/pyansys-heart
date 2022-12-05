"""
Use dynalib to create some commonly used material cards and their default values.

Note
----
E.g.:
Mat295
Mat077
MatNull

"""
from ansys.dyna.keywords import keywords
from ansys.heart.custom_logging import LOGGER

# import custom keywords in separate namespace
from ansys.heart.writer import custom_dynalib_keywords as custom_keywords
import numpy as np
import pandas as pd
import pkg_resources


class MaterialCap(keywords.MatNull):
    """Material of the closing cap/valves.

    Parameters
    ----------
    keywords : keywords.MatNull
        Inherits from Null type material
    """

    def __init__(self, mid: int = 1):
        super().__init__(mid=mid, ro=1.04e-6)


class MaterialAtrium(custom_keywords.Mat077H):
    """Material for the atrium.

    Parameters
    ----------
    Mat077H : Parent class
        Parent class from which this material is derived
    """

    def __init__(
        self,
        mid: int = 1,
        rho: float = 1e-6,
        poisson_ratio: float = 0.499,
        c10: float = 17.46,
    ):
        super().__init__(mid=mid, ro=rho, pr=poisson_ratio, n=0, c10=c10)
        return


class MaterialHGOMyocardium(keywords.Mat295):
    """HGO Material model - derived from Mat295."""

    def __init__(self, mid: int = 1, iso_user=None, anisotropy_user=None, active_user=None):

        # Default parameters
        base = {"aopt": 2.0, "itype": -3, "beta": 0.0, "nu": 0.499}
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
            anisotropy = {
                "theta": 0,
                "a": 0.08,
                "b": 0.76,
                "fcid": 0,
                "ftype": 1,
                "atype": -1,
                "intype": 0,
                "nf": 1,
            }
            self.atype = anisotropy["atype"]
            self.intype = anisotropy["intype"]
            self.nf = anisotropy["nf"]
            self.ftype = anisotropy["ftype"]

            # Update parameters with user's input
            anisotropy.update(anisotropy_user)
            # transfer into keywords
            self.anisotropic_settings = pd.DataFrame([anisotropy])

        if active_user is not None:
            if active_user["actype"] == 1:
                active = {
                    # "actype": 1,
                    "acdir": 1,
                    "acthr": active_user["ca2ionm"] / 2,
                    "sf": 1.0,
                    "ss": 0.0,
                    "sn": 0.0,
                    "n": 2,
                    "stf": 0.0,
                    "b": 4.75,
                    "l0": 1.58,
                    "l": 1.85,
                    "dtmax": 150,  # ms
                    "mr": 1048.9,  # ms*um^-1
                    "tr": -1429,  # ms
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
        file_path = pkg_resources.resource_filename("ansys.heart.writer", "calcium_from_EP.txt")
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
    dct_aniso = {"k1": 1, "k2": 2}
    dct_active = {"acid": 15, "taumax": 125, "ca2ionm": 4.35}
    kw = MaterialHGOMyocardium(
        mid=1, iso_user=dct_iso, anisotropy_user=dct_aniso, active_user=dct_active
    )
    print(kw)

    dct_iso2 = {
        "rho": 1e-6,
        "itype": -1,
        "mu1": 34.9,
        "alpha1": 2,
        "Comment": "Should be equivalent with MAT_077_H",
    }
    kw = MaterialHGOMyocardium(mid=1, iso_user=dct_iso2)
    print(kw)
