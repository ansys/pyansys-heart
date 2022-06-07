"""This uses dynalib to create some commonly used material cards and their 
default values. Uses inheritence.
E.g.:
Mat295
Mat077
MatNull
"""
import pandas as pd
import numpy as np

from ansys.heart.custom_logging import logger

from ansys.dyna.keywords import keywords

# import custom keywords. Overwrites classes defined in "keywords"
from ansys.heart.writer.custom_dynalib_keywords._custom_mat_077h import (
    Mat077H as _custom_Mat077H,
)


class MaterialCap(keywords.MatNull):
    """Material of the closing cap/valves

    Parameters
    ----------
    keywords : keywords.MatNull
        Inherits from Null type material
    """

    def __init__(self, mid: int = 1):
        super().__init__(mid=mid, ro=1.04e-6)


class MaterialAtrium(_custom_Mat077H):
    """Material for the atrium

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
        c10: float = 7.46,
    ):

        super().__init__(mid=mid, ro=rho, pr=poisson_ratio, n=0, c10=c10)
        return


class MaterialHGOMyocardium(keywords.Mat295):
    """Material for the myocardium. Either with active,
    passive, or both modules activated

    Parameters
    ----------
    keywords : keywords.Mat295
        Base keyword from which to derive this keyword

    Note
    ------
    How to properly parse material settings
    """

    def __init__(
        self,
        mid: int = 1,
        rho: float = 1e-6,
        aopt: float = 2.0,
        add_anisotropy: bool = False,
        add_active: bool = False,
    ):
        # some defaults for the isotropic module
        itype: int = -3
        beta: float = 0.0
        nu: float = 0.499
        k1: float = 2.36
        k2: float = 1.75

        # initialize isotropic module
        super().__init__(
            mid=mid,
            rho=rho,
            aopt=aopt,
            itype=itype,
            beta=beta,
            nu=nu,
            k1=k1,
            k2=k2,
        )

        if not add_anisotropy and add_active:
            raise ValueError(
                "Cannot add active module if not adding anisotropy"
            )

        # add defaults for anisotropic module
        if add_anisotropy:
            self.atype = -1
            self.intype = 0
            self.nf = 1  # number of fibers

            # activates fiber type of the anisotropic module
            ftype_aniso = 1
            self.ftype = ftype_aniso

            # prepare pandas dataframe:
            columns = self.anisotropic_settings.columns
            theta_aniso = 0.0
            a_aniso = 0.08
            b_aniso = 0.76

            fcid_aniso = 0
            k1_aniso = 0.49
            k2_aniso = 9.01
            aniso_data = np.array(
                [
                    [
                        theta_aniso,
                        a_aniso,
                        b_aniso,
                        ftype_aniso,
                        fcid_aniso,
                        k1_aniso,
                        k2_aniso,
                    ]
                ]
            )

            df = pd.DataFrame(data=aniso_data, columns=columns[0:7])

            # once set cannot change?
            self.anisotropic_settings = df

        if add_active:
            self._add_active_module(model_type="GuccioneWaldmanMcCulloch")

            self.acdir = 1

            self.acid = 15  # could be changed
            self.acthr = 0.1
            self.sf = 1.0
            self.ss = 0.02
            self.sn = 0.02

            # adds active fibers

    def _add_active_module(self, model_type: str = "GuccioneWaldmanMcCulloch"):
        """Adds the GuccioneWaldmanMcCulloch active module

        Parameters
        ----------
        model_type : str, optional
            Active model type, by default "GuccioneWaldmanMcCulloch"
        """
        if model_type == "GuccioneWaldmanMcCulloch":
            # activates actype 2
            self.actype = 2
            # parameters of actype 2
            self.ca2ionm = 4.35
            self.n = 2
            self.taumax = 125
            self.stf = 0.0
            self.b = 4.75
            self.l0 = 1.58
            self.l = 1.78
            self.eta = 1.45
        else:
            raise ValueError("Model type: %s not supported" % model_type)

        return


def active_curve(
    curve_type: str = "Strocchi2020",
    endtime: float = 15,
    timestep: float = 1e-2,
):
    """Computes various (normalized) curves used for the active module

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
        logger.warning("End-time set to 1000 ms: Strocchi uses 800 ms")
        t_end = 1000  # note in Strocchi this seems to be 800 ms actually
        t = np.arange(0, t_end, timestep * 1e3)
        # Tpeak = 125
        tau_r = 130
        tau_d = 100
        t_dur = 550

        # this is normalized in y: that is Tpeak is not used
        active_stress = (np.tanh(t / tau_r)) ** 2 * (
            np.tanh((t_dur - t) / tau_d)
        ) ** 2
        active_stress[t > t_dur] = 0

        # repeat dataset nCycles times:
        t = t / 1000  # to seconds
        t_end = t_end / 1000

        # number of cycles to return
        nCycles = int(np.ceil(endtime / t_end))

        time_array = t  # time array
        active_stress_array = active_stress  # active stress array
        for ii in range(1, nCycles):
            time_array = np.append(time_array, t + ii * t_end)
            active_stress_array = np.append(active_stress_array, active_stress)

        # import matplotlib

        # matplotlib.use(
        #     "Qt5Agg"
        # )  # note: solves unresponsive plotwindow in interactive mode
        # from matplotlib import pyplot as plt
        # plt.plot(T, TA, '.-')
        # plt.show()

    return time_array, active_stress_array


if __name__ == "__main__":

    active_curve()

    kw = MaterialCap()

    # Mat295 with anisotropy without active modue
    kw = MaterialHGOMyocardium(add_anisotropy=True, add_active=False)
    # Mat295 with anisotropy and active module
    kw = MaterialHGOMyocardium(add_anisotropy=True, add_active=True)

    kw = MaterialAtrium()

    print()
