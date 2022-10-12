"""
Methods for reading LS-DYNA binout file.

Note: depend on qd, can be replaced other Ansys's package like by pydpf, lsreader?
"""
import numpy as np

try:
    from qd.cae.Binout import Binout
except ImportError:
    raise ImportError("qd not found.")


class Elout:
    """Object of element output."""

    def __init__(self, fn):
        """
        Init Elout object.

        Parameters
        ----------
        fn: binout file from LS-DYNA
        """
        self.binout = Binout(fn)
        self.time = self.binout.read("elout", "solid_hist", "time")
        self.ids = self.binout.read("elout", "solid_hist", "ids")[0]
        self.solid_nb = len(self.ids)

    def get_stress(self):
        """
        Get stress.

        Returns
        -------
        numpy array of time * element * 6 stress.
        stress-xx   stress-yy   stress-zz   stress-xy   stress-yz   stress-zx
        """
        data = self.binout.read("elout", "solid_hist", "data").reshape(len(self.time), -1, 9)
        if not np.all(data[:, :, 0] == 4) or data.shape[1] != self.solid_nb:
            raise Exception("Cannot read stress correctly.")
        else:
            return data[:, :, 1:7]

    def get_history_variable(self):
        """
        Get history variables.

        Returns
        -------
        numpy array of time * element * 27 history variables.
        """
        data = self.binout.read("elout", "solid_hist", "hist").reshape(len(self.time), -1, 28)
        if not np.all(data[:, :, 0] == 4) or data.shape[1] != self.solid_nb:
            raise Exception("Cannot read hv correctly.")
        else:
            return data[:, :, 1:]


def get_coordinates(binfile: str):
    """
    Read binout file to extract node coordinate.

    Parameters
    ----------
    binfile

    Returns
    -------
    time
    coordinates
    """
    data = Binout(binfile)
    time = data.read("nodout", "time")
    x_coord = data.read("nodout", "x_coordinate")
    y_coord = data.read("nodout", "y_coordinate")
    z_coord = data.read("nodout", "z_coordinate")
    coords = np.stack((x_coord, y_coord, z_coord), axis=2)
    return time, coords


def get_icvout(binfile: str):
    """
    Read binout file to extract control volume information.

    Parameters
    ----------
    binfile

    Returns
    -------
    time, pressure, volume, flow
    """
    binout = Binout(binfile)
    time = binout.read("icvout", "time")  # s
    pressure = binout.read("icvout", "ICV_Pressure")
    volume = binout.read("icvout", "ICV_Volume")
    flow = binout.read("icvout", "ICVI_flow_rate")
    return time, pressure, volume, flow


if __name__ == "__main__":
    pass
