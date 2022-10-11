"""
Methods for reading LS-DYNA binout file.

Note: depend on qd, can be replaced other Ansys's package like by pydpf, lsreader?
"""
import numpy as np

try:
    from qd.cae.Binout import Binout
except ImportError:
    raise ImportError("qd not found.")


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
