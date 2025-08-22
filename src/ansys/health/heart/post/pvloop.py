# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Get pressure-volume (PV) loop from the LS-DYNA ICVOUT file."""

import os

import matplotlib.pyplot as plt
import numpy as np

from ansys.health.heart.post.dpf_utils import ICVoutReader

# NOTE: Assume control volume is defined in this order:
CV_NAME = ["LV", "RV", "LA", "RA"]


def write_pvloop_pngs(pressure: np.ndarray, volume: np.ndarray, save_to: str) -> None:
    """Write PV loop figures to a PNG file.

    Parameters
    ----------
    pressure : np.ndarray
        Pressure array.
    volume : np.ndarray
        Volume array.
    save_to : str
        Directory to save the file to.
    """
    n_cv = pressure.shape[0]
    for iframe in range(pressure.shape[1]):
        if iframe % 5 == 0:  # avoid too many figures
            fig: plt.Figure
            ax: plt.Axes
            fig, ax = plt.subplots(n_cv, figsize=(4, n_cv * 4))
            ax = [ax] if n_cv == 1 else ax

            for i in range(n_cv):  # loop for each cavity
                ax[i].set_xlim(volume[i].min() * 0.9, volume[i].max() * 1.1)
                ax[i].set_ylim(pressure[i].min() * 0.9, pressure[i].max() * 1.1)
                ax[i].set_xlabel("Volume (mL)")
                ax[i].set_ylabel("Pressure (kPa)")

                ax[i].plot(volume[i, :iframe], pressure[i, :iframe], label=CV_NAME[i])
                ax[i].legend()

            fig.savefig(os.path.join(save_to, f"pv_{iframe}.png"))

        plt.close()
    # NOTE: can combine into mp4 with:
    # ffmpeg -f image2 -i pv_%d.png output.mp4


def generate_pvloop(f: str, out_dir: str, t_to_keep: float = 800) -> None:
    """Generate PV loop figures from the ICVOUT file.

    Parameters
    ----------
    f : str
        Path to the binout file.
    out_dir : str
        Directory to save the file to.
    t_to_keep : float, optional
        Time to keep from the end. The default is ``800``, which is the last heart beat.
    """
    icvout = ICVoutReader(f)
    n_cv = len(icvout._icv_ids)
    t = icvout.get_time()
    pressure = icvout.get_pressure(1)
    volume = icvout.get_volume(1)

    for i in range(1, n_cv):
        pressure = np.vstack((pressure, icvout.get_pressure(1 + i)))
        volume = np.vstack((volume, icvout.get_volume(1 + i)))

    pressure = np.atleast_2d(pressure) * 1000  # to KPa
    volume = np.atleast_2d(volume) / 1000  # to mL

    np.savetxt(os.path.join(out_dir, "time_ms.txt"), t)
    np.savetxt(os.path.join(out_dir, "pressure_kPa.txt"), pressure.T)
    np.savetxt(os.path.join(out_dir, "volume_mL.txt"), volume.T)

    # last cycle
    ns = np.where(t >= t[-1] - t_to_keep)[0][0]
    pressure = pressure[:, ns:]
    volume = volume[:, ns:]

    #
    write_pvloop_pngs(pressure, volume, out_dir)
