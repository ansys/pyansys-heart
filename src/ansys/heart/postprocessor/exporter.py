# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
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

"""
Export information from d3plot.

Mostly related to the motion.
"""

import numpy as np
import pyvista as pv

from ansys.heart.postprocessor.dpf_utils import D3plotReader


class D3plotToVTKExporter:
    """Read d3plot and save deformed mesh."""

    def __init__(self, d3plot_file: str, t_to_keep: float = 10.0e10) -> None:
        """Init.

        Parameters
        ----------
        d3plot_file : str
            d3plot file path
        t_to_keep : float, optional
            time to be converted, by default 10.0e10
        """
        self.data = D3plotReader(d3plot_file)
        self.save_time = self.data.time[self.data.time >= self.data.time[-1] - t_to_keep]

    def convert_to_pvgrid_at_t(self, time: float, fname: str = None) -> pv.UnstructuredGrid:
        """Convert d3plot data into pyvista UnstructuredGrid.

        Parameters
        ----------
        time : float
            time to convert
        fname : str
            filename to be save save data, default is None

        Returns
        -------
        pv.UnstructuredGrid
            result in pyvista object
        """
        mesh = self.data.meshgrid.copy()
        i_frame = np.where(self.data.time == time)[0][0]
        dsp = self.data.get_displacement_at(time=time)
        mesh.points += dsp

        mesh.field_data["time"] = time
        mesh.cell_data["material_ids"] = self.data.get_material_ids()
        mesh.point_data["displacement"] = dsp

        tetra_ids = np.where(mesh.celltypes == 10)[0]

        mesh.cell_data["his16-18(fiber stretch)"] = np.empty((mesh.n_cells, 3))
        mesh.cell_data["his16-18(fiber stretch)"][tetra_ids] = self.data.get_history_variable(
            [15, 16, 17], at_step=i_frame
        ).T

        mesh.cell_data["his22-24(active stress)"] = np.empty((mesh.n_cells, 3))
        mesh.cell_data["his22-24(active stress)"][tetra_ids] = self.data.get_history_variable(
            [21, 22, 23], at_step=i_frame
        ).T

        mesh.cell_data["his25(ca2+)"] = np.empty(mesh.n_cells)
        mesh.cell_data["his25(ca2+)"][tetra_ids] = self.data.get_history_variable(
            [24], at_step=i_frame
        ).ravel()

        if fname is not None:
            mesh.save(fname)
        # NOTE: the returned pv_object seems corrupted, I suspect it's a bug of pyvista
        return mesh
