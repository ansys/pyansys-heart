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

import os
import pathlib
from pathlib import Path

import matplotlib.pyplot as plt
import meshio
import numpy as np
import pyvista as pv
import vtk

from ansys.heart.postprocessor.dpf_utils import D3plotReader
from ansys.heart.preprocessor.mesh.vtkmethods import vtk_cutter, write_vtkdata_to_vtkfile

# TODO: replace by v0_2
from ansys.heart.preprocessor.models import HeartModel, LeftVentricle


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
            filename to be saved, default is None

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


class LVContourExporter:
    """Export Left ventricle surface and Post-process."""

    def __init__(self, d3plot_file: str, model: HeartModel):
        """
        Init LVContourExporter.

        Parameters
        ----------
        d3plot_file:
        model:
        """
        self.model = model
        self.work_dir = os.path.join(Path(d3plot_file).parent.absolute(), "post")
        self.data = D3plotReader(d3plot_file)
        self.nb_frame = len(self.data.time)

        self.out_folder = "lv_surface"
        os.makedirs(os.path.join(self.work_dir, self.out_folder), exist_ok=True)

        # TODO: get ID dynamically by part.pid or mid
        if isinstance(self.model, LeftVentricle):
            keep_ids = [1]
        else:
            keep_ids = [1, 3]

        self.data.export_vtk(
            os.path.join(self.work_dir, self.out_folder), only_surface=True, keep_mat_ids=keep_ids
        )

        self.lv_surfaces = []
        for i in range(self.nb_frame):
            abspath = os.path.join(self.work_dir, self.out_folder, f"model_{i}.vtk")
            self.lv_surfaces.append(pv.read(abspath))

        # get ID of mesh
        for ap in self.model.left_ventricle.apex_points:
            if ap.name == "apex epicardium":
                self.apex_id = ap.node_id
        for cap in self.model.left_ventricle.caps:
            if cap.name == "mitral-valve":
                self.mv_ids = cap.global_node_ids_edge

    def export_contour_to_vtk(self, folder, cutter) -> [vtk.vtkPolyData]:
        """
        Cut and save the contour in vtk.

        Parameters
        ----------
        cutter: dict contain 'center' and 'normal'
        folder: output folder

        Returns
        -------
        List of contour data
        """
        w_dir = os.path.join(self.work_dir, folder)
        os.makedirs(w_dir, exist_ok=True)
        cut_surfaces = []
        for id, surface in enumerate(self.lv_surfaces):
            res = vtk_cutter(surface, cutter)
            cut_surfaces.append(res)

        for ic, contour in enumerate(cut_surfaces):
            write_vtkdata_to_vtkfile(contour, os.path.join(w_dir, f"contour_{ic}.vtk"))
        return cut_surfaces

    def _compute_lvls(self):
        """
        Compute left ventricle long axis shortening.

        Returns
        -------
        Coordinates of mitral valve center and apex.
        """
        ic = self.data.get_initial_coordinates()
        dsp = self.data.get_displacement()

        # get coordinates
        apex = np.zeros((self.nb_frame, 3))
        mv_center = np.zeros((self.nb_frame, 3))

        for iframe, d in enumerate(dsp):
            apex[iframe, :] = ic[self.apex_id] + d[self.apex_id]
            mv_center[iframe, :] = np.mean(ic[self.mv_ids] + d[self.mv_ids], axis=0)

        return mv_center, apex

    def plot_lvls(self):
        """
        Get lvls figure.

        Returns
        -------
        fig: figure handle
        """
        mv_center, apex = self._compute_lvls()

        fig, ax = plt.subplots(1)
        ax.plot(np.linalg.norm(apex - mv_center, axis=1))
        ax.set_xlabel("Time")
        ax.set_ylabel("(mm)")
        ax.set_xticklabels([])

        return fig

    def export_lvls_to_vtk(self, folder="lvls_vtk"):
        """
        Export lvls as vtk Polyline.

        Parameters
        ----------
        folder: out folder
        """
        mv_center, apex = self._compute_lvls()

        os.makedirs(os.path.join(self.work_dir, folder), exist_ok=True)

        for iframe in range(self.nb_frame):
            meshio.write_points_cells(
                os.path.join(self.work_dir, folder, f"lvls_{iframe}.vtk"),
                np.array([mv_center[iframe], apex[iframe]]),
                [("line", np.array([[0, 1]]))],
            )
        # export lvls as csv file
        dst = np.linalg.norm(apex - mv_center, axis=1)
        time = self.data.time
        np.savetxt(
            pathlib.Path(self.work_dir) / folder / "lvls.csv",
            np.array([time, dst]).T,
            header="time,mv_apex dst",
            delimiter=",",
            comments="",
        )

        return dst
