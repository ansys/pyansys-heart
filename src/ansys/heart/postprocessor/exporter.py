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
import vtk
from vtk.util.numpy_support import vtk_to_numpy  # noqa

from ansys.heart.postprocessor.dpf_utils import D3plotReader
from ansys.heart.preprocessor.mesh.vtkmethods import (
    read_vtk_polydata_file,
    vtk_cutter,
    write_vtkdata_to_vtkfile,
)
from ansys.heart.preprocessor.models.v0_1.models import HeartModel, LeftVentricle


class LVContourExporter:
    """Export Left ventricle surface and Post-process."""

    def __init__(self, d3plot_file: str, model: HeartModel):
        """
        Init.

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

        # todo get ID dynamically by part.pid or mid
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
            self.lv_surfaces.append(read_vtk_polydata_file(abspath))

        # get ID of mesh
        for ap in self.model.left_ventricle.apex_points:
            if ap.name == "apex epicardium":
                self.apex_id = ap.node_id
        for cap in self.model.left_ventricle.caps:
            if cap.name == "mitral-valve":
                self.mv_ids = cap.node_ids

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
