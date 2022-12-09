"""
Export information from d3plot.

Mostly related to the motion.
Depend on Paraview macro to convert d3plot to vtk.
"""
import os
from pathlib import Path
import subprocess

from ansys.heart.preprocessor.mesh.vtkmethods import (
    read_vtk_polydata_file,
    vtk_cutter,
    write_vtkdata_to_vtkfile,
)
from ansys.heart.preprocessor.models import HeartModel
import matplotlib.pyplot as plt
import meshio
import numpy as np
import pkg_resources
import vtk
from vtk.util.numpy_support import vtk_to_numpy  # noqa

# Paraview path
# todo
pvpython = r"C:\Program Files\ParaView 5.9.0-Windows-Python3.8-msvc2017-64bit\bin\pvpython.exe"
# script to convert d3plot to vtk
pv_script = pkg_resources.resource_filename("ansys.heart.misc.paraview_macro", "d3plot_to_vtk.pvpy")


class D3plotExporter:
    """D3plotExporter."""

    def __init__(self, d3plot_file: str, model: HeartModel):
        """
        Init.

        Parameters
        ----------
        d3plot_file: absolute path of d3plot
        model: HeartModel
        """
        self.d3plot = d3plot_file
        self.work_dir = Path(d3plot_file).parent.absolute()
        self.model = model

    def convert_to_vtk(self, parts: [str], out_folder: str, only_surface=True):
        """
        Convert d3plot to vtk via a Paraveiw Macro.

        Parameters
        ----------
        parts: parts to be kept
        out_folder: output folder
        only_surface: if extract surface
        """
        os.makedirs(os.path.join(self.work_dir, out_folder), exist_ok=True)
        out_file = os.path.join(self.work_dir, out_folder, "model.vtk")

        if only_surface:
            option = "1"
        else:
            option = "0"

        subprocess.call(
            [
                pvpython,
                pv_script,
                os.path.join(self.work_dir, self.d3plot),
                out_file,
                option,
                ";".join(parts),
            ]
        )


class LVContourExporter(D3plotExporter):
    """Export Left ventricle surface and Post-process."""

    def __init__(self, d3plot_file: str, model: HeartModel):
        """
        Init.

        Parameters
        ----------
        d3plot_file:
        model:
        """
        super().__init__(d3plot_file, model)

        self.out_folder = "lv_surface"
        parts = ["Left ventricle", "Septum"]
        self.convert_to_vtk(parts, self.out_folder, only_surface=True)
        self.nb_frame = len(os.listdir(os.path.join(self.work_dir, self.out_folder)))

        self.lv_surfaces = []
        for i in range(self.nb_frame):
            abspath = os.path.join(self.work_dir, self.out_folder, f"model_{i}.vtk")
            self.lv_surfaces.append(read_vtk_polydata_file(abspath))

        # get ID of mesh
        for ap in self.model.left_ventricle.apex_points:
            if ap.name == "apex epicardium":
                self.apex_id = ap.node_id + 1
        for cap in self.model.left_ventricle.caps:
            if cap.name == "mitral-valve":
                self.mv_ids = cap.node_ids + 1

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
        contours = self._get_contours(cutter)
        for ic, contour in enumerate(contours):
            write_vtkdata_to_vtkfile(contour, os.path.join(w_dir, f"contour_{ic}.vtk"))
        return contours

    def _get_contours(self, cut_plane):
        """
        Cut.

        Parameters
        ----------
        cut_plane: dict

        Returns
        -------
        contours
        """
        cut_surfaces = []
        for id, surface in enumerate(self.lv_surfaces):
            res = vtk_cutter(surface, cut_plane)
            cut_surfaces.append(res)
        return cut_surfaces

    def compute_lvls(self):
        """
        Compute left ventricle long axis shortening.

        Returns
        -------
        Coordinates of mitral valve center and apex.
        """
        # get node ID of vtk
        dyna_ids = vtk_to_numpy(self.lv_surfaces[0].GetPointData().GetArray("UserID"))
        apex_id = np.where(dyna_ids == self.apex_id)[0][0]
        sorter = np.argsort(dyna_ids)
        mv_ids = sorter[np.searchsorted(dyna_ids, self.mv_ids, sorter=sorter)]

        # get coordinates
        apex = np.zeros((self.nb_frame, 3))
        mv_center = np.zeros((self.nb_frame, 3))

        for iframe, surface in enumerate(self.lv_surfaces):
            apex[iframe, :] = np.array(surface.GetPoint(apex_id))
            coord = np.zeros(3)
            for point_id in mv_ids:
                coord += np.array(surface.GetPoint(point_id)) / len(mv_ids)
            mv_center[iframe, :] = coord

        return mv_center, apex

    def plot_lvls(self):
        """
        Get lvls figure.

        Returns
        -------
        fig: figure handle
        """
        mv_center, apex = self.compute_lvls()

        fig, ax = plt.subplots(1)
        ax.plot(np.linalg.norm(apex - mv_center, axis=1))
        ax.set_xlabel("Time")
        ax.set_ylabel("(mm)")
        ax.set_xticklabels([])

        return fig

    def export_lvls_to_vtk(self, fodler="lvls_vtk"):
        """
        Export lvls as vtk Polyline.

        Parameters
        ----------
        fodler: out folder
        """
        mv_center, apex = self.compute_lvls()

        os.makedirs(os.path.join(self.work_dir, fodler), exist_ok=True)

        for iframe in range(self.nb_frame):
            meshio.write_points_cells(
                os.path.join(self.work_dir, fodler, f"lvls_{iframe}.vtk"),
                np.array([mv_center[iframe], apex[iframe]]),
                [("line", np.array([[0, 1]]))],
            )
