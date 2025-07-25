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

"""Compute myocardial strain."""

import pathlib

import numpy as np
import pyvista as pv

from ansys.health.heart.models import BiVentricle, FourChamber, FullHeart, HeartModel, LeftVentricle
from ansys.health.heart.post.dpf_utils import D3plotReader
from ansys.health.heart.utils.landmark_utils import compute_aha17, compute_element_cs
from ansys.health.heart.utils.vtk_utils import find_corresponding_points, generate_thickness_lines


class AhaStrainCalculator:
    """Compute longitudinal, radial, and circumferential strain for the left ventricle."""

    def __init__(self, model: HeartModel, d3plot_file):
        """
        Initialize the AHA strain calculator.

        Parameters
        ----------
        model: HeartModel
            Heart model object.
        d3plot_file: Path.Path
            Path to the d3plot header file.
        """
        self.model = model

        self.aha_labels = compute_aha17(model, model.short_axis, model.l4cv_axis)
        self._aha_elements = np.where(~np.isnan(self.aha_labels))[0]

        self.d3plot = D3plotReader(d3plot_file)

    def _compute_thickness_lines(self, time_array: np.ndarray | list = None) -> list[pv.PolyData]:
        """Compute ventricular myocardium thickness.

        Parameters
        ----------
        time_array : np.ndarray | list, optional
           Time array to export. The d3plot time is exported by default.

        Returns
        -------
        list[pv.PolyData]
            Polydata that has lines from nodes on the endocardium to nodes on the epicardium.

        Notes
        -----
        Endocardium surfaces are supposed to be smooth.
        Artifact may occur on base (close to valves) region.
        """
        if time_array is None:
            time_array = self.d3plot.time

        surface_endo = self.model.left_ventricle.endocardium.copy()

        if isinstance(self.model, LeftVentricle):
            surface_epi = self.model.left_ventricle.epicardium.copy()
        elif isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
            sept_endo = self.model.right_ventricle.septum
            surface_epi = sept_endo.copy() + self.model.left_ventricle.epicardium.copy()

        lines_list = self._compute_thickness(time_array, surface_endo, surface_epi)

        # labeled as 1 for left ventricle
        for lines in lines_list:
            lines.cell_data["label"] = np.ones(lines.GetNumberOfCells())

        if isinstance(self.model, LeftVentricle):
            return lines_list

        elif isinstance(self.model, (BiVentricle, FourChamber, FullHeart)):
            # continue for right ventricle
            surface_endo = self.model.right_ventricle.endocardium.copy()
            surface_epi = self.model.right_ventricle.epicardium.copy()
            line_list2 = self._compute_thickness(time_array, surface_endo, surface_epi)

            # labeled as 2 for right ventricle
            for lines in line_list2:
                lines.cell_data["label"] = np.ones(lines.GetNumberOfCells()) * 2

            # merge polydata
            for i in range(len(lines_list)):
                lines_list[i] += line_list2[i]

            return lines_list

    def _compute_thickness(
        self, time_array: np.ndarray, surface_endo: pv.PolyData, surface_epi: pv.PolyData
    ) -> list[pv.PolyData]:
        """Compute thickness lines from endocardium to epicardium.

        Parameters
        ----------
        time_array : np.ndarray
            Time array to export.
        surface_endo : pv.PolyData
            Endocardium surface.
        surface_epi : pv.PolyData
            Epicardium surface.

        Returns
        -------
        list[pv.PolyData]
            Thickness lines.
        """
        res = []
        # assumes that corresponding points don't change in time
        pair = find_corresponding_points(surface_endo, surface_epi)
        for t in time_array:
            coordinates = (
                self.d3plot.model.results.coordinates.on_time_scoping(float(t)).eval()[0].data
            )
            # update surface coordinates
            surface_endo.points = coordinates[surface_endo["_global-point-ids"]]
            surface_epi.points = coordinates[surface_epi["_global-point-ids"]]

            thickness_lines = generate_thickness_lines(surface_endo, surface_epi, pair)
            thickness_lines.field_data["time"] = t

            res.append(thickness_lines)
        return res

    def compute_aha_strain(
        self, out_dir: str = None, write_vtk: bool = False, t_to_keep: float = 10e10
    ) -> np.ndarray:
        """Compute AHA 17 segment strain values from the deformation gradient.

        Parameters
        ----------
        out_dir : str, default: None
            Output folder.
        write_vtk : bool, default: False
            Whether to write to VTK files.
        t_to_keep : float, default: 10e10
            Time to stop.

        Returns
        -------
        np.ndarray
            Array of N_time * (1+17*3). Columns represent time and
            longitudinal, radial, and circumferential strain averaged of each segment.
        """
        save_time = self.d3plot.time[self.d3plot.time >= self.d3plot.time[-1] - t_to_keep]
        strain = np.zeros((len(save_time), 1 + 17 * 3))

        if write_vtk:
            vtk_dir = out_dir
        else:
            vtk_dir = None

        header = "time"
        for aha in range(1, 18):
            for dir in ["L", "R", "C"]:
                header = ",".join([header, "AHA{0:d}_{1:s}".format(aha, dir)])

        for i, t in enumerate(save_time):
            aha_lrc = self.compute_aha_strain_at(
                np.where(self.d3plot.time == t)[0][0], out_dir=vtk_dir
            )
            strain[i, 0] = t
            strain[i, 1:] = aha_lrc.ravel()

        if out_dir is not None:
            np.savetxt(
                pathlib.Path(out_dir) / "AHAstrain.csv",
                strain,
                header=header,
                delimiter=",",
                comments="",
            )

        return strain

    def compute_aha_strain_at(self, frame: int = 0, out_dir: pathlib.Path = None) -> np.ndarray:
        """
        Export AHA strain and/or save a VTK file for a given frame.

        Parameters
        ----------
        frame: int, default: 0
            Frame number to compute strain.
        out_dir: pathlib.Path, default: None
            Directory to save VTK file to. No VTK file is saved by default.

        Returns
        -------
        np.ndarry
            AHA LRC strain matrix (17 * 3).
        """
        element_lrc, aha_lrc, element_lrc_averaged = self._compute_myocardial_strain(frame)

        if out_dir is not None:
            aha_model = self.model.mesh.extract_cells(self._aha_elements)
            aha_model.cell_data["AHA"] = self.aha_labels[self._aha_elements]

            init_coord = self.d3plot.get_initial_coordinates()[
                aha_model.point_data["vtkOriginalPointIds"]
            ]

            dsp = self.d3plot.get_displacement_at(self.d3plot.time[frame])[
                aha_model.point_data["vtkOriginalPointIds"]
            ]
            aha_model.points = init_coord + dsp

            aha_model.cell_data.set_vectors(element_lrc, "LRC strain")
            aha_model.cell_data.set_vectors(element_lrc_averaged, "LRC averaged strain")
            aha_model.save(pathlib.Path(out_dir) / "LRC_{0:d}.vtk".format(frame))

        return aha_lrc

    def _compute_myocardial_strain(
        self, at_frame, reference=None
    ) -> tuple[list[float], list[float], list[float]]:
        """
        Compute left ventricle myocardial strain.

        Parameters
        ----------
        at_frame: int
        reference: not used

        Returns
        -------
        return1: [nelem * 3] elemental LRC strain.
        return2: [17 * 3] AHA17 LRC strain.
        return3: [nelem * 3] elemental LRC strain averaged from AHA17.
        """
        if reference is not None:
            raise NotImplementedError

        deformation_gradient = self.d3plot.get_history_variable(
            hv_index=list(range(9)), at_step=at_frame
        ).T
        def_grad = deformation_gradient[self._aha_elements]

        strain = np.zeros((len(self._aha_elements), 3))
        averaged_strain = np.zeros((len(self._aha_elements), 3))
        aha_strain = np.zeros((17, 3))

        # model info
        e_l, e_r, e_c = compute_element_cs(self.model, self.model.short_axis, self._aha_elements)
        for i_ele in range(len(self._aha_elements)):
            if reference is not None:
                pass
            else:
                right_cauchy_green = np.matmul(
                    def_grad[i_ele, :].reshape(3, 3),
                    def_grad[i_ele, :].reshape(3, 3).T,
                )

            # Green Lagrangian strain: E = 0.5*(lambda**2-1)
            # lambda = sqrt(e*right_cauchy_green*e)
            strain[i_ele, 0] = 0.5 * (
                np.matmul(np.matmul(e_l[i_ele].T, right_cauchy_green), e_l[i_ele]) - 1
            )
            strain[i_ele, 1] = 0.5 * (
                np.matmul(np.matmul(e_r[i_ele].T, right_cauchy_green), e_r[i_ele]) - 1
            )
            strain[i_ele, 2] = 0.5 * (
                np.matmul(np.matmul(e_c[i_ele].T, right_cauchy_green), e_c[i_ele]) - 1
            )

        # get aha17 label for left ventricle elements
        aha17_label = self.aha_labels[self._aha_elements]

        for i in range(1, 18):
            # get index in strain table
            indices = np.where(aha17_label == i)[0]
            # average
            aha_strain[i - 1] = np.mean(strain[indices, :], axis=0)
            averaged_strain[indices] = aha_strain[i - 1]
        return strain, aha_strain, averaged_strain

    @staticmethod
    def bullseye_17_segments(ax, data) -> None:
        """Bullseye representation of the 17 segments of the left ventricle.

        Parameters
        ----------
        ax : axes
        data : list of int and float
            Intensity values for each of the 17 segments.

        Notes
        -----
        This function creates the 17-segment model for the left ventricle according
        to the American Heart Association (AHA) [1]_.

        This method is modified from the ``Matplotlibs`` `bullseye <https://matplotlib.org/stable/gallery/specialty_plots/leftventricle_bulleye.html>`_
        example. Copyright |copy| 2012- Matplotlib Development Team; All Rights Reserved.
        Modifications consisted of removing colors and including the values for each segment.


        References
        ----------
        .. [1] M. D. Cerqueira, N. J. Weissman, V. Dilsizian, A. K. Jacobs,
            S. Kaul, W. K. Laskey, D. J. Pennell, J. A. Rumberger, T. Ryan,
            and M. S. Verani, "Standardized myocardial segmentation and
            nomenclature for tomographic imaging of the heart",
            Circulation, vol. 105, no. 4, pp. 539-542, 2002.
        """
        linewidth = 2

        data = np.ravel(data)

        if data.size != 17:
            raise ValueError(
                f"Data must be a 1D array of length 17. Received data of length {data.size}."
            )

        theta = np.linspace(0, 2 * np.pi, 768)
        radii = np.linspace(0.2, 1, 4)

        # Remove grid
        ax.grid(False)

        # Create the bounds for segment 17
        for ii in range(radii.shape[0]):
            ax.plot(theta, np.repeat(radii[ii], theta.shape), "-k", lw=linewidth)

        # Create the bounds for segments 1-12
        for ii in range(6):
            theta_i = np.deg2rad(ii * 60)
            ax.plot([theta_i, theta_i], [radii[1], 1], "-k", lw=linewidth)

        # Create the bounds for segments 13-16
        for ii in range(4):
            theta_i = np.deg2rad(ii * 90 - 45)
            ax.plot([theta_i, theta_i], [radii[0], radii[1]], "-k", lw=linewidth)

        # Fill segments 1-6
        r0 = radii[2:4]
        r0 = np.repeat(r0[:, np.newaxis], 128, axis=1).T
        for ii in range(6):
            # First segment starts at 60 degrees
            theta0 = theta[ii * 128 : ii * 128 + 128] + np.deg2rad(60)
            theta0 = np.repeat(theta0[:, np.newaxis], 2, axis=1)
            # Print value to segment.
            ax.text(theta0.mean(), r0.mean(), "{0:.2f}".format(data[ii]), fontsize=12)

        # Fill segments 7-12
        r0 = radii[1:3]
        r0 = np.repeat(r0[:, np.newaxis], 128, axis=1).T
        for ii in range(6):
            # First segment starts at 60 degrees
            theta0 = theta[ii * 128 : ii * 128 + 128] + np.deg2rad(60)
            theta0 = np.repeat(theta0[:, np.newaxis], 2, axis=1)
            # Print value to segment.
            ax.text(theta0.mean(), r0.mean(), "{0:.2f}".format(data[ii + 6]), fontsize=12)

        # Fill segments 13-16
        r0 = radii[0:2]
        r0 = np.repeat(r0[:, np.newaxis], 192, axis=1).T
        for ii in range(4):
            # First segment start at 45 degrees
            theta0 = theta[ii * 192 : ii * 192 + 192] + np.deg2rad(45)
            theta0 = np.repeat(theta0[:, np.newaxis], 2, axis=1)
            # Print value to segment.
            ax.text(theta0.mean(), r0.mean(), "{0:.2f}".format(data[ii + 12]), fontsize=12)

        # Fill segment 17
        r0 = np.array([0, radii[0]])
        r0 = np.repeat(r0[:, np.newaxis], theta.size, axis=1).T
        theta0 = np.repeat(theta[:, np.newaxis], 2, axis=1)
        # Print value to segment.
        ax.text(theta0.mean(), r0.mean(), "{0:.2f}".format(data[16]), fontsize=12)

        ax.set_ylim([0, 1])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
