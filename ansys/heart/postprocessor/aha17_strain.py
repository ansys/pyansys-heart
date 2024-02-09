"""Compute myocardial strain."""
import pathlib

from ansys.heart import LOG as LOGGER
from ansys.heart.postprocessor.dpf_utils import D3plotReader
from ansys.heart.preprocessor.models.v0_1.models import HeartModel
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np


class AhaStrainCalculator:
    """Compute Longitudinal, Radial, Circumferential strain for left ventricle."""

    def __init__(self, model: HeartModel, d3plot_file):
        """
        Init.

        Parameters
        ----------
        model: HeartModel
            Heart model object.
        d3plot_file: Path.Path
            d3plot header file path.
        """
        self.model = model
        self._aha_elements = np.where(~np.isnan(model.aha_ids))[0]

        self.d3plot = D3plotReader(d3plot_file)

    def compute_aha_strain(self, out_dir: str, with_vtk: bool = False) -> np.ndarray:
        """Compute and save AHA 17 segment strain values.

        Parameters
        ----------
        out_dir : str
            Output directory where the segments are saved.
        with_vtk : bool, optional
            Flag indicating whether to save the VTK file, by default False

        Returns
        -------
        np.ndarray
            Average strain values for each of the 17 segments.
        """
        strain = np.zeros((len(self.d3plot.time), 1 + 17 * 3))

        if with_vtk:
            vtk_dir = out_dir
        else:
            vtk_dir = None

        header = "time"
        for aha in range(1, 18):
            for dir in ["L", "R", "C"]:
                header = ",".join([header, "AHA{0:d}_{1:s}".format(aha, dir)])

        for i, t in enumerate(self.d3plot.time):
            aha_lrc = self.compute_aha_strain_once(i, out_dir=vtk_dir)
            strain[i, 0] = t
            strain[i, 1:] = aha_lrc.ravel()

        np.savetxt(
            pathlib.Path(out_dir) / "AHAstrain.csv",
            strain,
            header=header,
            delimiter=",",
            comments="",
        )

        return strain

    def compute_aha_strain_once(self, frame: int = 0, out_dir: pathlib.Path = None) -> np.ndarray:
        """
        Export AHA strain and/or save vtk file for a given frame.

        Parameters
        ----------
        frame: int
            at this frame, by default 0.
        out_dir: pathlib.Path
            folder where vtk files are saved, by default not save.

        Returns
        -------
        np.ndarry
            AHA LRC strain matrix (17 * 3)
        """
        element_lrc, aha_lrc, element_lrc_averaged = self._compute_myocardial_strain(frame)

        if out_dir is not None:
            aha_model = self.model.mesh.extract_cells(self._aha_elements)
            aha_model.cell_data["AHA"] = self.model.aha_ids[self._aha_elements]

            init_coord = self.d3plot.get_initial_coordinates()[
                aha_model.point_data["vtkOriginalPointIds"]
            ]
            dsp = self.d3plot.get_displacement()[frame][aha_model.point_data["vtkOriginalPointIds"]]
            aha_model.points = init_coord + dsp

            aha_model.cell_data.set_vectors(element_lrc, "LRC strain")
            aha_model.cell_data.set_vectors(element_lrc_averaged, "LRC averaged strain")
            aha_model.save(pathlib.Path(out_dir) / "LRC_{0:d}.vtk".format(frame))

        return aha_lrc

    def _compute_myocardial_strain(self, at_frame, reference=None):
        """
        Compute left ventricle myocardial strain.

        Parameters
        ----------
        at_frame: int
        reference: not used

        Returns
        -------
        return1: [nelem * 3] elemental LRC strain
        return2: [17 * 3] AHA17 LRC strain
        return3: [nelem * 3] elemental LRC strain averaged from AHA17
        """
        if reference is not None:
            LOGGER.warning("Not implemented")
            exit()

        deformation_gradient = self.d3plot.get_history_variable(
            hv_index=list(range(9)), at_frame=at_frame
        ).T
        def_grad = deformation_gradient[self._aha_elements]

        strain = np.zeros((len(self._aha_elements), 3))
        averaged_strain = np.zeros((len(self._aha_elements), 3))
        aha_strain = np.zeros((17, 3))

        # model info
        e_l, e_r, e_c = self.model.compute_left_ventricle_element_cs()
        # todo: vectorization
        for i_ele in range(len(self._aha_elements)):
            if reference is not None:
                # todo
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
        aha17_label = self.model.aha_ids[self._aha_elements]

        for i in range(1, 18):
            # get index in strain table
            indices = np.where(aha17_label == i)[0]
            # average
            aha_strain[i - 1] = np.mean(strain[indices, :], axis=0)
            averaged_strain[indices] = aha_strain[i - 1]
        return strain, aha_strain, averaged_strain
        # return {"strain":strain,"aha_strain":aha_strain,"averaged_strain":averaged_strain}

    @staticmethod
    def bullseye_plot(ax, data, seg_bold=None, cmap=None, norm=None):
        """
        Bullseye representation for the left ventricle.

        Based on:
        https://matplotlib.org/stable/gallery/specialty_plots/leftventricle_bulleye.html#sphx-glr-gallery-specialty-plots-leftventricle-bulleye-py

        Parameters
        ----------
        ax : axes
        data : list of int and float
            The intensity values for each of the 17 segments
        seg_bold : list of int, optional
            A list with the segments to highlight
        cmap : ColorMap or None, optional
            Optional argument to set the desired colormap
        norm : Normalize or None, optional
            Optional argument to normalize data into the [0.0, 1.0] range

        Notes
        -----
        This function creates the 17 segment model for the left ventricle according
        to the American Heart Association (AHA) [1]_

        References
        ----------
        .. [1] M. D. Cerqueira, N. J. Weissman, V. Dilsizian, A. K. Jacobs,
            S. Kaul, W. K. Laskey, D. J. Pennell, J. A. Rumberger, T. Ryan,
            and M. S. Verani, "Standardized myocardial segmentation and
            nomenclature for tomographic imaging of the heart",
            Circulation, vol. 105, no. 4, pp. 539-542, 2002.
        """
        if seg_bold is None:
            seg_bold = []

        linewidth = 2
        data = np.ravel(data)

        if cmap is None:
            cmap = plt.cm.viridis

        if norm is None:
            norm = mpl.colors.Normalize(vmin=data.min(), vmax=data.max())

        theta = np.linspace(0, 2 * np.pi, 768)
        r = np.linspace(0.2, 1, 4)

        # Remove grid
        ax.grid(False)

        # Create the bound for the segment 17
        for i in range(r.shape[0]):
            ax.plot(theta, np.repeat(r[i], theta.shape), "-k", lw=linewidth)

        # Create the bounds for the segments 1-12
        for i in range(6):
            theta_i = np.deg2rad(i * 60)
            ax.plot([theta_i, theta_i], [r[1], 1], "-k", lw=linewidth)

        # Create the bounds for the segments 13-16
        for i in range(4):
            theta_i = np.deg2rad(i * 90 - 45)
            ax.plot([theta_i, theta_i], [r[0], r[1]], "-k", lw=linewidth)

        # Fill the segments 1-6
        r0 = r[2:4]
        r0 = np.repeat(r0[:, np.newaxis], 128, axis=1).T
        for i in range(6):
            # First segment start at 60 degrees
            theta0 = theta[i * 128 : i * 128 + 128] + np.deg2rad(60)
            theta0 = np.repeat(theta0[:, np.newaxis], 2, axis=1)
            z = np.ones((128, 2)) * data[i]
            # ax.pcolormesh(theta0, r0, z, cmap=cmap, norm=norm, shading="auto")

            ax.text(theta0.mean(), r0.mean(), "{0:.2f}".format(data[i]), fontsize=12)

            if i + 1 in seg_bold:
                ax.plot(theta0, r0, "-k", lw=linewidth + 2)
                ax.plot(theta0[0], [r[2], r[3]], "-k", lw=linewidth + 1)
                ax.plot(theta0[-1], [r[2], r[3]], "-k", lw=linewidth + 1)

        # Fill the segments 7-12
        r0 = r[1:3]
        r0 = np.repeat(r0[:, np.newaxis], 128, axis=1).T
        for i in range(6):
            # First segment start at 60 degrees
            theta0 = theta[i * 128 : i * 128 + 128] + np.deg2rad(60)
            theta0 = np.repeat(theta0[:, np.newaxis], 2, axis=1)
            z = np.ones((128, 2)) * data[i + 6]
            # ax.pcolormesh(theta0, r0, z, cmap=cmap, norm=norm, shading="auto")

            ax.text(theta0.mean(), r0.mean(), "{0:.2f}".format(data[i + 6]), fontsize=12)

            if i + 7 in seg_bold:
                ax.plot(theta0, r0, "-k", lw=linewidth + 2)
                ax.plot(theta0[0], [r[1], r[2]], "-k", lw=linewidth + 1)
                ax.plot(theta0[-1], [r[1], r[2]], "-k", lw=linewidth + 1)

        # Fill the segments 13-16
        r0 = r[0:2]
        r0 = np.repeat(r0[:, np.newaxis], 192, axis=1).T
        for i in range(4):
            # First segment start at 45 degrees
            theta0 = theta[i * 192 : i * 192 + 192] + np.deg2rad(45)
            theta0 = np.repeat(theta0[:, np.newaxis], 2, axis=1)
            z = np.ones((192, 2)) * data[i + 12]
            # ax.pcolormesh(theta0, r0, z, cmap=cmap, norm=norm, shading="auto")

            ax.text(theta0.mean(), r0.mean(), "{0:.2f}".format(data[i + 12]), fontsize=12)

            if i + 13 in seg_bold:
                ax.plot(theta0, r0, "-k", lw=linewidth + 2)
                ax.plot(theta0[0], [r[0], r[1]], "-k", lw=linewidth + 1)
                ax.plot(theta0[-1], [r[0], r[1]], "-k", lw=linewidth + 1)

        # Fill the segments 17
        if data.size == 17:
            r0 = np.array([0, r[0]])
            r0 = np.repeat(r0[:, np.newaxis], theta.size, axis=1).T
            theta0 = np.repeat(theta[:, np.newaxis], 2, axis=1)
            z = np.ones((theta.size, 2)) * data[16]
            # ax.pcolormesh(theta0, r0, z, cmap=cmap, norm=norm, shading="auto")

            ax.text(theta0.mean(), r0.mean(), "{0:.2f}".format(data[16]), fontsize=12)

            if 17 in seg_bold:
                ax.plot(theta0, r0, "-k", lw=linewidth + 2)

        ax.set_ylim([0, 1])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
