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

"""D3plot parser using Ansys-dpf."""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

from ansys.heart.core.models import HeartModel
from ansys.heart.postprocessor.dpf_utils import D3plotReader


class EPpostprocessor:
    """Postprocess Electrophysiology results."""

    def __init__(self, results_path: Path, model: HeartModel = None):
        """Postprocess EP results.

        Parameters
        ----------
        results_path : Path
            Path to results.
        model : HeartModel
            Heart model.
        """
        self.reader = D3plotReader(results_path)
        self.fields = None
        self.model = model

    def load_ep_fields(self):
        """Load all EP fields."""
        if self.fields is None:
            self.fields = self.reader.get_ep_fields()

    def get_activation_times(self, at_step: int = None):
        """Get activation times field."""
        step = (
            self.reader.model.metadata.time_freq_support.time_frequencies.scoping.ids[-1]
            if at_step is None
            else [at_step]
        )
        field = self.reader.get_ep_fields(at_step=step).get_field({"variable_id": 129})
        return field

    def get_transmembrane_potential(self, node_id=None, plot: bool = False):
        """Get transmembrane potential."""
        phi, times = self._get_ep_field(node_id=node_id, plot=plot, variable_id=126)
        return phi, times

    def get_extracellular_potential(self, node_id=None, plot: bool = False):
        """Get extracellular potential."""
        phi, times = self._get_ep_field(variable_id=127, node_id=node_id, plot=plot)
        return phi, times

    def get_intracellular_potential(self, node_id=None, plot: bool = False):
        """Get intracellular potential."""
        phi, times = self._get_ep_field(variable_id=128, node_id=node_id, plot=plot)
        return phi, times

    def get_calcium(self, node_id=None, plot: bool = False):
        """Get calcium concentration."""
        phi, times = self._get_ep_field(variable_id=130, node_id=node_id, plot=plot)
        return phi, times

    def _get_ep_field(self, variable_id: int, node_id=None, plot: bool = False):
        """Get EP field."""
        self.load_ep_fields()
        times = self.reader.time
        if node_id is None:
            nnodes = len(self.reader.meshgrid.points)
            node_id = np.int64(np.linspace(0, nnodes - 1, nnodes))
        phi = np.zeros((len(times), len(node_id)))

        for time_id in range(1, len(times) + 1):
            phi[time_id - 1, :] = self.fields.get_field(
                {"variable_id": variable_id, "time": time_id}
            ).data[node_id]
        if plot:
            plt.plot(times, phi, label="node 0")
            plt.xlabel("time (ms)")
            plt.ylabel("phi (mV)")
            plt.show(block=True)
        return phi, times

    def read_ep_nodout(self):
        """Read Electrophysiology results."""
        em_nodout_path = os.path.join(self.results_path, "em_nodout_EP_001.dat")
        with open(em_nodout_path, "r") as f:
            lines = f.readlines()

        times = []
        line_indices = []
        nodes_list = []

        # Get times
        for index, line in enumerate(lines):
            if " at time " in line:
                times.append(float((line.split())[-2]))
                line_indices.append(index)
        self.times = times

        # Get node ids
        nodes_list = map(
            lambda x: int(x.split()[0]),
            (lines[int(line_indices[0]) + 3 : int(line_indices[1]) - 3]),
        )
        node_ids = list(nodes_list)

        self.node_ids = np.array(node_ids)

        # Get node activation times
        act_t = map(
            lambda x: float(x.split()[8]),
            (lines[int(line_indices[-1]) + 3 : int(len(lines))]),
        )
        self.activation_time = np.array(list(act_t))
        self._assign_pointdata(pointdata=self.activation_time, node_ids=self.node_ids)

    def create_post_folder(self, path: Path = None):
        """Create Postprocessing folder."""
        if path is None:
            post_path = os.path.join(os.path.dirname(self.reader.ds.result_files[0]), "post")
        else:
            post_path = path
        path_exists = os.path.exists(post_path)
        if not path_exists:
            # Create a new directory because it does not exist
            os.makedirs(post_path)
        return post_path

    def animate_transmembrane(self):
        """Animate transmembrane potentials and export to vtk."""
        vm, times = self.get_transmembrane_potential()
        # Creating scene and loading the mesh
        post_path = self.create_post_folder()
        grid = self.reader.meshgrid.copy()
        p = pv.Plotter()
        p.add_mesh(grid, scalars=vm[0, :])
        p.show(interactive_update=True)

        for i in range(vm.shape[0]):
            grid.point_data["transemembrane_potential"] = vm[i, :]
            grid.save(post_path + "\\vm_" + str(i) + ".vtk")
            p.update_scalars(vm[i, :])
            p.update()

        return

    def export_transmembrane_to_vtk(self):
        """Export transmembrane potentials to vtk."""
        vm, times = self.get_transmembrane_potential()
        post_path = self.create_post_folder()
        grid = self.reader.meshgrid.copy()

        for ii in range(vm.shape[0]):
            # TODO: vtk is not optimal for scalar fields with
            # non moving meshes, consider using ROM format
            grid.point_data["transmembrane_potential"] = vm[ii, :]
            grid.save(os.path.join(post_path, "vm_{0}.vtk".format(ii)))
        return

    def compute_ECGs(self, electrodes: np.ndarray):  # noqa: N802
        """Compute ECGs."""
        grid = self.reader.meshgrid
        grid = grid.compute_cell_sizes(length=False, area=False, volume=True)
        cell_volumes = grid.cell_data["Volume"]
        centroids = grid.cell_centers()
        vm, times = self.get_transmembrane_potential()
        # NOTE: ecgs: Electrocardiograms
        ecgs = np.zeros([vm.shape[0], electrodes.shape[0]])

        for time_step in range(vm.shape[0]):
            grid.point_data["vmi"] = vm[time_step, :]
            grid = grid.compute_derivative(scalars="vmi")
            grid = grid.point_data_to_cell_data()

            for electrode_id in range(electrodes.shape[0]):
                electrode = electrodes[electrode_id, :]
                r_vector = centroids.points - electrode
                distances = np.linalg.norm(r_vector, axis=1)
                # TODO: add conductivity tensor in the calculation (necessary?)
                # TODO: add method to handle beam gradients as well
                integral = sum(
                    sum(
                        np.transpose(
                            r_vector
                            * grid.cell_data["gradient"]
                            / (np.power(distances[:, None], 3) * 4 * np.pi)
                        )
                    )
                    * cell_volumes
                )
                ecgs[time_step, electrode_id] = integral
        return ecgs, times

    def read_ECGs(self, path: Path):  # noqa: N802
        """Read ECG text file produced by LS-DYNA simulation."""
        data = np.loadtxt(path, skiprows=4)
        times = data[:, 0]
        ecgs = data[:, 1:11]
        return ecgs, times

    # TODO: @KarimElHouari can we avoid capital letters in variable names somehow?
    def compute_12_lead_ECGs(  # noqa: N802
        self,
        ECGs: np.ndarray,  # noqa: N803
        times: np.ndarray,
        plot: bool = True,
    ) -> np.ndarray:
        """Compute 12-Lead ECGs from 10 electrodes.

        Parameters
        ----------
        ECGs : np.ndarray
            mxn array containing ECGs, where m is the number of time steps
            and n the 10 electrodes in this order:
            "V1" "V2" "V3" "V4" "V5" "V6" "RA" "LA" "RL" "LL"
        plot : bool, optional
            plot option, by default True

        Returns
        -------
        np.ndarray
            12-Lead ECGs in this order:
            "I" "II" "III" "aVR" "aVL" "aVF" "V1" "V2" "V3" "V4" "V5" "V6"
        """
        right_arm = ECGs[:, 6]
        left_arm = ECGs[:, 7]
        left_leg = ECGs[:, 9]
        lead1 = left_arm - right_arm  # noqa: E741
        lead2 = left_leg - right_arm
        lead3 = left_leg - left_arm
        lead_avr = right_arm - (left_arm + left_leg) / 2
        lead_avl = left_arm - (left_leg + right_arm) / 2
        lead_avf = left_leg - (right_arm + left_arm) / 2
        Vwct = (left_arm + right_arm + left_leg) / 3  # noqa: N806
        lead_v1 = ECGs[:, 0] - Vwct
        lead_v2 = ECGs[:, 1] - Vwct
        lead_v3 = ECGs[:, 2] - Vwct
        lead_v4 = ECGs[:, 3] - Vwct
        lead_v5 = ECGs[:, 4] - Vwct
        lead_v6 = ECGs[:, 5] - Vwct
        right_arm = ECGs[:, 6] - Vwct
        left_arm = ECGs[:, 7] - Vwct
        # RL = ECGs[8, :] - Vwct
        left_leg = ECGs[:, 9] - Vwct
        ecg_12lead = np.vstack(
            (
                lead1,
                lead2,
                lead3,
                lead_avr,
                lead_avl,
                lead_avf,
                lead_v1,
                lead_v2,
                lead_v3,
                lead_v4,
                lead_v5,
                lead_v6,
            )
        )
        if plot:
            t = times
            fig, axes = plt.subplots(nrows=3, ncols=4, layout="tight")
            # Major ticks every 20, minor ticks every 5
            major_xticks = np.arange(0, int(max(times)), 200)
            minor_xticks = np.arange(0, int(max(times)), 40)
            for i, ax in enumerate(fig.axes):
                ax.set_xticks(major_xticks)
                ax.set_xticks(minor_xticks, minor=True)
                ax.set_yticks(major_xticks)
                ax.set_yticks(minor_xticks, minor=True)
                ax.grid(which="both")
            axes[0, 0].plot(t, lead1)
            axes[0, 0].set_ylabel("I")
            axes[1, 0].plot(t, lead2)
            axes[1, 0].set_ylabel("II")
            axes[2, 0].plot(t, lead3)
            axes[2, 0].set_ylabel("III")
            axes[0, 1].plot(t, lead_avr)
            axes[0, 1].set_ylabel("aVR")
            axes[1, 1].plot(t, lead_avl)
            axes[1, 1].set_ylabel("aVL")
            axes[2, 1].plot(t, lead_avf)
            axes[2, 1].set_ylabel("aVF")
            axes[0, 2].plot(t, lead_v1)
            axes[0, 2].set_ylabel("V1")
            axes[1, 2].plot(t, lead_v2)
            axes[1, 2].set_ylabel("V2")
            axes[2, 2].plot(t, lead_v3)
            axes[2, 2].set_ylabel("V3")
            axes[0, 3].plot(t, lead_v4)
            axes[0, 3].set_ylabel("V4")
            axes[1, 3].plot(t, lead_v5)
            axes[1, 3].set_ylabel("V5")
            axes[2, 3].plot(t, lead_v6)
            axes[2, 3].set_ylabel("V6")
            plt.setp(plt.gcf().get_axes(), xticks=[0, 200, 400, 600, 800], yticks=[])
            # fig.add_subplot(111, frameon=False)
            # plt.tick_params(
            #     labelcolor="none", which="both", top=False, bottom=False, left=False, right=False
            # )
            # plt.xlabel("time (ms)")
            post_path = self.create_post_folder()
            filename = os.path.join(post_path, "12LeadECGs.png")
            plt.savefig(fname=filename, format="png")
            plt.show(block=True)
        return ecg_12lead

    def _assign_pointdata(self, pointdata: np.ndarray, node_ids: np.ndarray):
        """Assign point data to mesh."""
        result = np.zeros(self.mesh.n_points)
        result[node_ids - 1] = pointdata
        self.mesh.point_data["activation_time"] = result
