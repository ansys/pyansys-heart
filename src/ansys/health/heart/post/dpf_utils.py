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

"""D3plot parser using Ansys DPF."""

import os
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

from ansys.dpf import core as dpf
from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.exceptions import (
    MissingEnvironmentVariableError,
    SupportedDPFServerNotFoundError,
)
from ansys.health.heart.models import HeartModel

_SUPPORTED_DPF_SERVERS = ["2025.2", "2025.2rc0", "2024.1", "2024.1rc1", "2024.2rc0"]
"""List of supported DPF servers."""
#! NOTE:
#! 2024.1rc0: not supported due to missing ::tf operator
#! => 2024.2rc1: not supported due to bug in dpf server when reading d3plots mixed with EM results


def _check_accept_dpf():
    if not os.getenv("ANSYS_DPF_ACCEPT_LA", None) == "Y":
        LOGGER.error(
            """DPF requires you to accept the license agreement.
            Set the environment variable "ANSYS_DPF_ACCEPT_LA" to "Y"."""
        )
        raise MissingEnvironmentVariableError("ANSYS_DPF_ACCEPT_LA not set to 'Y'.")
    return


def _get_dpf_server():
    """Get the DPF server."""
    server = None

    # sort available servers from latest to oldest version.
    available_dpf_servers = dict(reversed(dpf.server.available_servers().items()))
    LOGGER.info(f"Available DPF Servers: {available_dpf_servers.keys()}")

    for version, available_server in available_dpf_servers.items():
        if version in _SUPPORTED_DPF_SERVERS:
            LOGGER.info(f"Trying to launch DPF server {version}.")
            server = available_server
            break

    if server is None:
        mess = f"""Failed to launch supported DPF server:
                    Make sure one of {_SUPPORTED_DPF_SERVERS} is installed."""
        LOGGER.error(mess)
        raise SupportedDPFServerNotFoundError(mess)

    return server, version


class D3plotReader:
    """Use DPF to parse the d3plot."""

    def __init__(self, path: Path):
        """
        Initialize D3plotReader.

        Parameters
        ----------
        path : Path
            Path to the d3plot file.
        """
        _check_accept_dpf()

        # TODO: retrieve version from docker
        self._server, self._dpf_version = _get_dpf_server()

        self.ds = dpf.DataSources()
        self.ds.set_result_file_path(path, "d3plot")

        self.model = dpf.Model(self.ds)

        self.meshgrid: pv.UnstructuredGrid = self.model.metadata.meshed_region.grid
        self.time = self.model.metadata.time_freq_support.time_frequencies.data

    def get_initial_coordinates(self) -> np.ndarray:
        """Get initial coordinates."""
        return self.model.results.initial_coordinates.eval()[0].data

    def get_ep_fields(self, at_step: int = None) -> dpf.FieldsContainer:
        """Get EP fields container."""
        fields = dpf.FieldsContainer()

        time_ids = (
            self.model.metadata.time_freq_support.time_frequencies.scoping.ids
            if at_step is None
            else [at_step]
        )

        time_scoping = dpf.Scoping(ids=time_ids, location=dpf.locations.time_freq)
        # NOTE: to get time steps:
        # self.model.metadata.time_freq_support.time_frequencies.data_as_list

        op = dpf.Operator("lsdyna::ms::results")  # LS-DYNA EP operator
        op.inputs.data_sources(self.ds)
        op.inputs.time_scoping(time_scoping)
        fields = op.eval()
        return fields
        # activation_time_field = fields_container[10]

        # use to know which variable to use:
        # lsdyna::ms::result_info_provider

        # sub_fields_container: dpf.FieldsContainer = dpf.operators.utility.extract_sub_fc(
        #     fields_container=full_fields_container,
        #     label_space={"variable_id": 129},
        # ).eval()
        # sub_fields_container.animate()
        # print(self.model.operator())
        return

    # def get_transmembrane_potentials_fc(self, fc: dpf.FieldsContainer) -> dpf.FieldsContainer:
    #     """Get sub field container."""
    #     op = dpf.operators.utility.extract_sub_fc(
    #         fields_container=fc,
    #         label_space={"variable_id": 126},
    #     )
    #     return op.eval()

    #     # activation_time_field = fields_container[10]

    #     # use to know which variable to use:
    #     # lsdyna::ms::result_info_provider

    #     # sub_fields_container: dpf.FieldsContainer = dpf.operators.utility.extract_sub_fc(
    #     #     fields_container=full_fields_container,
    #     #     label_space={"variable_id": 129},
    #     # ).eval()
    #     # sub_fields_container.animate()
    #     # print(self.model.operator())
    #     return

    def print_lsdyna_ms_results(self) -> None:
        """Print available ms results."""
        # NOTE: map between variable id and variable name.
        #  Elemental Electrical Conductivity(domain Id: 1, Variable Id: 33)
        #  Elemental Scalar Potential(domain Id: 2, Variable Id: 32)
        #  Elemental Current Density(domain Id: 2, Variable Id: 1013)
        #  Elemental Electric Field(domain Id: 2, Variable Id: 1014)
        #  Elemental Ohm Heating Power(domain Id: 2, Variable Id: 35)
        #  Elemental Volumic Ohm Power(domain Id: 2, Variable Id: 100)
        #  Elemental Electrical Conductivity(domain Id: 2, Variable Id: 33)
        #  Nodal Ep Transmembrane Pot(domain Id: 3, Variable Id: 126)
        #  Nodal Ep Extra Cell Pot(domain Id: 3, Variable Id: 127)
        #  Nodal Ep Intra Cell Pot(domain Id: 3, Variable Id: 128)
        #  Nodal Ep Active. Time(domain Id: 3, Variable Id: 129)
        #  Nodal Ep Ca2+ Concentration(domain Id: 3, Variable Id: 130)
        #  Nodal (Domain Id: 3, Variable Id: 139)
        op = dpf.Operator("lsdyna::ms::result_info_provider")  # ls dyna EP operator
        op.inputs.data_sources(self.ds)
        print(op.eval())

        return

    def get_displacement_at(self, time: float) -> np.ndarray:
        """Get the displacement field.

        Parameters
        ----------
        time : float
            Time to get the displacement field at.

        Returns
        -------
        np.ndarray
            Displacement array.
        """
        if time not in self.time:
            LOGGER.warning("No data available at given time. Results are from interpolation.")
        return self.model.results.displacement.on_time_scoping(float(time)).eval()[0].data

    def get_material_ids(self) -> np.ndarray:
        """Get a list of the material IDs."""
        return self.model.metadata.meshed_region.elements.materials_field.data

    def get_history_variable(
        self,
        hv_index: List[int],
        at_step: int = 0,
    ) -> np.ndarray:
        """
        Get history variables in the d3plot.

        Parameters
        ----------
        hv_index: List[int]
            History variables index.
        at_step: int, default: 0
            Step to get the history variables at.

        Returns
        -------
        np.ndarray
            History variables data.

        Notes
        -----
        ``d3plot.get_history_variable(hv_index=list(range(9)), at_frame=at_frame)``. To
        get the deformation gradient (column-wise storage), see MAT_295 in the LS-DYNA manuals.

        """
        if at_step > self.model.metadata.time_freq_support.n_sets:
            LOGGER.warning("Frame ID doesn't exist.")
            return np.empty()

        hist_op = dpf.Operator("lsdyna::d3plot::history_var")
        time_scoping = dpf.Scoping(ids=[at_step], location=dpf.locations.time_freq)
        hist_op.connect(4, self.ds)  # why 4?
        hist_op.connect(0, time_scoping)  # why 0
        hist_vars = hist_op.eval()

        res = []
        for i in hv_index:
            if self._dpf_version.startswith("2025.2"):
                # Skip duplicate values.
                data = hist_vars[i].data[0::3]
            else:
                data = hist_vars[i].data

            res.append(data)

        return np.array(res)

    def get_heatflux(self, step: int = 2) -> np.ndarray:
        """Get nodal heat flux vector from the d3plot.

        Parameters
        ----------
        step : int, default: 2
            Time step

        Returns
        -------
        np.ndarray
            Heat flux.
        """
        op = dpf.Operator("lsdyna::d3plot::TF")
        op.inputs.data_sources(self.ds)
        time_scoping = dpf.Scoping(ids=[step], location=dpf.locations.time_freq)
        op.inputs.time_scoping(time_scoping)
        return op.eval()[0].data


class ICVoutReader:
    """Read control volume data from the binout file."""

    def __init__(self, fn: str) -> None:
        """Initialize reader.

        Parameters
        ----------
        fn : str
            Path to the binout file.
        """
        _check_accept_dpf()
        self._ds = dpf.DataSources()
        self._ds.set_result_file_path(fn, "binout")
        try:
            self._get_available_ids()
        except IndexError as error:
            LOGGER.error(f"{fn} does not contain icvout. {error}")
            raise IndexError(
                f"File {fn} does not contain control volume data. "
                "Make sure the binout file contains ICVOUT results."
            ) from error

    def _get_available_ids(self) -> np.ndarray:
        """Get available CV IDs and CVI IDs."""
        icvout_op = dpf.Operator("lsdyna::binout::ICV_ICVIID")
        icvout_op.inputs.data_sources(self._ds)
        fields1 = icvout_op.outputs.results()
        # available ICVI id
        self._icvi_ids = fields1[0].data.astype(int)

        icvout_op = dpf.Operator("lsdyna::binout::ICV_ICVID")
        icvout_op.inputs.data_sources(self._ds)
        fields2 = icvout_op.outputs.results()
        # available ICV id
        self._icv_ids = fields2[0].data.astype(int)

        return self._icv_ids

    def get_time(self) -> np.ndarray:
        """Get time array.

        Returns
        -------
        np.ndarray
            time array
        """
        # see pydpf examples, lsdyna-operators
        icvout_op = dpf.Operator("lsdyna::binout::ICV_P")
        icvout_op.inputs.data_sources(self._ds)
        p_fc = icvout_op.eval()
        rescope_op = dpf.operators.scoping.rescope()
        rescope_op.inputs.fields.connect(p_fc.time_freq_support.time_frequencies)
        rescope_op.inputs.mesh_scoping.connect(p_fc[0].scoping)
        t_field = rescope_op.outputs.fields_as_field()
        return t_field.data

    def get_pressure(self, icv_id: int) -> np.ndarray:
        """Get pressure array.

        Parameters
        ----------
        icv_id : int
            Control volume ID.

        Returns
        -------
        np.ndarray
            Pressure array.
        """
        if icv_id not in self._icv_ids:
            raise ValueError("'icv_id' is not found.")

        return self._get_field(icv_id, "ICV_P")

    def get_volume(self, icv_id: int) -> np.ndarray:
        """Get volume array.

        Parameters
        ----------
        icv_id : int
            Control volume ID.

        Returns
        -------
        np.ndarray
            Volume array.
        """
        if icv_id not in self._icv_ids:
            raise ValueError("'icv_id' not found.")

        v = self._get_field(icv_id, "ICV_V")
        # MPP bug: volume is zero at t0
        if v[0] == 0:
            v[0] = v[1]

        return v

    def get_flowrate(self, icvi_id: int) -> np.ndarray:
        """Get flow rate array.

        Parameters
        ----------
        icvi_id : int
            Control volume interaction ID.

        Returns
        -------
        np.ndarray
            Flow rate array.
        """
        if icvi_id not in self._icvi_ids:
            raise ValueError("'icvi_id' is not found.")
        # area is obtained by 'ICVI_A'
        return self._get_field(icvi_id, "ICVI_FR")

    def _get_field(self, id: int, operator_name: str) -> np.ndarray:
        icvout_op = dpf.Operator(f"lsdyna::binout::{operator_name}")
        icvout_op.inputs.data_sources(self._ds)

        my_scoping = dpf.Scoping()
        my_scoping.location = "interface"
        my_scoping.ids = [id]
        icvout_op.connect(6, my_scoping)
        fields3 = icvout_op.outputs.results()

        return fields3[0].data


class EPpostprocessor:
    """Postprocess EP (plectrophysiology) results."""

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
        """Get the field with activation times."""
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
        """Read EP results."""
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
        """Create postprocessing folder."""
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
        """Animate transmembrane potentials and export to VTK."""
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
        """Export transmembrane potentials to VTK."""
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
        """Read ECG text file produced by the LS-DYNA simulation."""
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
        """Compute 12-lead ECGs from 10 electrodes.

        Parameters
        ----------
        ECGs : np.ndarray
            mxn array containing ECGs, where m is the number of time steps
            and n is 10 electrodes in this order:
            ''V1'' ''V2'' ''V3'' ''V4'' ''V5'' ''V6'' ''RA'' ''LA'' ''RL'' ''LL''
        plot : bool, default: True
            Whether to plot.

        Returns
        -------
        np.ndarray
            12-lead ECGs in this order:
            ``I`` ``II`` ``III`` ``aVR`` ``aVL`` ``aVF`` ``V1`` ``V2`` ``V3`` ``V4`` ``V5`` ``V6``
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
        """Assign point data to the mesh."""
        result = np.zeros(self.mesh.n_points)
        result[node_ids - 1] = pointdata
        self.mesh.point_data["activation_time"] = result


class D3plotToVTKExporter:
    """Read d3plot and save the deformed mesh."""

    def __init__(self, d3plot_file: str, t_to_keep: float = 10.0e10) -> None:
        """Initialize.

        Parameters
        ----------
        d3plot_file : str
            Path to the d3plot file.
        t_to_keep : float, default: 10.0e10
            Time to convert.
        """
        self.data = D3plotReader(d3plot_file)
        self.save_time = self.data.time[self.data.time >= self.data.time[-1] - t_to_keep]

    def convert_to_pvgrid_at_t(self, time: float, fname: str = None) -> pv.UnstructuredGrid:
        """Convert d3plot data into a PyVista ``UnstructuredGrid`` object.

        Parameters
        ----------
        time : float
            Time to convert.
        fname : str, default: None
            Name of file to save data to.

        Returns
        -------
        pv.UnstructuredGrid
            Result in PyVista object.
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
        return mesh.copy()
