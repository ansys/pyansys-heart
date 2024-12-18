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

"""D3plot parser using Ansys-dpf."""

import os
from pathlib import Path
from typing import List

import numpy as np
import psutil
import pyvista as pv

from ansys.dpf import core as dpf
from ansys.heart.core import LOG as LOGGER

_KILL_ANSYSCL_ON_DEL: bool = False
"""Flag indicating whether to kill the ansys licence client upon deleting the D3PlotReader."""


def _check_env():
    if "ANSYS_DPF_ACCEPT_LA" in os.environ and os.environ["ANSYS_DPF_ACCEPT_LA"] == "Y":
        pass
    else:
        LOGGER.warning(
            """DPF requires you to accept the license agreement.
            Set the environment variable "ANSYS_DPF_ACCEPT_LA" to "Y"."""
        )
        exit()
    return


class SupportedDPFServerNotFoundError(Exception):
    """SupportedDPFServerNotFoundError."""

    pass


class D3plotReader:
    """Use DPF to parse d3plot."""

    def __init__(self, path: Path):
        """
        Initialize D3plotReader.

        Parameters
        ----------
        path : Path
            d3plot file path.
        """
        _check_env()

        _running_ansyscl_pids = set([p.pid for p in psutil.process_iter() if "ansyscl" in p.name()])

        try:
            LOGGER.info("Trying to launch DPF Server 2024.1")
            self._server = dpf.server.available_servers()["2024.1"]()
        except KeyError as e:
            LOGGER.error(f"Failed to launch DPF Server 2024.1. {e}")
            raise SupportedDPFServerNotFoundError(f"Failed to launch DPF Server 2024.1. {e}")

        self.ds = dpf.DataSources()
        self.ds.set_result_file_path(path, "d3plot")

        self.model = dpf.Model(self.ds)

        self._ansyscl_pid = [
            p.pid
            for p in psutil.process_iter()
            if "ansyscl" in p.name() and p.pid not in _running_ansyscl_pids
        ]
        """ansyscl process id triggered by (Py)DPF."""

        self.meshgrid: pv.UnstructuredGrid = self.model.metadata.meshed_region.grid
        self.time = self.model.metadata.time_freq_support.time_frequencies.data

    def __del__(self):
        """Force shutdown ansyscl after use of dpf."""
        # NOTE: kills the ansys client triggered by (py)dpf.
        # May resolve issues when using tools using other versions of the license client.
        if _KILL_ANSYSCL_ON_DEL:
            try:
                for pid in self._ansyscl_pid:
                    psutil.Process(pid).kill()
            except Exception as e:
                LOGGER.info(f"Failed to kill ansyscl upon delete. {e}")

    def get_initial_coordinates(self):
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

        op = dpf.Operator("lsdyna::ms::results")  # ls dyna EP operator
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

    def print_lsdyna_ms_results(self):
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

    def get_displacement_at(self, time: float) -> np.ndarray:
        """Get displacement field.

        Parameters
        ----------
        time : float
            at which time

        Returns
        -------
        np.ndarray
            displacement
        """
        if time not in self.time:
            LOGGER.warning("No data at given time, results are from interpolation.")
        return self.model.results.displacement.on_time_scoping(float(time)).eval()[0].data

    def get_material_ids(self):
        """Get list of material id."""
        return self.model.metadata.meshed_region.elements.materials_field.data

    def get_history_variable(
        self,
        hv_index: List[int],
        at_step: int = 0,
    ):
        """
        Get history variables in d3plot.

        Parameters
        ----------
        hv_index: List[int]
            History variables index.
        at_step: int, optional
            At this frame, by default 0.

        Returns
        -------
        np.ndarray
            History variables data.

        Notes
        -----
        d3plot.get_history_variable(hv_index=list(range(9)), at_frame=at_frame) to
        get Deformation gradient (column-wise storage),see MAT_295 in LS-DYNA manual.

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
            res.append(hist_vars[i].data)

        return np.array(res)


class ICVoutReader:
    """Read control volume data from binout."""

    def __init__(self, fn: str) -> None:
        """Init reader.

        Parameters
        ----------
        fn : str
            binout file path
        """
        _check_env()
        self._ds = dpf.DataSources()
        self._ds.set_result_file_path(fn, "binout")
        try:
            self._get_available_ids()
        except IndexError:
            LOGGER.error(f"{fn} do not contain icvout.")
            exit()

    def _get_available_ids(self):
        """Get available CV ids and CVI ids."""
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
            control volume id

        Returns
        -------
        np.ndarray
            pressure array
        """
        if icv_id not in self._icv_ids:
            raise ValueError("icv_id not found.")

        return self._get_field(icv_id, "ICV_P")

    def get_volume(self, icv_id: int) -> np.ndarray:
        """Get volume array.

        Parameters
        ----------
        icv_id : int
            control volume id

        Returns
        -------
        np.ndarray
            volume array
        """
        if icv_id not in self._icv_ids:
            raise ValueError("icv_id not found.")

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
            control volume interaction id

        Returns
        -------
        np.ndarray
            flowrate array
        """
        if icvi_id not in self._icvi_ids:
            raise ValueError("icvi_id not found.")
        # area is obtained by 'ICVI_A'
        return self._get_field(icvi_id, "ICVI_FR")

    def _get_field(self, id: int, operator_name: str):
        icvout_op = dpf.Operator(f"lsdyna::binout::{operator_name}")
        icvout_op.inputs.data_sources(self._ds)

        my_scoping = dpf.Scoping()
        my_scoping.location = "interface"
        my_scoping.ids = [id]
        icvout_op.connect(6, my_scoping)
        fields3 = icvout_op.outputs.results()

        return fields3[0].data


if __name__ == "__main__":
    pass
