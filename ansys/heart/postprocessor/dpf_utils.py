"""D3plot parser using Ansys-dpf."""
import os
import pathlib as Path
from typing import List

from ansys.dpf import core as dpf
import numpy as np

# from ansys.dpf.core.dpf_operator import available_operator_names


class D3plotReader:
    """Use DPF to parse d3plot."""

    def __init__(self, path: Path.Path):
        """
        Initialize D3plotReader.

        Parameters
        ----------
        path : Path
            d3plot file path.
        """
        os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"
        # os.environ["ANSYS_DPF_SERVER_CONTEXT"] = "PREMIUM"
        # dpf.set_default_server_context(dpf.AvailableServerContexts.premium)

        # server = dpf.start_local_server()

        self.ds = dpf.DataSources()
        self.ds.set_result_file_path(path, "d3plot")

        self.model = dpf.Model(self.ds)

        # common
        self.meshgrid = self.model.metadata.meshed_region.grid
        self.time = self.model.metadata.time_freq_support.time_frequencies.data

    def get_initial_coordinates(self):
        """Get initial coordinates."""
        return self.model.results.initial_coordinates.eval()[0].data

    def get_timesteps(self):
        """Get list of timesteps."""
        return self.model.metadata.time_freq_support.time_frequencies.data_as_list

    def get_ep_fields(self, at_step: int = None) -> dpf.FieldsContainer:
        """Get EP fields container."""
        fields = dpf.FieldsContainer()

        time_ids = (
            self.model.metadata.time_freq_support.time_frequencies.scoping.ids
            if at_step == None
            else [at_step]
        )

        time_scoping = dpf.Scoping(ids=time_ids, location=dpf.locations.time_freq)
        # to get time steps:
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

    def get_transmembrane_potentials_fc(self, fc: dpf.FieldsContainer) -> dpf.FieldsContainer:
        """Get sub field container."""
        op = dpf.operators.utility.extract_sub_fc(
            fields_container=fc,
            label_space={"variable_id": 126},
        )
        return op.eval()

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

    def print_lsdyna_ms_results(self):
        """Print available ms results."""
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

    def get_displacement(self):
        """Get displacement."""
        displacements = self.model.results.displacement.on_all_time_freqs
        fields = displacements.eval()
        res = []
        for f in fields:
            res.append(f.data)
        return res

    def get_history_variable(
        self,
        hv_index: List[int],
        at_frame: int = 0,
    ):
        """
                Get history variables in d3plot.

        Parameters
        ----------
        hv_index: List[int]
            History variables index.
        at_frame: int, optional
            At this frame, by default 0.

        Returns
        -------
        np.ndarray
            History variables data.

        Notes
        -----
        d3plot.get_history_variable(hv_index=list(range(9)), at_frame=at_frame)
        to get Deformation gradient (column-wise storage),see *MAT_295 in LS-DYNA manual.
        """
        if at_frame > self.model.metadata.time_freq_support.n_sets:
            print("Frame ID too big")
            exit()

        hist_op = dpf.Operator("lsdyna::d3plot::history_var")
        time_scoping = dpf.Scoping(ids=[at_frame], location=dpf.locations.time_freq)
        hist_op.connect(4, self.ds)  # why 4?
        hist_op.connect(0, time_scoping)  # why 0
        hist_vars = hist_op.eval()

        res = []
        for i in hv_index:
            res.append(hist_vars[i].data)

        return np.array(res)

    def export_vtk(
        self,
        file_path: str,
        prefix: str = "model",
        only_surface: bool = False,
        keep_mat_ids: List[int] = None,
    ):
        """
        Convert d3plot to vtk.

        Parameters
        ----------
        keep_mat_ids: keep cells by material IDs.
        only_surface: extract surface.
        file_path: Path
        prefix: vtk file prefix, optional

        """
        mat_ids = self.model.metadata.meshed_region.elements.materials_field.data

        if not np.all(np.isin(keep_mat_ids, np.unique(mat_ids))):
            print("Invalid material IDs, all parts will be saved.")
            keep_mat_ids = None

        if keep_mat_ids is None:
            # keep all cells
            keep_cells = np.ones((len(mat_ids)), dtype=bool)
        else:
            # keep cells by material ID
            keep_cells = np.zeros((len(mat_ids)), dtype=bool)
            for id in keep_mat_ids:
                keep_cells = keep_cells | (mat_ids == id)

        # displacements = self.model.results.displacement.on_all_time_freqs
        # fields = displacements.eval()
        fields = self.get_displacement()

        for i, field in enumerate(fields):
            # get pyvista grid
            vista_grid = self.meshgrid.copy()

            # deform mesh by displacement
            vista_grid.points += field.data

            # keep cells
            vista_grid = vista_grid.extract_cells(keep_cells)

            # extract surface
            if only_surface:
                vista_grid = vista_grid.extract_surface()

            # export
            vista_grid.save(os.path.join(file_path, f"{prefix}_{i}.vtk"))
            # pyvista.save_meshio(os.path.join(file_path, f"{prefix}_{i}.vtk"),vista_grid)

        # This needs Premium dpf licence
        # op = dpf.operators.serialization.vtk_export(
        #     export_type=0,
        #     file_path=os.path.join(file_path,f"{prefix}_{i}.vtk"),
        #     mesh=mesh,
        #     # fields1=field,
        #     # fields2=my_fields2,
        # ).eval()

        return


if __name__ == "__main__":
    pass
