"""D3plot parser using Ansys-dpf."""
import os
import pathlib as Path
from typing import List

from ansys.dpf import core as dpf
import numpy as np

# from ansys.dpf.core.dpf_operator import available_operator_names


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
        os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"
        # os.environ["ANSYS_DPF_SERVER_CONTEXT"] = "PREMIUM"
        # dpf.set_default_server_context(dpf.AvailableServerContexts.premium)

        # server = dpf.start_local_server()

        ds = dpf.DataSources()
        ds.set_result_file_path(path, "d3plot")

        self.model = dpf.Model(ds)
        self.results = self.model.results
        self.time = self.model.metadata.time_freq_support.time_frequencies.data

    # todo extract history output: F, active stress...

    def get_displacement(self):
        """Get displacement."""
        displacements = self.results.displacement.on_all_time_freqs
        fields = displacements.eval()
        res = []
        for f in fields:
            res.append(f.data)
        return res

    def get_initial_coordinates(self):
        """Get initial coordinates."""
        return self.results.initial_coordinates.eval()[0].data

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
        displacements = self.results.displacement.on_all_time_freqs
        fields = displacements.eval()

        mesh = self.model.metadata.meshed_region
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

        for i, field in enumerate(fields):
            # convert to pyvista grid
            vista_grid = mesh.grid.copy()

            # deform mesh by displacement
            vista_grid.points += field.data

            # keep cells
            vista_grid = vista_grid.extract_cells(keep_cells)

            # extract surface
            if only_surface:
                vista_grid = vista_grid.extract_surface()

            # export
            vista_grid.save(os.path.join(file_path, f"{prefix}_{i}.vtk"))

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
    data = D3plotReader("D:\Heart20\healthy20\health03_BV_2mm\simulation\zeropressure\iter3.d3plot")
    # export to vtk:
    data.export_vtk(".", keep_mat_ids=[1, 3])

# dpf.operators.serialization.field_to_csv(
#     file_path="d:\\development\\temp6\\displacements.csv",
#     field_or_fields_container=fields,
#     # my_storage_type=0
# ).eval()
