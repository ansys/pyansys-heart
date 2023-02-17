"""D3plot parser using Ansys-dpf."""
import os
import pathlib as Path

from ansys.dpf import core as dpf

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
        server = dpf.start_local_server()

        ds = dpf.DataSources()
        ds.set_result_file_path(path, "d3plot")

        self.model = dpf.Model(ds)
        self.results = self.model.results

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


# # export to vtk:
# op = dpf.operators.serialization.vtk_export(
#     export_type=0,
#     file_path="test.vtk",
#     mesh=model.metadata.meshed_region,
#     fields1=fields,
#     # fields2=my_fields2,
# ).eval()

# dpf.operators.serialization.field_to_csv(
#     file_path="d:\\development\\temp6\\displacements.csv",
#     field_or_fields_container=fields,
#     # my_storage_type=0
# ).eval()
