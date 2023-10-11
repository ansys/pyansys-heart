"""D3plot parser using Ansys-dpf."""
import os
import pathlib as Path

from ansys.heart.postprocessor.dpf_utils import D3plotReader
from ansys.heart.preprocessor.models import HeartModel
import numpy as np


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

    def load_all_ep_fields(self):
        """Load all EP fields."""
        self.fields = self.reader.get_ep_fields()

    def get_activation_times(self, at_step: int = None):
        """Get activation times field."""
        step = (
            self.reader.model.metadata.time_freq_support.time_frequencies.scoping.ids[-1]
            if at_step == None
            else [at_step]
        )
        field = self.reader.get_ep_fields(at_step=step)[10]
        return field

    # def get_transmembrane_potentials(self, at_step: int = 0):
    #     # field = self.fields.get_field({"domain_id": 3, "time": at_step, "variable_id": 126})
    #     return field

    # def animate_transmembrane(self, start_time, endtime):
    #     return

    def read_EP_nodout(self):
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

    def _assign_pointdata(self, pointdata: np.ndarray, node_ids: np.ndarray):
        result = np.zeros(self.mesh.n_points)
        result[node_ids - 1] = pointdata
        self.mesh.point_data["activation_time"] = result
