"""D3plot parser using Ansys-dpf."""
import os
import pathlib as Path

from ansys.heart.preprocessor.models import HeartModel
import numpy as np


class EPpostprocessor:
    """Postprocess Electrophysiology results."""

    def __init__(self, results_path: Path, model: HeartModel):
        """Initialize EP postprocessor.

        Parameters
        ----------
        results_path : Path
            Path to results files.
        model : HeartModel
            Heartmodel to be post-processed.
        """
        self.mesh = model.mesh
        self.results_path = results_path
        self.times = None
        self.node_ids: np.ndarray = None

    def read_EP_results(self, type: str = "nodout"):
        """Read Electrophysiology results.

        Parameters
        ----------
        type : str, optional
            type of EP result, by default "nodout"
        """
        if type == "nodout":
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
