"""D3plot parser using Ansys-dpf."""
import os
import pathlib as Path

from ansys.heart.postprocessor.dpf_utils import D3plotReader
from ansys.heart.preprocessor.models import HeartModel
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv


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
        if self.fields == None:
            self.fields = self.reader.get_ep_fields()

    def get_activation_times(self, at_step: int = None):
        """Get activation times field."""
        step = (
            self.reader.model.metadata.time_freq_support.time_frequencies.scoping.ids[-1]
            if at_step == None
            else [at_step]
        )
        field = self.reader.get_ep_fields(at_step=step).get_field({"variable_id": 129})
        return field

    def get_transmembrane_potential(self, node_id=None, plot: bool = False):
        """Get transmembrane potential."""
        self.load_ep_fields()
        times = self.reader.get_timesteps()
        if node_id == None:
            nnodes = len(self.reader.meshgrid.points)
            node_id = np.int64(np.linspace(0, nnodes - 1, nnodes))
        vm = np.zeros((len(times), len(node_id)))

        for time_id in range(1, len(times) + 1):
            vm[time_id - 1, :] = self.fields.get_field({"variable_id": 126, "time": time_id}).data[
                node_id
            ]
        if plot == True:
            plt.plot(times, vm, label="node 0")
            plt.xlabel("time (ms)")
            plt.ylabel("vm (mV)")
            plt.show(block=True)
        return vm, times

    # def animate_transmembrane_potentials(self):
    #     """Animate transmembrane potential."""
    #     self.load_ep_fields()
    #     tmp_fc = self.reader.get_transmembrane_potentials_fc(self.fields)
    #     tmp_fc.animate()

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

    def create_post_folder(self, path: Path = None):
        """Create Postprocessing folder."""
        if path == None:
            post_path = os.path.join(os.path.dirname(self.reader.ds.result_files[0]), "post")
        else:
            post_path = path
        isExist = os.path.exists(post_path)
        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(post_path)
        return post_path

    def animate_transmembrane(self):
        """Animate transmembrane potentials and export to vtk."""
        vm, times = self.get_transmembrane_potential()
        # Creating scene and loading the mesh
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

        for i in range(vm.shape[0]):
            # TODO vtk is not optimal for scalar fields with
            # non moving meshes, consider using ROM format
            grid.point_data["transemembrane_potential"] = vm[i, :]
            grid.save(post_path + "\\vm_" + str(i) + ".vtk")
        return

    def compute_ECGs(self, electrodes: np.ndarray):
        """Compute ECGs."""
        grid = self.reader.meshgrid
        grid = grid.compute_cell_sizes(length=False, area=False, volume=True)
        cell_volumes = grid.cell_data["Volume"]
        centroids = grid.cell_centers()
        vm, times = self.get_transmembrane_potential()
        ECGs = np.zeros([vm.shape[0], electrodes.shape[0]])

        for time_step in range(vm.shape[0]):
            grid.point_data["vmi"] = vm[time_step, :]
            grid = grid.compute_derivative(scalars="vmi")
            grid = grid.point_data_to_cell_data()

            for electrode_id in range(electrodes.shape[0]):
                electrode = electrodes[electrode_id, :]
                r_vector = centroids.points - electrode
                distances = np.linalg.norm(r_vector, axis=1)
                # TODO add conductivity tensor in the calculation (necessary?)
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
                ECGs[time_step, electrode_id] = integral
            # testing:
            # grid.point_data["testgrad"] = grid.points[:, 0]
            # grid = grid.compute_derivative(scalars="testgrad")
            # grid = grid.point_data_to_cell_data()

            # gradients = get vm gradient on elem centroids
            # volumes = get element volumes
            # q1 = r / (np.power(distances, 3) * 4 * np.pi)
            # ECGi = q1 .gradients(t) x volumes
        return ECGs, times

    def read_ECGs(self, path: Path):
        """Read ECG text file produced by LS-DYNA simulation."""
        data = np.loadtxt(path, skiprows=4)
        times = data[:, 0]
        ECGs = data[:, 1:11]
        return ECGs, times

    def _assign_pointdata(self, pointdata: np.ndarray, node_ids: np.ndarray):
        result = np.zeros(self.mesh.n_points)
        result[node_ids - 1] = pointdata
        self.mesh.point_data["activation_time"] = result
