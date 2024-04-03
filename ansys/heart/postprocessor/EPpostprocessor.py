"""D3plot parser using Ansys-dpf."""
from math import ceil
import os
import pathlib as Path

from ansys.heart.postprocessor.dpf_utils import D3plotReader

heart_version = os.getenv("ANSYS_HEART_MODEL_VERSION")
if heart_version == "v0.2":
    from ansys.heart.preprocessor.models.v0_2.models import HeartModel
elif heart_version == "v0.1" or not heart_version:
    from ansys.heart.preprocessor.models.v0_1.models import (
        HeartModel,
    )

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
        times = self.reader.get_timesteps()
        if node_id == None:
            nnodes = len(self.reader.meshgrid.points)
            node_id = np.int64(np.linspace(0, nnodes - 1, nnodes))
        phi = np.zeros((len(times), len(node_id)))

        for time_id in range(1, len(times) + 1):
            phi[time_id - 1, :] = self.fields.get_field(
                {"variable_id": variable_id, "time": time_id}
            ).data[node_id]
        if plot == True:
            plt.plot(times, phi, label="node 0")
            plt.xlabel("time (ms)")
            plt.ylabel("phi (mV)")
            plt.show(block=True)
        return phi, times

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
            grid.point_data["transmembrane_potential"] = vm[i, :]
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
                # TODO add method to handle beam gradients as well
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

    def compute_12_lead_ECGs(
        self, ECGs: np.ndarray, times: np.ndarray, plot: bool = True
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
        """
        RA = ECGs[:, 6]
        LA = ECGs[:, 7]
        LL = ECGs[:, 9]
        I = LA - RA
        II = LL - RA
        III = LL - LA
        aVR = RA - (LA + LL) / 2
        aVL = LA - (LL + RA) / 2
        aVF = LL - (RA + LA) / 2
        Vwct = (LA + RA + LL) / 3
        V1 = ECGs[:, 0] - Vwct
        V2 = ECGs[:, 1] - Vwct
        V3 = ECGs[:, 2] - Vwct
        V4 = ECGs[:, 3] - Vwct
        V5 = ECGs[:, 4] - Vwct
        V6 = ECGs[:, 5] - Vwct
        RA = ECGs[:, 6] - Vwct
        LA = ECGs[:, 7] - Vwct
        # RL = ECGs[8, :] - Vwct
        LL = ECGs[:, 9] - Vwct
        ECGs = np.vstack((I, II, III, aVR, aVL, aVF, V1, V2, V3, V4, V5, V6))
        if plot:
            t = times
            plt.subplot(3, 4, 1)
            plt.plot(t, I)
            plt.ylabel("I")
            plt.subplot(3, 4, 5)
            plt.plot(t, II)
            plt.ylabel("II")
            plt.subplot(3, 4, 9)
            plt.plot(t, III)
            plt.ylabel("III")
            plt.subplot(3, 4, 2)
            plt.plot(t, aVR)
            plt.ylabel("aVR")
            plt.subplot(3, 4, 6)
            plt.plot(t, aVL)
            plt.ylabel("aVL")
            plt.subplot(3, 4, 10)
            plt.plot(t, aVF)
            plt.ylabel("aVF")
            plt.subplot(3, 4, 3)
            plt.plot(t, V1)
            plt.ylabel("V1")
            plt.subplot(3, 4, 7)
            plt.plot(t, V2)
            plt.ylabel("V2")
            plt.subplot(3, 4, 11)
            plt.plot(t, V3)
            plt.ylabel("V3")
            plt.subplot(3, 4, 4)
            plt.plot(t, V4)
            plt.ylabel("V4")
            plt.subplot(3, 4, 8)
            plt.plot(t, V5)
            plt.ylabel("V5")
            plt.subplot(3, 4, 12)
            plt.plot(t, V6)
            plt.ylabel("V6")
            plt.show(block=True)
            # self._plot_12LeadECGs(ECGs)
        return ECGs

    def _assign_pointdata(self, pointdata: np.ndarray, node_ids: np.ndarray):
        result = np.zeros(self.mesh.n_points)
        result[node_ids - 1] = pointdata
        self.mesh.point_data["activation_time"] = result
