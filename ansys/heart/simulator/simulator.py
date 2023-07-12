"""Simulator module.

Options for simulation:
- EP-only
    with/without fbers
    with/without purkinje
- Electro-mechanics
    simplified EP (imposed activation)
    coupled electro-mechanics
"""
import copy
import glob as glob
import json
import os
import pathlib as Path
import shutil
import subprocess
from typing import Literal

from ansys.heart.misc.element_orth import read_orth_element_kfile
from ansys.heart.postprocessor.D3plotPost import LVContourExporter
from ansys.heart.postprocessor.Klotz_curve import EDPVR
from ansys.heart.postprocessor.SystemModelPost import SystemModelPost
from ansys.heart.postprocessor.dpf_d3plot import D3plotReader
from ansys.heart.preprocessor.mesh.objects import Cavity, SurfaceMesh
from ansys.heart.preprocessor.models import HeartModel, LeftVentricle
from ansys.heart.simulator.settings.settings import SimulationSettings
import ansys.heart.writer.dynawriter as writers
import numpy as np


def which(program):
    """Return path if program exists, else None."""

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


class BaseSimulator:
    """Base class for the simulator."""

    def __init__(
        self,
        model: HeartModel,
        lsdynapath: Path,
        dynatype: Literal["smp", "intelmpi", "platformmpi", "msmpi"],
        num_cpus: int = 1,
        simulation_directory: Path = "",
    ) -> None:
        """Initialize BaseSimulator.

        Parameters
        ----------
        model : HeartModel
            Heart model to simulate.
        lsdynapath : Path
            Path to LS-DYNA executable.
        dynatype : Literal[&quot;smp&quot;, &quot;intelmpi&quot;, &quot;platformmpi&quot;]
            Type of LS-DYNA executable. shared memory parallel or massively parallel process.
        num_cpus : int, optional
            Number of cpu's to use for simulation, by default 1
        simulation_directory : Path, optional
            Directory of simulation, by default defined by information in HeartModel.
        """
        self.model: HeartModel = model
        """HeartModel to simulate."""
        self.lsdynapath = lsdynapath
        """Path to LS-DYNA executable."""
        self.dynatype = dynatype
        """LS-DYNA Type"""
        self.num_cpus = num_cpus
        """Number of cpus to use for simulation."""
        self.directories: dict = {}
        """Dictionary of all defined directories."""

        if which(lsdynapath) is None:
            print(f"{lsdynapath} not exist")
            exit()

        if simulation_directory == "":
            simulation_directory = os.path.join(self.model.info.workdir, "simulation")

        self.root_directory = simulation_directory
        """Root simulation directory."""

        self.settings: SimulationSettings = SimulationSettings()
        """Simulation settings."""

        pass

    def load_default_settings(self) -> SimulationSettings:
        """Load default simulation settings."""
        self.settings.load_defaults()
        return self.settings

    def compute_fibers(self):
        """Compute the fiber direction on the model."""
        directory = self._write_fibers()

        print("Computing fiber orientation...")

        # self.settings.save(os.path.join(directory, "simulation_settings.yml"))
        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)

        print("done.")

        # interpolate new fibers onto model.mesh
        # Number of cells/points or element/node ordering may not be the same
        # especially in the case of a full-heart model where we do not use
        # the full heart to compute the fibers. Hence, interpolate using the cell
        # centers.
        # NOTE: How to handle null values?

        # # read results.
        # print("Interpolating fibers onto model.mesh")
        # vtk_with_fibers = os.path.join(directory, "vtk_FO_ADvectors.vtk")
        # vtk_with_fibers = pyvista.UnstructuredGrid(vtk_with_fibers)
        #
        # cell_centers_target = vtk_with_fibers.cell_centers()
        # cell_centers_source = self.model.mesh.cell_centers()
        #
        # cell_centers_source = cell_centers_source.interpolate(cell_centers_target)
        #
        # self.model.mesh.cell_data["fiber"] = cell_centers_source.point_data["aVector"]
        # self.model.mesh.cell_data["sheet"] = cell_centers_source.point_data["dVector"]
        # print("Done.")

        print("Assigning fiber orientation to model...")
        elem_ids, part_ids, connect, fib, sheet = read_orth_element_kfile(
            os.path.join(directory, "element_solid_ortho.k")
        )

        self.model.mesh.cell_data["fiber"][elem_ids - 1] = fib
        self.model.mesh.cell_data["sheet"][elem_ids - 1] = sheet

        # dump the model to reuse fiber information
        self.model.dump_model(os.path.join(self.root_directory, "model_with_fiber.pickle"))
        return

    def _run_dyna(self, path_to_input: Path, options: str = ""):
        """Run LS-DYNA with path and options.

        Parameters
        ----------
        path_to_input : Path
            Path to the LS-DYNA simulation file.
        options : str, optional
            Additional options to pass to command line, by default ""
        """
        os.chdir(os.path.dirname(path_to_input))

        if self.dynatype in ["intelmpi", "platformmpi", "msmpi"]:
            commands = [
                "mpirun",
                "-np",
                str(self.num_cpus),
                self.lsdynapath,
                "i=" + path_to_input,
                options,
            ]
        elif self.dynatype in ["smp"]:
            commands = [
                self.lsdynapath,
                "i=" + path_to_input,
                "ncpu=" + str(self.num_cpus),
                options,
            ]

        # launch LS-DYNA
        with subprocess.Popen(commands, stdout=subprocess.PIPE, text=True) as p:
            for line in p.stdout:
                print(line.rstrip())

        os.chdir(self.root_directory)
        return

    def _write_fibers(
        self,
        alpha_endocardium: float = -60,
        alpha_epicardium: float = 60,
        beta_endocardium: float = 25,
        beta_epicardium: float = -65,
    ) -> Path:
        """Write LS-DYNA files for fiber generation."""
        export_directory = os.path.join(self.root_directory, "fibergeneration")
        self.directories["fibergeneration"] = export_directory

        dyna_writer = writers.FiberGenerationDynaWriter(copy.deepcopy(self.model), self.settings)
        dyna_writer.update()
        dyna_writer.export(export_directory)

        return export_directory


class EPSimulator(BaseSimulator):
    """EP Simulator."""

    def __init__(
        self,
        model: HeartModel,
        lsdynapath: Path,
        dynatype: Literal["smp", "intelmpi", "platformmpi"],
        num_cpus: int = 1,
        simulation_directory: Path = "",
    ) -> None:
        super().__init__(model, lsdynapath, dynatype, num_cpus, simulation_directory)

        return

    def simulate(self, folder_name="main-ep"):
        """Launch the main simulation."""
        directory = os.path.join(self.root_directory, folder_name)
        self._write_main_simulation_files(folder_name)

        print("Launching main EP simulation...")

        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)

        print("done.")

        return

    def compute_purkinje(self):
        """Compute the purkinje network."""
        directory = self._write_purkinje_files()

        print("Computing the Purkinje network...")

        # self.settings.save(os.path.join(directory, "simulation_settings.yml"))
        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)

        print("done.")

        print("Assign the Purkinje network to the model...")
        purkinje_files = glob.glob(os.path.join(directory, "purkinjeNetwork_*.k"))
        for purkinje_file in purkinje_files:
            self.model.mesh.add_purkinje_from_kfile(purkinje_file)

    def compute_conduction_system(self):
        """Compute the conduction system."""
        self.model.mesh.compute_av_conduction()

    def _write_main_simulation_files(self, folder_name):
        """Write LS-DYNA files that are used to start the main simulation."""
        export_directory = os.path.join(self.root_directory, folder_name)
        self.directories["main-ep"] = export_directory

        dyna_writer = writers.ElectrophysiologyDynaWriter(self.model, self.settings)
        dyna_writer.update()
        dyna_writer.export(export_directory)

        return export_directory

    def _write_purkinje_files(
        self,
        pointstx: float = 0,  # TODO instantiate this
        pointsty: float = 0,  # TODO instantiate this
        pointstz: float = 0,  # TODO instantiate this
        inodeid: int = 0,  # TODO instantiate this
        iedgeid: int = 0,  # TODO instantiate this
        edgelen: float = 2,  # TODO instantiate this
        ngen: float = 50,
        nbrinit: int = 8,
        nsplit: int = 2,
    ) -> Path:
        """Write purkinje files.

        Parameters
        ----------
        pointstx : float, optional
            _description_, by default 0
        pointsty : float, optional
            _description_, by default 0
        pointstz : float, optional
            _description_, by default 0
        nbrinit : int, optional
            _description_, by default 8
        nsplit : int, optional
            _description_, by default 2
        """
        export_directory = os.path.join(self.root_directory, "purkinjegeneration")
        self.directories["purkinjegeneration"] = export_directory

        dyna_writer = writers.PurkinjeGenerationDynaWriter(copy.deepcopy(self.model), self.settings)
        dyna_writer.update()
        dyna_writer.export(export_directory)
        return export_directory


class MechanicsSimulator(BaseSimulator):
    """Mechanics simulator with imposed active stress."""

    def __init__(
        self,
        model: HeartModel,
        lsdynapath: Path,
        dynatype: Literal["smp", "intelmpi", "platformmpi"],
        num_cpus: int = 1,
        simulation_directory: Path = "",
        initial_stress: bool = True,
    ) -> None:
        super().__init__(model, lsdynapath, dynatype, num_cpus, simulation_directory)

        """If stress free computation is taken into considered."""
        # include initial stress by default
        self.initial_stress = initial_stress

        """A dictionary save stress free computation information"""
        self.stress_free_report = None

        return

    def simulate(self, folder_name="main-mechanics", zerop_folder=None, auto_post=True):
        """
        Launch the main simulation.

        Parameters
        ----------
        zerop_folder : str
            folder contains stress free simulation.
            Default is "zeropressure" under roo_directory.
        auto_post : bool
            if run post-process scripts.
        folder_name: str
            main simulation folder name.

        """
        directory = os.path.join(self.root_directory, folder_name)
        os.makedirs(directory, exist_ok=True)

        if zerop_folder is None:
            zerop_folder = os.path.join(self.root_directory, "zeropressure")

        if self.initial_stress:
            try:
                # get dynain.lsda file from
                dynain_file = glob.glob(os.path.join(zerop_folder, "iter*.dynain.lsda"))[-1]

                shutil.copy(dynain_file, os.path.join(directory, "dynain.lsda"))
                shutil.copy(
                    os.path.join(zerop_folder, "post", "Post_report.json"),
                    os.path.join(directory, "Post_report.json"),
                )
            except IndexError:
                # handle if lsda file not exist.
                print(
                    "Cannot find initial stress file, simulation will run without initial stress."
                )
                self.initial_stress = False

        self._write_main_simulation_files(folder_name=folder_name)

        print("Launching main simulation...")

        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)

        print("done.")

        def _post():
            os.makedirs(os.path.join(directory, "post"), exist_ok=True)

            system = SystemModelPost(directory)
            fig = system.plot_pv_loop()
            fig.savefig(os.path.join(directory, "post", "pv.png"))

            fig = system.plot_pressure_flow_volume(system.lv_system)
            fig.savefig(os.path.join(directory, "post", "lv.png"))

            if not isinstance(self.model, LeftVentricle):
                fig = system.plot_pressure_flow_volume(system.rv_system)
                fig.savefig(os.path.join(directory, "post", "rv.png"))

            exporter = LVContourExporter(os.path.join(directory, "d3plot"), self.model)

            self.model.compute_left_ventricle_anatomy_axis(first_cut_short_axis=0.2)
            exporter.export_contour_to_vtk("l4cv", self.model.l4cv_axis)
            exporter.export_contour_to_vtk("l2cv", self.model.l2cv_axis)
            normal = self.model.short_axis["normal"]
            p_start = self.model.short_axis["center"]
            for ap in self.model.left_ventricle.apex_points:  # use next()?
                if ap.name == "apex epicardium":
                    p_end = ap.xyz

            for icut in range(2):
                p_cut = p_start + (p_end - p_start) * icut / 2
                cutter = {"center": p_cut, "normal": normal}
                exporter.export_contour_to_vtk(f"shor_{icut}", cutter)

            exporter.export_lvls_to_vtk("lvls")

        if auto_post:
            _post()
        return

    def compute_stress_free_configuration(self):
        """Compute the stress-free configuration of the model."""
        directory = self._write_stress_free_configuration_files()

        print("Computing stress-free configuration...")

        self.settings.save(Path.Path(directory) / "simulation_settings.yml")
        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file, options="case")

        print("done.")

        def _post():
            """Post process zeropressure folder."""
            os.makedirs(os.path.join(directory, "post"), exist_ok=True)
            data = D3plotReader(glob.glob(os.path.join(directory, "iter*.d3plot"))[-1])
            expected_time = self.settings.stress_free.analysis.end_time.to("millisecond").m

            if data.time[-1] != expected_time:
                Warning("Stress free computation is not converged, skip post process.")
                return None

            stress_free_coord = data.get_initial_coordinates()
            displacements = data.get_displacement()

            if len(self.model.cap_centroids) == 0:
                nodes = self.model.mesh.nodes
            else:
                # a center node for each cap has been created, add them into create the cavity
                nodes = np.vstack((self.model.mesh.nodes, np.zeros((len(self.model.cap_centroids), 3))))
                for cap_center in self.model.cap_centroids:
                    nodes[cap_center.node_id] = cap_center.xyz

            # convergence information
            dst = np.linalg.norm(
                stress_free_coord + displacements[-1] - nodes, axis=1
            )
            error_mean = np.mean(dst)
            error_max = np.max(dst)

            # geometry information
            self.model.mesh.save(os.path.join(directory, "post", "True_ED.vtk"))
            tmp = copy.deepcopy(self.model)
            tmp.mesh.points = stress_free_coord
            tmp.mesh.save(os.path.join(directory, "post", "zerop.vtk"))
            tmp.mesh.points = stress_free_coord + displacements[-1]
            tmp.mesh.save(os.path.join(directory, "post", "Simu_ED.vtk"))

            def _post_cavity(name: str):
                """Extract cavity volume."""
                try:
                    faces = (
                        np.loadtxt(
                            os.path.join(directory, name + ".segment"), delimiter=",", dtype=int
                        )
                        - 1
                    )
                except FileExistsError:
                    print(f"Cannot find {name}.segment")

                volumes = []
                for i, dsp in enumerate(displacements):
                    cavity_surface = SurfaceMesh(
                        name=name, triangles=faces, nodes=stress_free_coord + dsp
                    )
                    cavity_surface.save(os.path.join(directory, "post", f"{name}_{i}.stl"))
                    volumes.append(Cavity(cavity_surface).volume)

                return volumes

            # cavity information
            lv_volumes = _post_cavity("left_ventricle")
            # todo: not safe to use list index
            true_lv_ed_volume = self.model.cavities[0].volume
            volume_error = [(lv_volumes[-1] - true_lv_ed_volume) / true_lv_ed_volume]

            # Klotz curve information
            # unit is mL and mmHg
            lv_pr_mmhg = (
                self.settings.mechanics.boundary_conditions.end_diastolic_cavity_pressure[
                    "left_ventricle"
                ]
                .to("mmHg")
                .m
            )

            klotz = EDPVR(true_lv_ed_volume / 1000, lv_pr_mmhg)
            sim_vol_ml = [v / 1000 for v in lv_volumes]
            sim_pr = lv_pr_mmhg * data.time / data.time[-1]

            fig = klotz.plot_EDPVR(simulation_data=[sim_vol_ml, sim_pr])
            fig.savefig(os.path.join(directory, "post", "klotz.png"))

            dct = {
                "Simulation output time (ms)": data.time.tolist(),
                "Left ventricle EOD pressure (mmHg)": lv_pr_mmhg,
                "True left ventricle volume (mm3)": true_lv_ed_volume,
                "Simulation Left ventricle volume (mm3)": lv_volumes,
                "Convergence": {
                    "max_error (mm)": error_max,
                    "mean_error (mm)": error_mean,
                    "relative volume error (100%)": volume_error,
                },
            }

            # right ventricle exist
            if not isinstance(self.model, LeftVentricle):
                rv_volumes = _post_cavity("right_ventricle")
                true_rv_ed_volume = self.model.cavities[1].volume
                volume_error.append((rv_volumes[-1] - true_rv_ed_volume) / true_rv_ed_volume)

                rv_pr_mmhg = (
                    self.settings.mechanics.boundary_conditions.end_diastolic_cavity_pressure[
                        "right_ventricle"
                    ]
                    .to("mmHg")
                    .m
                )
                dct["Right ventricle EOD pressure (mmHg)"] = rv_pr_mmhg
                dct["True right ventricle volume"] = true_rv_ed_volume
                dct["Simulation Right ventricle volume"] = rv_volumes

            return dct

        self.stress_free_report = _post()
        with open(os.path.join(directory, "post", "Post_report.json"), "w") as f:
            json.dump(self.stress_free_report, f)

        return

    def _write_main_simulation_files(self, folder_name):
        """Write LS-DYNA files that are used to start the main simulation."""
        export_directory = os.path.join(self.root_directory, folder_name)
        self.directories["main-mechanics"] = export_directory

        dyna_writer = writers.MechanicsDynaWriter(
            self.model,
            self.settings,
        )
        dyna_writer.update(with_dynain=self.initial_stress)
        dyna_writer.export(export_directory)

        return export_directory

    def _write_stress_free_configuration_files(self) -> Path:
        """Write LS-DYNA files to compute stress-free configuration."""
        export_directory = os.path.join(self.root_directory, "zeropressure")
        self.directories["zeropressure"] = export_directory

        dyna_writer = writers.ZeroPressureMechanicsDynaWriter(self.model, self.settings)
        dyna_writer.update()
        dyna_writer.export(export_directory)

        return export_directory


class EPMechanicsSimulator(EPSimulator, MechanicsSimulator):
    """Coupled EP-mechanics simulator with computed Electrophysiology."""

    def __init__(
        self,
        model: HeartModel,
        lsdynapath: Path,
        dynatype: Literal["smp", "intelmpi", "platformmpi"],
    ) -> None:
        super().__init__(model, lsdynapath, dynatype)
        raise NotImplementedError("Simulator EPMechanicsSimulator not implemented.")
