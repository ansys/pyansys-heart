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
import os
import pathlib as Path
import shutil
import subprocess
from typing import Literal

from ansys.heart.postprocessor.auto_process import mech_post, zerop_post
from ansys.heart.preprocessor.models_new import FourChamber, FullHeart, HeartModel
from ansys.heart.simulator.settings.settings import SimulationSettings
import ansys.heart.writer.dynawriter as writers
import pyvista as pv


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
        dynatype: Literal["smp", "intelmpi", "platformmpi", "msmpi"] = "smp",
        num_cpus: int = 1,
        simulation_directory: Path = "",
        platform: Literal["windows", "wsl", "linux"] = "windows",
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
        platform : str, by default "windows"
            Operating system.
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
        self.platform = platform
        """Operating System."""
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
        self._run_dyna(path_to_input=input_file)

        print("done.")

        # interpolate new fibers onto model.mesh
        # Number of cells/points or element/node ordering may not be the same
        # especially in the case of a full-heart model where we do not use
        # the full heart to compute the fibers. Hence, interpolate using the cell
        # centers.
        # NOTE: How to handle null values?

        # read results.
        print("Interpolating fibers onto model.mesh")
        vtk_with_fibers = os.path.join(directory, "vtk_FO_ADvectors.vtk")
        vtk_with_fibers = pv.UnstructuredGrid(vtk_with_fibers)

        cell_centers_target = vtk_with_fibers.cell_centers()
        cell_centers_source = self.model.mesh.cell_centers()

        cell_centers_source = cell_centers_source.interpolate(cell_centers_target)

        self.model.mesh.cell_data["fiber"] = cell_centers_source.point_data["aVector"]
        self.model.mesh.cell_data["sheet"] = cell_centers_source.point_data["dVector"]
        print("Done.")

        # from ansys.heart.misc.element_orth import read_orth_element_kfile
        # print("Assigning fiber orientation to model...")
        # elem_ids, part_ids, connect, fib, sheet = read_orth_element_kfile(
        #     os.path.join(directory, "element_solid_ortho.k")
        # )

        # self.model.mesh.cell_data.set_vectors(fib, name="fiber", deep_copy=True)
        # self.model.mesh.cell_data.set_vectors(sheet, name="sheet", deep_copy=True)

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
        run_lsdyna(
            lsdynapath=self.lsdynapath,
            path_to_input=path_to_input,
            simulation_directory=self.root_directory,
            options=options,
            dynatype=self.dynatype,
            num_cpus=self.num_cpus,
            platform=self.platform,
        )

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
        platform: Literal["windows", "wsl", "linux"] = "windows",
    ) -> None:
        super().__init__(
            model, lsdynapath, dynatype, num_cpus, simulation_directory, platform=platform
        )

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
        if isinstance(self.model, (FourChamber, FullHeart)):
            self.model.compute_av_conduction()
            self.model.compute_His_conduction()
            self.model.compute_bundle_branches()

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
        platform: Literal["windows", "wsl", "linux"] = "windows",
    ) -> None:
        super().__init__(
            model, lsdynapath, dynatype, num_cpus, simulation_directory, platform=platform
        )

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

        if auto_post:
            mech_post(directory, self.model)
        return

    def compute_stress_free_configuration(self, folder_name="zeropressure"):
        """Compute the stress-free configuration of the model."""
        directory = os.path.join(self.root_directory, folder_name)
        os.makedirs(directory, exist_ok=True)

        self._write_stress_free_configuration_files(folder_name)
        self.settings.save(Path.Path(directory) / "simulation_settings.yml")

        print("Computing stress-free configuration...")
        self._run_dyna(os.path.join(directory, "main.k"), options="case")
        print("Simulation done.")

        self.stress_free_report = zerop_post(directory, self)

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

    def _write_stress_free_configuration_files(self, folder_name) -> Path:
        """Write LS-DYNA files to compute stress-free configuration."""
        export_directory = os.path.join(self.root_directory, folder_name)
        self.directories["zeropressure"] = export_directory

        dyna_writer = writers.ZeroPressureMechanicsDynaWriter(self.model, self.settings)
        dyna_writer.update()
        dyna_writer.export(export_directory)

        return


class EPMechanicsSimulator(EPSimulator, MechanicsSimulator):
    """Coupled EP-mechanics simulator with computed Electrophysiology."""

    def __init__(
        self,
        model: HeartModel,
        lsdynapath: Path,
        dynatype: Literal["smp", "intelmpi", "platformmpi"],
        platform: Literal["windows", "wsl", "linux"] = "windows",
    ) -> None:
        super().__init__(model, lsdynapath, dynatype, platform=platform)
        raise NotImplementedError("Simulator EPMechanicsSimulator not implemented.")


def run_lsdyna(
    lsdynapath: Path,
    path_to_input: Path,
    simulation_directory: Path = "",
    options: str = "",
    dynatype: str = "intelmpi",
    num_cpus: int = 1,
    platform="windows",
):
    """Standalone function for running LS-DYNA.

    Parameters
    ----------
    lsdynapath : Path
        Lsdyna path.
    path_to_input : Path
        Path to input.
    simulation_directory : Path, optional
        Simulation directory, by default"".
    options : str, optional
        Options, by default "".
    dynatype : str, optional
        Dynatype, by default "intelmpi".
    num_cpus : int, optional
        number of cpus, by default 1.
    platform : str
        Operating system, by default "windows".
    """
    os.chdir(os.path.dirname(path_to_input))
    if platform == "windows" or platform == "linux":
        if dynatype in ["intelmpi", "platformmpi", "msmpi"]:
            commands = [
                "mpirun",
                "-np",
                str(num_cpus),
                lsdynapath,
                "i=" + path_to_input,
                options,
            ]
        elif dynatype in ["smp"]:
            commands = [
                lsdynapath,
                "i=" + path_to_input,
                "ncpu=" + str(num_cpus),
                options,
            ]
    elif platform == "wsl":
        path_to_input = (
            subprocess.run(["wsl", "wslpath", os.path.basename(path_to_input)], capture_output=1)
            .stdout.decode()
            .strip()
        )
        lsdynapath = (
            subprocess.run(["wsl", "wslpath", str(lsdynapath).replace("\\", "/")], capture_output=1)
            .stdout.decode()
            .strip()
        )

        if dynatype in ["intelmpi", "platformmpi", "msmpi"]:
            commands = [
                "mpirun",
                "-np",
                str(num_cpus),
                lsdynapath,
                "i=" + path_to_input,
                options,
            ]
        elif dynatype in ["smp"]:
            commands = [
                lsdynapath,
                "i=" + path_to_input,
                "ncpu=" + str(num_cpus),
                options,
            ]

        with open("run_lsdyna.sh", "w", newline="\n") as f:
            f.write("#!/usr/bin/env sh\n")
            f.write("echo start lsdyna in wsl...\n")
            f.write(" ".join([i.strip() for i in commands]))

        commands = ["powershell", "-Command", "wsl", "-e", "bash", "-lic", "./run_lsdyna.sh"]

    with subprocess.Popen(commands, stdout=subprocess.PIPE, text=True) as p:
        for line in p.stdout:
            print(line.rstrip())

    os.chdir(simulation_directory)
    return
