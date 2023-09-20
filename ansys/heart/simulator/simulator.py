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
from typing import List, Literal

from ansys.heart.custom_logging import LOGGER
from ansys.heart.misc.element_orth import read_orth_element_kfile
from ansys.heart.postprocessor.auto_process import mech_post, read_uvc, zerop_post
from ansys.heart.preprocessor.models import FourChamber, FullHeart, HeartModel
from ansys.heart.simulator.settings.settings import SimulationSettings
import ansys.heart.writer.dynawriter as writers
import yaml


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


class DynaSettings:
    """Class for collecting, managing and validating LS-DYNA settings."""

    def __init__(
        self,
        lsdyna_path: Path,
        dynatype: Literal["smp", "intelmpi", "platformmpi", "msmpi"] = "intelmpi",
        num_cpus: int = 1,
        platform: Literal["windows", "wsl", "linux"] = "windows",
        dyna_options: str = "",
        mpi_options: str = "",
    ):
        """Initialize Dyna settings.

        Parameters
        ----------
        lsdyna_path : Path
            Path to LS-DYNA
        dynatype : Literal[&quot;smp&quot;, &quot;intelmpi&quot;, &quot;platformmpi&quot;]
            Type of LS-DYNA executable. Shared Memory Parallel or Massively Parallel Processing
        num_cpus : int, optional
            Number of CPU's requested, by default 1
        platform : Literal[&quot;windows&quot;, &quot;wsl&quot;, &quot;linux&quot;], optional
            Platform, by default "windows"
        dyna_options : str, optional
            Additional command line options, by default ""
        mpi_options : str, optional
            Additional mpi options, by default ""
        """
        self.lsdyna_path: Path = lsdyna_path
        """Path to LS-DYNA executable."""
        self.dynatype: str = dynatype
        """Type of dyna executable."""
        self.num_cpus: int = num_cpus
        """Number of CPU's requested."""
        self.platform: str = platform
        """Platform on which dyna is executed."""

        self.dyna_options = dyna_options
        """Additional command line options for dyna."""

        if dynatype in ["intelmpi", "platformmpi", "msmpi"]:
            self.mpi_options = mpi_options
            """additional mpi options."""
        elif dynatype == "smp":
            self.mpi_options = ""

        return

    def get_commands(self, path_to_input: Path) -> List[str]:
        # """Get command line arguments from the defined settings."""
        """Get command line arguments from the defined settings.

        Parameters
        ----------
        path_to_input : Path
            Path to the LS-DYNA input file.

        Returns
        -------
        List[str]
            List of strings of each of the commands.
        """
        lsdyna_path = self.lsdyna_path

        if self.platform == "windows" or self.platform == "linux":
            if self.dynatype in ["intelmpi", "platformmpi"]:
                commands = [
                    "mpirun",
                    self.mpi_options,
                    "-np",
                    str(self.num_cpus),
                    lsdyna_path,
                    "i=" + path_to_input,
                    self.dyna_options,
                ]
            elif self.dynatype in ["smp"]:
                commands = [
                    lsdyna_path,
                    "i=" + path_to_input,
                    "ncpu=" + str(self.num_cpus),
                    self.dyna_options,
                ]
        if self.platform == "windows" and self.dynatype == "msmpi":
            commands = [
                "mpiexec",
                self.mpi_options,
                "-np",
                str(self.num_cpus),
                lsdyna_path,
                "i=" + path_to_input,
                self.dyna_options,
            ]

        elif self.platform == "wsl":
            path_to_input_wsl = (
                subprocess.run(
                    ["wsl", "wslpath", os.path.basename(path_to_input)], capture_output=1
                )
                .stdout.decode()
                .strip()
            )
            # redefines LS-DYNA path.
            lsdyna_path = (
                subprocess.run(
                    ["wsl", "wslpath", str(lsdyna_path).replace("\\", "/")], capture_output=1
                )
                .stdout.decode()
                .strip()
            )

            if self.dynatype in ["intelmpi", "platformmpi", "msmpi"]:
                commands = [
                    "mpirun",
                    self.mpi_options,
                    "-np",
                    str(self.num_cpus),
                    lsdyna_path,
                    "i=" + path_to_input_wsl,
                    self.dyna_options,
                ]
            elif self.dynatype in ["smp"]:
                commands = [
                    lsdyna_path,
                    "i=" + path_to_input_wsl,
                    "ncpu=" + str(self.num_cpus),
                    self.dyna_options,
                ]

            with open("run_lsdyna.sh", "w", newline="\n") as f:
                f.write("#!/usr/bin/env sh\n")
                f.write("echo start lsdyna in wsl...\n")
                f.write(" ".join([i.strip() for i in commands]))

            commands = ["powershell", "-Command", "wsl", "-e", "bash", "-lic", "./run_lsdyna.sh"]

        # remove empty strings from commands
        commands = [c for c in commands if c != ""]

        return commands

    def _set_env_variables(self):
        r"""Try to set environment variables for MPI run using Ansys installation root directories.

        Notes
        -----
        This is based on lsdynaintelvar.bat and lsdynamsvar.bat in
        ANSYS Inc\v231\ansys\bin\winx64\lsprepost48\LS-Run 1.0
        and requires you to install the MPI libraries with the Ansys installer.
        """
        # get latest installed Ansys version
        env_variables = list(os.environ.keys())
        ansys_env_roots = sorted([k for k in env_variables if "AWP_ROOT" in k])
        ansys_env_roots.reverse()

        if self.platform == "windows":
            platform = "winx64"
        elif self.platform == "linux":
            platform = "linx64"

        if self.dynatype in "intelmpi":
            INTEL_CMP_REV = "2019.5.281"
            INTEL_MKL_REV = "2020.0.166"
            INTEL_MPI_REV = "2018.3.210"

            if os.getenv("MPI_ROOT"):
                LOGGER.warning(
                    "MPI_ROOT already defined. Not trying to automatically set MPI env variables."
                )
                return

            for root in ansys_env_roots:
                mpi_root = os.path.join(
                    os.getenv(root), "commonfiles", "MPI", "Intel", INTEL_MPI_REV, platform
                )
                mpi_path = os.path.join(mpi_root, "bin")
                mkl_path = os.path.join(os.getenv(root), "tp", "IntelMKL", INTEL_MKL_REV, platform)
                cmp_path = os.path.join(
                    os.getenv(root), "tp", "IntelCompiler", INTEL_CMP_REV, platform
                )
                for p in [mpi_root, mpi_path, mkl_path, cmp_path]:
                    if not os.path.isdir(p):
                        LOGGER.debug("Failed to set env variables with %s " % root)
                        continue

                LOGGER.info("Setting MPI environment variable MPI_ROOT to %s" % mpi_root)
                LOGGER.info("Adding %s to PATH" % [mpi_path, mkl_path, cmp_path])

                os.environ["MPI_ROOT"] = mpi_root
                os.environ["PATH"] = (
                    ";".join([mpi_path, mkl_path, cmp_path]) + ";" + os.environ["PATH"]
                )
                os.environ["I_MPI_AUTH_METHOD"] = "delegate"
                os.environ["KMP_AFFINITY"] = "verbose"
                break

        elif self.dynatype == "msmpi" and self.platform == "windows":
            INTEL_CMP_REV = "2019.5.281"
            INTEL_MKL_REV = "2020.0.166"
            MS_MPI_REV = "10.1.12498.18"

            if os.getenv("MPI_ROOT"):
                LOGGER.warning(
                    "MPI_ROOT already defined. Not trying to automatically set MPI env variables."
                )
                return

            for root in ansys_env_roots:
                mpi_root = os.path.join(
                    os.getenv(root), "commonfiles", "MPI", "Microsoft", MS_MPI_REV, platform
                )
                mpi_path = os.path.join(mpi_root, "bin")
                mkl_path = os.path.join(os.getenv(root), "tp", "IntelMKL", INTEL_MKL_REV, platform)
                cmp_path = os.path.join(
                    os.getenv(root), "tp", "IntelCompiler", INTEL_CMP_REV, platform
                )
                for p in [mpi_root, mpi_path, mkl_path, cmp_path]:
                    if not os.path.isdir(p):
                        LOGGER.debug("Failed to set env variables with %s " % root)
                        continue

                LOGGER.info("Setting MPI environment variable MPI_ROOT to %s" % mpi_root)
                LOGGER.info("Adding %s to PATH" % [mpi_path, mkl_path, cmp_path])

                os.environ["MPI_ROOT"] = mpi_root
                os.environ["PATH"] = (
                    ";".join([mpi_path, mkl_path, cmp_path]) + ";" + os.environ["PATH"]
                )
                break

        elif self.dynatype == "platformmpi":
            LOGGER.error("Automatically setting env variables for platform mpi not implemented yet")
            return

        return

    def __repr__(self):
        """Represent self as string."""
        return yaml.dump(vars(self), allow_unicode=True, default_flow_style=False)


class BaseSimulator:
    """Base class for the simulator."""

    def __init__(
        self,
        model: HeartModel,
        dynasettings: DynaSettings = None,
        simulation_directory: Path = "",
    ) -> None:
        """Initialize BaseSimulator.

        Parameters
        ----------
        model : HeartModel
            Heart model to simulate.
        dynasettings : DynaSettings
            Settings used for launching LS-DYNA.
        simulation_directory : Path, optional
            Directory in which to start the simulation, by default ""
        """
        self.model: HeartModel = model
        """HeartModel to simulate."""
        if not dynasettings:
            LOGGER.warning("Setting default LS-DYNA settings.")
            self.dyna_settings = DynaSettings()
        else:
            self.dyna_settings: DynaSettings = dynasettings
            """Contains the settings to launch LS-DYNA."""

        self.directories: dict = {}
        """Dictionary of all defined directories."""

        """Operating System."""
        if which(self.dyna_settings.lsdyna_path) is None:
            LOGGER.error(f"{self.dyna_settings.lsdyna_path} not exist")
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

        LOGGER.info("Computing fiber orientation...")

        # self.settings.save(os.path.join(directory, "simulation_settings.yml"))
        input_file = os.path.join(directory, "main.k")
        self._run_dyna(path_to_input=input_file)

        LOGGER.info("done.")

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

        LOGGER.info("Assigning fiber orientation to model...")
        elem_ids, part_ids, connect, fib, sheet = read_orth_element_kfile(
            os.path.join(directory, "element_solid_ortho.k")
        )

        self.model.mesh.cell_data["fiber"][elem_ids - 1] = fib
        self.model.mesh.cell_data["sheet"][elem_ids - 1] = sheet

        # dump the model to reuse fiber information
        self.model.dump_model(os.path.join(self.root_directory, "model_with_fiber.pickle"))
        return

    def compute_uvc(
        self,
    ):
        """Compute universal 'heart' coordinates system."""
        if isinstance(self.model, FullHeart):
            raise NotImplementedError("Not yet tested for the full heart")

        LOGGER.info("Computing universal ventricular coordinates...")

        dirname = "uvc"
        export_directory = os.path.join(self.root_directory, dirname)
        self.directories[dirname] = export_directory

        # Dyna writer
        dyna_writer = writers.UHCWriter(
            copy.deepcopy(self.model),
        )
        dyna_writer.update()
        dyna_writer.export(export_directory)

        input_file = os.path.join(export_directory, "main.k")
        self._run_dyna(path_to_input=input_file, options="case")

        LOGGER.info("done.")

        grid = read_uvc(export_directory)

    def _run_dyna(self, path_to_input: Path, options: str = ""):
        """Run LS-DYNA with path and options.

        Parameters
        ----------
        path_to_input : Path
            Path to the LS-DYNA simulation file.
        options : str, optional
            Additional options to pass to command line, by default ""
        """
        if options != "":
            old_options = copy.deepcopy(self.dyna_settings.dyna_options)
            self.dyna_settings.dyna_options = self.dyna_settings.dyna_options + " " + options

        run_lsdyna(
            path_to_input=path_to_input,
            settings=self.dyna_settings,
            simulation_directory=self.root_directory,
        )

        if options != "":
            self.dyna_settings.dyna_options = old_options

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
        dyna_settings: DynaSettings,
        simulation_directory: Path = "",
    ) -> None:
        super().__init__(model, dyna_settings, simulation_directory)

        return

    def simulate(self, folder_name="main-ep"):
        """Launch the main simulation."""
        directory = os.path.join(self.root_directory, folder_name)
        self._write_main_simulation_files(folder_name)

        LOGGER.info("Launching main EP simulation...")

        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)

        LOGGER.info("done.")

        return

    def compute_purkinje(self):
        """Compute the purkinje network."""
        directory = self._write_purkinje_files()

        LOGGER.info("Computing the Purkinje network...")

        # self.settings.save(os.path.join(directory, "simulation_settings.yml"))
        input_file = os.path.join(directory, "main.k")

        LOGGER.debug("Compute Purkinje network on 1 cpu.")
        orig_num_cpus = self.dyna_settings.num_cpus
        self.dyna_settings.num_cpus = 1
        self._run_dyna(input_file)
        self.dyna_settings.num_cpus = orig_num_cpus

        LOGGER.info("done.")

        LOGGER.info("Assign the Purkinje network to the model...")
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
        dyna_settings: DynaSettings,
        simulation_directory: Path = "",
        initial_stress: bool = True,
    ) -> None:
        super().__init__(model, dyna_settings, simulation_directory)

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
                LOGGER.warning(
                    "Cannot find initial stress file, simulation will run without initial stress."
                )
                self.initial_stress = False

        self._write_main_simulation_files(folder_name=folder_name)

        LOGGER.info("Launching main simulation...")

        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)

        LOGGER.info("done.")

        if auto_post:
            mech_post(directory, self.model)
        return

    def compute_stress_free_configuration(self, folder_name="zeropressure"):
        """Compute the stress-free configuration of the model."""
        directory = os.path.join(self.root_directory, folder_name)
        os.makedirs(directory, exist_ok=True)

        self._write_stress_free_configuration_files(folder_name)
        self.settings.save(Path.Path(directory) / "simulation_settings.yml")

        LOGGER.info("Computing stress-free configuration...")
        self._run_dyna(os.path.join(directory, "main.k"), options="case")
        LOGGER.info("Simulation done.")

        self.stress_free_report = zerop_post(directory, self.model)

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

    def __init__(self, model: HeartModel, dyna_settings: DynaSettings) -> None:
        super().__init__(model, dyna_settings)
        raise NotImplementedError("Simulator EPMechanicsSimulator not implemented.")


def run_lsdyna(
    path_to_input: Path,
    settings: DynaSettings = None,
    simulation_directory: Path = None,
):
    """Standalone function for running LS-DYNA.

    Parameters
    ----------
    path_to_input : Path
        Input file for LS-DYNA.
    settings : DynaSettings, optional
        LS-DYNA settings, such as path to executable, executable type, platform, by default None
    simulation_directory : Path, optional
        Directory where to simulate, by default None
    """
    if not settings:
        LOGGER.info("Using default LS-DYNA settings.")
        settings = DynaSettings()

    commands = settings.get_commands(path_to_input)

    os.chdir(os.path.dirname(path_to_input))

    with subprocess.Popen(commands, stdout=subprocess.PIPE, text=True) as p:
        for line in p.stdout:
            LOGGER.info(line.rstrip())

    os.chdir(simulation_directory)
    return
