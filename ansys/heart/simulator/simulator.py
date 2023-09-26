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

from ansys.heart.custom_logging import LOGGER
from ansys.heart.misc.element_orth import read_orth_element_kfile
from ansys.heart.postprocessor.auto_process import mech_post, read_uvc, zerop_post
from ansys.heart.preprocessor.models_new import FourChamber, FullHeart, HeartModel
from ansys.heart.simulator.settings.settings import DynaSettings, SimulationSettings
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

        LOGGER.info("Assigning fiber orientation to model...")
        elem_ids, part_ids, connect, fib, sheet = read_orth_element_kfile(
            os.path.join(directory, "element_solid_ortho.k")
        )

        # self.model.mesh.cell_data.set_vectors(fib, name="fiber", deep_copy=True)
        # self.model.mesh.cell_data.set_vectors(sheet, name="sheet", deep_copy=True)

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
