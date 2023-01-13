"""Simulator module.

Options for simulation:
- EP-only
    with/without fbers
    with/without purkinje
- Electro-mechanics
    simplified EP (imposed activation)
    coupled electro-mechanics
"""
import glob as glob
import os
import pathlib as Path
import shutil
import subprocess
from typing import Literal

from ansys.heart.preprocessor.models import HeartModel
import ansys.heart.writer.dynawriter as writers
import numpy as np
import pyvista


class BaseSimulator:
    """Base class for the simulator."""

    def __init__(
        self,
        model: HeartModel,
        lsdynapath: Path,
        dynatype: Literal["smp", "intelmpi", "platformmpi"],
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

        if simulation_directory == "":
            simulation_directory = os.path.join(self.model.info.workdir, "simulation")

        self.root_directory = simulation_directory
        """Root simulation directory."""
        pass

    def compute_fibers(self):
        """Compute the fiber direction on the model."""
        print("Computing fiber orientation...")

        directory = self._write_fibers()
        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)

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
        vtk_with_fibers = pyvista.UnstructuredGrid(vtk_with_fibers)

        cell_centers_target = vtk_with_fibers.cell_centers()
        cell_centers_source = self.model.mesh.cell_centers()

        cell_centers_source = cell_centers_source.interpolate(cell_centers_target)

        self.model.mesh.cell_data["fiber"] = cell_centers_source.point_data["aVector"]
        self.model.mesh.cell_data["sheet"] = cell_centers_source.point_data["dVector"]
        print("Done.")

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

        if self.dynatype in ["intelmpi"]:
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
        p = subprocess.run(commands, stdout=subprocess.PIPE)

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

        dyna_writer = writers.FiberGenerationDynaWriter(self.model)
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

    def simulate(self):
        """Launch the main simulation."""
        print("Launching main simulation.")
        directory = self._write_main_simulation_files()
        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)
        print("done.")
        return

    def compute_purkinje(self):
        """Compute the purkinje network."""
        print("Computing the Purkinje network...")

        directory = self._write_purkinje_files()
        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)

        purkinje_files = glob.glob(os.path.join(directory, "purkinjeNetwork_*.k"))
        for purkinje_file in purkinje_files:
            self.model.mesh.add_purkinje_from_kfile(purkinje_file)

        print("done.")

    def _write_main_simulation_files(self):
        """Write LS-DYNA files that are used to start the main simulation."""
        export_directory = os.path.join(self.root_directory, "main-ep")
        self.directories["main-ep"] = export_directory

        dyna_writer = writers.ElectrophysiologyDynaWriter(self.model)
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

        dyna_writer = writers.PurkinjeGenerationDynaWriter(self.model)
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
    ) -> None:
        super().__init__(model, lsdynapath, dynatype, num_cpus, simulation_directory)

        return

    def simulate(self):
        """Launch the main simulation."""
        print("Launching main simulation.")
        directory = self._write_main_simulation_files()
        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)
        print("done.")
        return

    def compute_stress_free_configuration(self):
        """Compute the stress-free configuration of the model."""
        print("Computing stress-free configuration...")

        directory = self._write_stress_free_configuration_files()
        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file, options="case")

        print("done.")

        Warning("Replace by methods from postprocessing module.")
        binout_files = glob.glob(os.path.join(directory, "iter*.binout"))
        iter_files = glob.glob(os.path.join(directory, "iter*.guess"))
        # read nodes of last iteration. (avoid qd dependency for now.)
        stress_free_nodes = np.loadtxt(
            iter_files[-1], skiprows=2, max_rows=self.model.mesh.nodes.shape[0]
        )
        self.model.mesh.points = stress_free_nodes[:, 1:]
        return

    def _write_main_simulation_files(self):
        """Write LS-DYNA files that are used to start the main simulation."""
        export_directory = os.path.join(self.root_directory, "main-mechanics")
        self.directories["main-mechanics"] = export_directory

        dyna_writer = writers.MechanicsDynaWriter(self.model, "ConstantPreloadWindkesselAfterload")
        dyna_writer.update()
        dyna_writer.export(export_directory)

        return export_directory

    def _write_stress_free_configuration_files(self) -> Path:
        """Write LS-DYNA files to compute stress-free configuration."""
        export_directory = os.path.join(self.root_directory, "zeropressure")
        self.directories["zeropressure"] = export_directory

        dyna_writer = writers.ZeroPressureMechanicsDynaWriter(self.model)
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


class Simulator:
    """
    Perform pre-simulation steps.

    Some extra info here

    """

    def __init__(
        self,
        model: HeartModel,
        lsdynapath: Path,
        lsdyna_type: Literal["smp", "intelmpi"],
        num_cpus: int = 1,
    ) -> None:
        DeprecationWarning("This class is deprecated")
        self.model = model
        """Heart model."""
        self.lsdynapath = lsdynapath
        """Path of the lsdyna executable."""

    def _write_fibers(
        self,
        workdir: str,
        alpha_endocardium: float = -60,
        alpha_eepicardium: float = 60,
        beta_endocardium: float = 25,
        beta_epicardium: float = -65,
    ):
        dyna_writer = writers.FiberGenerationDynaWriter(self.model)
        dyna_writer.update()
        dyna_writer.export(workdir)

        return

    def _write_zeropressureconfiguration(self, workdir: str = ""):
        dyna_writer = writers.ZeroPressureMechanicsDynaWriter(self.model)
        dyna_writer.update()
        dyna_writer.export(workdir)
        return

    def _write_purkinje(
        self,
        workdir: str,
        pointstx: float = 0,  # TODO instantiate this
        pointsty: float = 0,  # TODO instantiate this
        pointstz: float = 0,  # TODO instantiate this
        inodeid: int = 0,  # TODO instantiate this
        iedgeid: int = 0,  # TODO instantiate this
        edgelen: float = 2,  # TODO instantiate this
        ngen: float = 50,
        nbrinit: int = 8,
        nsplit: int = 2,
    ):
        """Write purkinje files.

        Parameters
        ----------
        workdir : str
            path where files are dumped
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
        dyna_writer = writers.PurkinjeGenerationDynaWriter(self.model)
        dyna_writer.update()
        dyna_writer.export(workdir)
        return

    def __get_stressfreenodes(workdir: str):
        """Get stres free nodes."""
        # TODO check if converged
        guess_files = []
        for file in os.listdir(workdir):
            if file[-5:] == "guess":
                guess_files.append(file)

        return guess_files[-1]

    def build(
        self,
        path_lsdyna: str,
        path_simulation: str,
        fibers: bool = False,
        purkinje: bool = False,
        zeropressure: bool = False,
        ep: bool = False,
        mechanics: bool = False,
    ):
        """Build LS-DYNA Heart simulation."""
        # TODO add getters and setters for fiber angles, purkinje properties, simulation times and
        # other parameters to expose to the user
        path_simulation = os.path.join(self.model.info.path_to_model, "simulation")
        simulationdynawriter = writers.BaseDynaWriter(self.model)

        if fibers:
            path_fibers = os.path.join(self.model.info.path_to_model, "fiber_generation")
            self._write_fibers(path_fibers)
            # TODO run dyna

        if zeropressure:
            path_zeropressure = os.path.join(
                self.model.info.path_to_model, "zeropressure_generation"
            )
            self._write_zeropressureconfiguration(path_zeropressure)
            if fibers:
                shutil.copy2(
                    os.path.join(path_fibers, "element_solid_ortho.k"),
                    os.path.join(path_zeropressure, "solid_elements.k"),
                )

            # TODO run dyna
            nodes_stressfree = self.__get_stressfreenodes(path_zeropressure)
            shutil.copy(
                os.path.join(path_zeropressure, "nodes.k"),
                os.path.join(path_zeropressure, "nodes_endofdiastole.k"),
            )
            shutil.copy(
                os.path.join(path_zeropressure, nodes_stressfree),
                os.path.join(path_zeropressure, "nodes.k"),
            )

        if purkinje:
            path_purkinje = os.path.join(self.model.info.path_to_model, "purkinje_generation")
            # if exist left ventricle:
            simulationdynawriter.model.add_part("Left Purkinje")
            simulationdynawriter.include_files.append("purkinjenetworkLEFT.k")
            # if exist right ventricle:
            simulationdynawriter.model.add_part("Right Purkinje")
            simulationdynawriter.include_files.append("purkinjenetworkRIGHT.k")
            self._write_purkinje(path_purkinje)
            if zeropressure:
                shutil.copy2(
                    os.path.join(path_zeropressure, "nodes.k"),
                    os.path.join(path_purkinje, "nodes.k"),
                )
            # if exist right ventricle:
            # TODO run dyna mainLeft
            shutil.copy2(
                os.path.join(path_purkinje, "purkinjenetwork.k"),
                os.path.join(path_purkinje, "purkinjenetworkLEFT.k"),
            )
            shutil.copy2(
                os.path.join(path_purkinje, "purkinjenetworkLEFT.k"),
                os.path.join(path_simulation, "purkinjenetworkLEFT.k"),
            )
            # TODO if exist right ventricle:
            # TODO run dyna mainRight
            shutil.copy2(
                os.path.join(path_purkinje, "purkinjenetwork.k"),
                os.path.join(path_purkinje, "purkinjenetworkRIGHT.k"),
            )
            shutil.copy2(
                os.path.join(path_purkinje, "purkinjenetworkRIGHT.k"),
                os.path.join(path_simulation, "purkinjenetworkRIGHT.k"),
            )

        # if (ep) and not (mechanics):
        # if not (ep) and (mechanics):
        # write mechanics
        # if ep and mechanics:
        # TODO add coupling stuff and ignore default
        # Ca2+ mechanics active stress/replaced by EP simulation

        return
