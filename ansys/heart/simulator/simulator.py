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
from ansys.heart.postprocessor.dpf_d3plot import D3plotReader
from ansys.heart.preprocessor.mesh.objects import Cavity, SurfaceMesh
from ansys.heart.preprocessor.models import HeartModel, LeftVentricle
from ansys.heart.simulator.settings.settings import SimulationSettings
import ansys.heart.writer.dynawriter as writers
import numpy as np


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

        self.settings.save(os.path.join(directory, "simulation_settings.yml"))
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

    def simulate(self):
        """Launch the main simulation."""
        directory = self._write_main_simulation_files()

        print("Launching main EP simulation...")

        self.settings.save(os.path.join(directory, "simulation_settings.yml"))
        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)

        print("done.")
        return

    def compute_purkinje(self):
        """Compute the purkinje network."""
        directory = self._write_purkinje_files()

        print("Computing the Purkinje network...")

        self.settings.save(os.path.join(directory, "simulation_settings.yml"))
        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)

        print("done.")

        print("Assign the Purkinje network to the model...")
        purkinje_files = glob.glob(os.path.join(directory, "purkinjeNetwork_*.k"))
        for purkinje_file in purkinje_files:
            self.model.mesh.add_purkinje_from_kfile(purkinje_file)

    def _write_main_simulation_files(self):
        """Write LS-DYNA files that are used to start the main simulation."""
        export_directory = os.path.join(self.root_directory, "main-ep")
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

        return

    def simulate(self):
        """Launch the main simulation."""
        directory = self._write_main_simulation_files()
        input_file = os.path.join(directory, "main.k")

        if self.initial_stress:
            try:
                # get dynain.lsda file from
                dynain_file = glob.glob(
                    os.path.join(self.root_directory, "zeropressure", "iter*.dynain.lsda")
                )[-1]

                shutil.copy(dynain_file, os.path.join(directory, "dynain.lsda"))
            except IndexError:
                # handle if lsda file not exist.
                print(
                    "Cannot find initial stress file, simulation will run without initial stress."
                )
                self.initial_stress = False
                self._write_main_simulation_files()

        print("Launching main simulation...")

        self.settings.save(os.path.join(directory, "simulation_settings.yml"))
        self._run_dyna(input_file)

        print("done.")
        return

    def compute_stress_free_configuration(self):
        """Compute the stress-free configuration of the model."""
        directory = self._write_stress_free_configuration_files()

        print("Computing stress-free configuration...")

        self.settings.save(os.path.join(directory, "simulation_settings.yml"))
        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file, options="case")

        print("done.")

        def _post():
            """Post process zeropressure folder."""
            data = D3plotReader(glob.glob(os.path.join(directory, "iter*.d3plot"))[-1])
            stress_free_coord = data.get_initial_coordinates()
            displacements = data.get_displacement()

            # convergence information
            dst = np.linalg.norm(
                stress_free_coord + displacements[-1] - self.model.mesh.points, axis=1
            )
            error_mean = np.mean(dst)
            error_max = np.max(dst)

            # geometry information
            self.model.mesh.save(os.path.join(directory, "True_ED.vtk"))
            tmp = copy.deepcopy(self.model)
            tmp.mesh.points = stress_free_coord
            tmp.mesh.save(os.path.join(directory, "zerop.vtk"))
            tmp.mesh.points = stress_free_coord + displacements[-1]
            tmp.mesh.save(os.path.join(directory, "Simu_ED.vtk"))

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
                    cavity_surface.save(os.path.join(directory, f"{name}_{i}.stl"))
                    volumes.append(Cavity(cavity_surface).volume)

                return volumes

            # cavity information
            lv_volumes = _post_cavity("left_ventricle")
            true_lv_ed_volume = self.model.cavities[0].volume
            volume_error = [(lv_volumes[-1] - true_lv_ed_volume) / true_lv_ed_volume]

            dct = {
                "True left ventricle volume": true_lv_ed_volume,
                "Simulation Left ventricle volume": lv_volumes,
                "Convergence": {
                    "max_error": error_max,
                    "mean_error": error_mean,
                    "relative volume error": volume_error,
                },
            }

            if not isinstance(self.model, LeftVentricle):
                rv_volumes = _post_cavity("right_ventricle")
                true_rv_ed_volume = self.model.cavities[1].volume
                volume_error.append((rv_volumes[-1] - true_rv_ed_volume) / true_rv_ed_volume)

                dct["True right ventricle volume"] = true_rv_ed_volume
                dct["Simulation Right ventricle volume"] = rv_volumes

            # dump info
            with open(os.path.join(directory, "Post_report.json"), "w") as f:
                json.dump(dct, f)

        _post()

        return

    def _write_main_simulation_files(self):
        """Write LS-DYNA files that are used to start the main simulation."""
        export_directory = os.path.join(self.root_directory, "main-mechanics")
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
