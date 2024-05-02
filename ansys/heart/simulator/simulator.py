# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Simulator module.

Options for simulation:

- EP-only
    with/without fbers.
    with/without purkinje.
- Electro-mechanics
    simplified EP (imposed activation).
    coupled electro-mechanics.
"""

import copy
import os
import pathlib as Path
import shutil
import subprocess
from typing import List, Literal

from ansys.heart.core import LOG as LOGGER
from ansys.heart.misc.element_orth import read_orth_element_kfile
from ansys.heart.postprocessor.auto_process import (
    compute_la_fiber_cs,
    compute_ra_fiber_cs,
    mech_post,
    read_uvc,
    zerop_post,
)
from ansys.heart.preprocessor.mesh.objects import Part
from ansys.heart.simulator.settings.material.material import NeoHookean

global heart_version
heart_version = os.getenv("ANSYS_HEART_MODEL_VERSION")
if heart_version == "v0.2":
    from ansys.heart.preprocessor.models.v0_2.models import FourChamber, HeartModel, LeftVentricle
elif heart_version == "v0.1" or not heart_version:
    from ansys.heart.preprocessor.models.v0_1.models import FourChamber, HeartModel, LeftVentricle

    heart_version = "v0.1"
from ansys.heart.preprocessor.models.conduction_beam import ConductionSystem
from ansys.heart.simulator.settings.settings import DynaSettings, SimulationSettings
import ansys.heart.writer.dynawriter as writers
import numpy as np
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
        dyna_settings: DynaSettings = None,
        simulation_directory: Path = "",
    ) -> None:
        """Initialize BaseSimulator.

        Parameters
        ----------
        model : HeartModel
            Heart model to simulate.
        dyna_settings : DynaSettings
            Settings used for launching LS-DYNA.
        simulation_directory : Path, optional
            Directory in which to start the simulation, by default ""

        """
        self.model: HeartModel = model
        """HeartModel to simulate."""
        if not dyna_settings:
            LOGGER.warning("Setting default LS-DYNA settings.")
            self.dyna_settings = DynaSettings()
        else:
            self.dyna_settings: DynaSettings = dyna_settings
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

    def compute_uhc(self) -> pv.UnstructuredGrid:
        """Compute universal 'heart' coordinates system."""
        LOGGER.info("Computing universal ventricular coordinates...")

        type = "uvc"
        export_directory = os.path.join(self.root_directory, type)

        target = self.run_laplace_problem(export_directory, type)
        grid = read_uvc(export_directory)

        grid["cell_ids"] = target["cell_ids"]
        grid["point_ids"] = target["point_ids"]

        LOGGER.info("Assigning data to full model...")

        ids = grid["point_ids"]
        self.model.mesh["apico-basal"] = np.zeros(self.model.mesh.n_points)
        self.model.mesh["apico-basal"][:] = np.nan
        self.model.mesh["transmural"] = np.zeros(self.model.mesh.n_points)
        self.model.mesh["transmural"][:] = np.nan
        self.model.mesh["rotational"] = np.zeros(self.model.mesh.n_points)
        self.model.mesh["rotational"][:] = np.nan
        self.model.mesh["apico-basal"][ids] = grid["apico-basal"]
        self.model.mesh["transmural"][ids] = grid["transmural"]
        self.model.mesh["rotational"][ids] = grid["rotational"]

        return grid

    def compute_right_atrial_fiber(self, appendage: List[float]) -> pv.UnstructuredGrid:
        """
        Compute right atrium fiber with LDRBD method.

        Parameters
        ----------
        appendage
            Right atrium appendage apex coordinates.

        Returns
        -------
            right atrium UnstructuredGrid with related information.

        """
        LOGGER.info("Computing RA fiber...")
        export_directory = os.path.join(self.root_directory, "ra_fiber")

        target = self.run_laplace_problem(export_directory, "ra_fiber", raa=np.array(appendage))

        endo_surface = self.model.right_atrium.endocardium
        ra_pv = compute_ra_fiber_cs(
            export_directory, self.settings.atrial_fibers, endo_surface=endo_surface
        )
        LOGGER.info("Generating fibers done.")

        # arrays that save ID map to full model
        ra_pv["cell_ids"] = target["cell_ids"]
        ra_pv["point_ids"] = target["point_ids"]

        LOGGER.info("Assigning fibers to full model...")

        # cell IDs in full model mesh
        ids = ra_pv["cell_ids"]
        self.model.mesh.cell_data["fiber"][ids] = ra_pv["e_l"]
        self.model.mesh.cell_data["sheet"][ids] = ra_pv["e_t"]

        return ra_pv

    def compute_left_atrial_fiber(self, appendage: List[float] = None) -> pv.UnstructuredGrid:
        """
        Compute left atrium fiber with LDRBD method.

        Returns
        -------
            right atrium UnstructuredGrid with related information.

        """
        LOGGER.info("Computing LA fiber...")
        export_directory = os.path.join(self.root_directory, "la_fiber")

        target = self.run_laplace_problem(export_directory, "la_fiber", laa=appendage)

        endo_surface = self.model.left_atrium.endocardium
        la_pv = compute_la_fiber_cs(
            export_directory, self.settings.atrial_fibers, endo_surface=endo_surface
        )

        LOGGER.info("Generating fibers done.")

        # arrays that save ID map to full model
        la_pv["cell_ids"] = target["cell_ids"]
        la_pv["point_ids"] = target["point_ids"]

        LOGGER.info("Assigning fibers to full model...")

        # cell IDs in full model mesh
        ids = la_pv["cell_ids"]
        self.model.mesh.cell_data["fiber"][ids] = la_pv["e_l"]
        self.model.mesh.cell_data["sheet"][ids] = la_pv["e_t"]

        return la_pv

    def run_laplace_problem(
        self, export_directory, type: Literal["uvc", "la_fiber", "ra_fiber"], **kwargs
    ):
        """
        Run Laplace-Dirichlet (thermal) problem in LSDYNA.

        Parameters
        ----------
        export_directory: str
            LSDYNA directory
        type: str
            Simulation type.
        kwargs : dict
            Additional arguments.

        Returns
        -------
            UnstructuredGrid with array to map data back to full mesh.

        """
        if type == "ra_fiber":
            for key, value in kwargs.items():
                if key == "raa":
                    break
            dyna_writer = writers.UHCWriter(copy.deepcopy(self.model), type, raa=value)
        elif type == "la_fiber":
            for key, value in kwargs.items():
                if key == "laa":
                    break
            if value is None:
                dyna_writer = writers.UHCWriter(copy.deepcopy(self.model), type)
            else:
                dyna_writer = writers.UHCWriter(copy.deepcopy(self.model), type, laa=value)
        else:
            dyna_writer = writers.UHCWriter(copy.deepcopy(self.model), type)

        dyna_writer.update()
        dyna_writer.export(export_directory)

        input_file = os.path.join(export_directory, "main.k")
        self._run_dyna(path_to_input=input_file, options="case")

        LOGGER.info("Solving laplace-dirichlet done.")

        return dyna_writer.target

    def _run_dyna(self, path_to_input: Path, options: str = ""):
        """Run LS-DYNA with path and options.

        Parameters
        ----------
        path_to_input : Path
            Path to the LS-DYNA simulation file.
        options : str, optional
            Additional options to pass to command line, by default "".

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
        """Initialize EP Simulator."""
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

    def _simulate_conduction(self, folder_name="main-ep-onlybeams"):
        """Launch the main simulation."""
        directory = os.path.join(self.root_directory, folder_name)
        self._write_main_conduction_simulation_files(folder_name)

        LOGGER.info("Launching main EP simulation...")

        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)

        LOGGER.info("done.")

        return

    def compute_purkinje(self):
        """Compute the purkinje network."""
        directory = os.path.join(self.root_directory, "purkinjegeneration")
        self.directories["purkinjegeneration"] = directory

        self._write_purkinje_files(directory)

        LOGGER.info("Computing the Purkinje network...")

        # self.settings.save(os.path.join(directory, "simulation_settings.yml"))

        LOGGER.debug("Compute Purkinje network on 1 cpu.")
        orig_num_cpus = self.dyna_settings.num_cpus
        self.dyna_settings.num_cpus = 1

        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)

        self.dyna_settings.num_cpus = orig_num_cpus
        LOGGER.debug(f"Set number of cpus back to {orig_num_cpus}.")

        LOGGER.info("done.")

        LOGGER.info("Assign the Purkinje network to the model...")

        purkinje_k_file = os.path.join(directory, "purkinjeNetwork_001.k")
        self.model.add_purkinje_from_kfile(purkinje_k_file, "Left-purkinje")

        if not isinstance(self.model, LeftVentricle):
            purkinje_k_file = os.path.join(directory, "purkinjeNetwork_002.k")
            self.model.add_purkinje_from_kfile(purkinje_k_file, "Right-purkinje")

    def compute_conduction_system(self):
        """Compute the conduction system."""
        if isinstance(self.model, FourChamber):
            beam_length = self.settings.purkinje.edgelen.m

            cs = ConductionSystem(self.model)
            cs.compute_SA_node()
            cs.compute_AV_node()
            cs.compute_av_conduction(beam_length=beam_length)
            left, right = cs.compute_His_conduction(beam_length=beam_length)
            cs.compute_left_right_bundle(
                left.xyz, left.node_id, side="Left", beam_length=beam_length
            )
            cs.compute_left_right_bundle(
                right.xyz, right.node_id, side="Right", beam_length=beam_length
            )

            # # TODO define end point by uhc, or let user choose
            # Note: must on surface after zerop if coupled with meca
            # cs.compute_Bachman_bundle(
            #     start_coord=self.model.right_atrium.get_point("SA_node").xyz,
            #     end_coord=np.array([-34, 163, 413]),
            # )
            return cs

    def _write_main_simulation_files(self, folder_name):
        """Write LS-DYNA files that are used to start the main simulation."""
        export_directory = os.path.join(self.root_directory, folder_name)
        self.directories["main-ep"] = export_directory
        model = copy.deepcopy(self.model)
        dyna_writer = writers.ElectrophysiologyDynaWriter(model, self.settings)
        dyna_writer.update()
        dyna_writer.export(export_directory)

        return export_directory

    def _write_main_conduction_simulation_files(self, folder_name):
        """Write LS-DYNA files that are used to start the main simulation."""
        export_directory = os.path.join(self.root_directory, folder_name)
        self.directories["main-ep"] = export_directory
        model = copy.deepcopy(self.model)
        dyna_writer = writers.ElectrophysiologyBeamsDynaWriter(model, self.settings)
        dyna_writer.update()
        dyna_writer.export(export_directory)

        return export_directory

    def _write_purkinje_files(
        self,
        export_directory,
    ) -> Path:
        """Write purkinje files."""
        model = copy.deepcopy(self.model)
        dyna_writer = writers.PurkinjeGenerationDynaWriter(model, self.settings)
        dyna_writer.update()
        dyna_writer.export(export_directory)
        return


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

        self.initial_stress = initial_stress
        """If stress free computation is taken into considered."""
        # include initial stress by default

        self.stress_free_report = None
        """A dictionary save stress free computation information"""

        return

    def create_stiff_ventricle_base(
        self, threshold: float = 0.9, stiff_material=NeoHookean(rho=0.001, c10=0.1, nu=0.499)
    ) -> Part:
        """Create a stiff base part from uvc longitudinal value.

        Parameters
        ----------
        threshold : float, optional
            uvc_l larger than threshold will be set as stiff base, by default 0.9
        stiff_material : _type_, optional
            material to assign, by default NeoHookean(rho=0.001, c10=0.1, nu=0.499)

        Returns
        -------
        Part
            stiff base part
        """
        try:
            v = self.model.mesh.point_data_to_cell_data()["apico-basal"]
        except:
            self.compute_uhc()
            v = self.model.mesh.point_data_to_cell_data()["apico-basal"]

        eids = np.intersect1d(np.where(v > threshold)[0], self.model.left_ventricle.element_ids)
        if not isinstance(self.model, LeftVentricle):
            # uvc-L of RV is generally smaller, *1.05 to be comparable with LV
            eid_r = np.intersect1d(
                np.where(v > threshold * 1.05)[0], self.model.right_ventricle.element_ids
            )
            eids = np.hstack((eids, eid_r))

        part: Part = self.model.create_part_by_ids(eids, "base")
        part.part_type = "ventricle"
        part.has_fiber = False
        part.is_active = False
        part.meca_material = stiff_material

        return part

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
                dynain_file = os.path.join(zerop_folder, "iter3.dynain.lsda")

                shutil.copy(dynain_file, os.path.join(directory, "dynain.lsda"))
                shutil.copy(
                    os.path.join(zerop_folder, "post", "Post_report.json"),
                    os.path.join(directory, "Post_report.json"),
                )
            except IndexError:
                # handle if lsda file not exist.
                LOGGER.warning("Cannot find initial stress file iter3.dynain.lsda")
                exit()

        self._write_main_simulation_files(folder_name=folder_name)

        LOGGER.info("Launching main simulation...")

        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)

        LOGGER.info("done.")

        if auto_post:
            mech_post(Path.Path(directory), self.model)
        return

    def compute_stress_free_configuration(self, folder_name="zeropressure", overwrite: bool = True):
        """Compute the stress-free configuration of the model."""
        directory = os.path.join(self.root_directory, folder_name)
        os.makedirs(
            directory,
            exist_ok=True,
        )

        if overwrite or len(os.listdir(directory)) == 0:
            self._write_stress_free_configuration_files(folder_name)
            self.settings.save(Path.Path(directory) / "simulation_settings.yml")

            LOGGER.info("Computing stress-free configuration...")
            self._run_dyna(os.path.join(directory, "main.k"), options="case")
            LOGGER.info("Simulation done.")
        else:
            LOGGER.info(f"Re-using existing results in {directory}")

        self.stress_free_report = zerop_post(directory, self.model)

        # replace node coordinates by computed ED geometry
        LOGGER.info("Updating nodes after stress-free.")

        # Note: cap center node will be added into mesh.points
        if heart_version == "v0.1":
            n_caps = len(self.model.cap_centroids)
            guess_ed_coords = np.array(self.stress_free_report["guess_ed_coord"])[:-n_caps]
        elif heart_version == "v0.2":
            guess_ed_coords = np.array(self.stress_free_report["guess_ed_coord"])
        self.model.mesh.nodes = guess_ed_coords

        # Note: synchronization
        for part in self.model.parts:
            for surface in part.surfaces:
                surface.nodes = guess_ed_coords

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

        model = copy.deepcopy(self.model)
        if isinstance(model, FourChamber) and type(self) == EPMechanicsSimulator:
            model._create_isolation_part()

        dyna_writer = writers.ZeroPressureMechanicsDynaWriter(model, self.settings)
        dyna_writer.update()
        dyna_writer.export(export_directory)

        return


class EPMechanicsSimulator(EPSimulator, MechanicsSimulator):
    """Coupled EP-mechanics simulator with computed Electrophysiology."""

    def __init__(
        self,
        model: HeartModel,
        dyna_settings: DynaSettings,
        simulation_directory: Path = "",
    ) -> None:
        MechanicsSimulator.__init__(self, model, dyna_settings, simulation_directory)

        return

    def simulate(self, folder_name="ep_meca"):
        """Launch the main simulation."""
        # MechanicalSimulator handle dynain file from zerop
        MechanicsSimulator.simulate(self, folder_name=folder_name)

        return

    def _write_main_simulation_files(self, folder_name):
        """Write LS-DYNA files that are used to start the main simulation."""
        export_directory = os.path.join(self.root_directory, folder_name)
        self.directories["main-coupling"] = export_directory

        dyna_writer = writers.ElectroMechanicsDynaWriter(self.model, self.settings)
        dyna_writer.update(with_dynain=self.initial_stress)
        dyna_writer.export(export_directory)

        return export_directory


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
        LS-DYNA settings, such as path to executable, executable type,
        platform, by default ``None``.
    simulation_directory : Path, optional
        Directory where to simulate, by default ``None``.

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
