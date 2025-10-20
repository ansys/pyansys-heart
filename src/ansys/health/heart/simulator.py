# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
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

  - With/without fibers
  - With/without Purkinje

- Electro-mechanics

  - Simplified EP (imposed activation)
  - Coupled electro-mechanics
"""

import copy
import glob
import os
import pathlib
import shutil
import subprocess
from typing import Literal

import natsort
import numpy as np
import psutil
import pyvista as pv

from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.exceptions import LSDYNANotFoundError, LSDYNATerminationError
import ansys.health.heart.models as models
from ansys.health.heart.models_utils import HeartModelUtils
from ansys.health.heart.objects import SurfaceMesh
from ansys.health.heart.post.auto_process import mech_post, zerop_post
from ansys.health.heart.post.laplace_post import (
    compute_la_fiber_cs,
    compute_ra_fiber_cs,
    compute_ventricle_fiber_by_drbm,
    read_laplace_solution,
)
from ansys.health.heart.pre.conduction_path import ConductionPath, ConductionPathType
from ansys.health.heart.settings.settings import DynaSettings, SimulationSettings
from ansys.health.heart.utils.misc import _read_orth_element_kfile
import ansys.health.heart.writer as writers

_KILL_ANSYSCL_PRIOR_TO_RUN = True
"""Flag indicating whether to kill all Ansys license clients prior to an LS-DYNA run."""


class BaseSimulator:
    """Base class for the simulator."""

    def __init__(
        self,
        model: models.FullHeart | models.FourChamber | models.BiVentricle | models.LeftVentricle,
        dyna_settings: DynaSettings = None,
        simulation_directory: pathlib = "",
    ) -> None:
        """Initialize BaseSimulator.

        Parameters
        ----------
        model : HeartModel
            Heart model to simulate.
        dyna_settings : DynaSettings
            Settings used for launching LS-DYNA.
        simulation_directory : Path, default: ""
            Directory to start the simulation in.

        """
        self.model = model
        """Heart model to simulate."""
        if not dyna_settings:
            LOGGER.warning("Setting default LS-DYNA settings.")
            self.dyna_settings = DynaSettings()
        else:
            self.dyna_settings: DynaSettings = dyna_settings
            """Contains the settings to launch LS-DYNA."""

        if self.dyna_settings.platform != "wsl":
            if shutil.which(self.dyna_settings.lsdyna_path) is None:
                LOGGER.error(f"{self.dyna_settings.lsdyna_path} does not exist.")
                raise LSDYNANotFoundError(
                    f"LS-DYNA executable {self.dyna_settings.lsdyna_path} file is not found."
                )

        if simulation_directory == "":
            simulation_directory = os.path.join(self.model.workdir, "simulation")

        self.root_directory = simulation_directory
        """Root simulation directory."""

        self.settings: SimulationSettings = SimulationSettings()
        """Simulation settings."""

        pass

    def load_default_settings(self) -> SimulationSettings:
        """Load default simulation settings."""
        self.settings.load_defaults()
        return self.settings

    def compute_fibers(
        self, method: Literal["LSDYNA", "D-RBM"] = "LSDYNA", rotation_angles: dict = None
    ):
        """Compute the fiber sheet directions on the ventricles.

        Parameters
        ----------
        method : Literal["LSDYNA", "D-RBM"], default: "LSDYNA"
            Method for computing the fiber orientation.
        rotation_angles : dict, default: None
            Rotation angle alpha and beta.
        """
        LOGGER.info("Computing fiber orientation...")

        if method == "LSDYNA":
            if rotation_angles is None:
                # find default settings
                rotation_angles = self.settings.get_ventricle_fiber_rotation(method="LSDYNA")

            for name in ["alpha", "beta", "beta_septum"]:
                if name not in rotation_angles.keys():
                    LOGGER.error(f"Must provide key {name} for D-RBM method.")
                    raise KeyError(
                        f"Must provide key {name} for D-RBM method. "
                        "Please check the settings or provide the rotation angles."
                    )

            self._compute_fibers_lsdyna(rotation_angles)

        elif method == "D-RBM":
            if rotation_angles is None:
                # find default settings
                rotation_angles = self.settings.get_ventricle_fiber_rotation(method="D-RBM")

            for a, b in zip(["alpha", "beta"], ["_left", "_right", "_ot"]):
                if a + b not in rotation_angles.keys():
                    LOGGER.error(f"Must provide key {name} for D-RBM method.")
                    raise KeyError(
                        f"Must provide key {name} for D-RBM method. "
                        "Please check the settings or provide the rotation angles."
                    )
            self._compute_fibers_drbm(rotation_angles)

        else:
            LOGGER.error(f"Method {method} is not recognized.")
            raise ValueError(
                f"Method {method} is not recognized. "
                "Please use 'LSDYNA' or 'D-RBM' as the method for computing fibers."
            )

        return

    def _compute_fibers_drbm(self, rotation_angles: dict):
        """Use D-RBM fiber method."""
        export_directory = os.path.join(self.root_directory, "D-RBM")
        target = self.run_laplace_problem(export_directory, type="D-RBM")
        grid = compute_ventricle_fiber_by_drbm(
            export_directory,
            settings=rotation_angles,
            left_only=isinstance(self.model, models.LeftVentricle),
        )
        grid.save(os.path.join(export_directory, "drbm_fibers.vtu"))

        # arrays that save ID map to full model
        grid["cell_ids"] = target["cell_ids"]

        LOGGER.info("Assigning fibers to full model...")

        # cell IDs in full model mesh
        ids = grid["cell_ids"]
        self.model.mesh.cell_data["fiber"][ids] = grid["fiber"]
        self.model.mesh.cell_data["sheet"][ids] = grid["sheet"]

    def _compute_fibers_lsdyna(self, rotation_angles: dict):
        """Use LSDYNA native fiber method."""
        directory = os.path.join(self.root_directory, "fibergeneration")

        dyna_writer = writers.FiberGenerationDynaWriter(copy.deepcopy(self.model), self.settings)
        dyna_writer.update(rotation_angles)
        dyna_writer.export(directory)

        input_file = os.path.join(directory, "main.k")
        self._run_dyna(path_to_input=input_file)

        # TODO: May want to replace by ansys.dyna.core.keywords
        LOGGER.info("Assigning fiber orientation to model...")
        elem_ids, part_ids, connect, fib, sheet = _read_orth_element_kfile(
            os.path.join(directory, "element_solid_ortho.k")
        )
        self.model.mesh.cell_data["fiber"][elem_ids - 1] = fib
        self.model.mesh.cell_data["sheet"][elem_ids - 1] = sheet

        return

    def compute_uhc(self) -> pv.UnstructuredGrid:
        """Compute universal heart coordinates system."""
        LOGGER.info("Computing universal ventricular coordinates...")

        type = "uvc"
        export_directory = os.path.join(self.root_directory, type)

        target = self.run_laplace_problem(export_directory, type)
        grid = read_laplace_solution(
            export_directory, field_list=["apico-basal", "transmural", "rotational"]
        )

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

    def compute_right_atrial_fiber(
        self, appendage: list[float], top: list[list[float]] = None
    ) -> pv.UnstructuredGrid:
        """
        Compute right atrium fiber with the LDRBD method.

        Parameters
        ----------
        appendage: list[float]
            Coordinates of appendage.
        top : list[list[float]], default: None
            List of nodal coordinates to define the top path.

            The top path is a set of nodes connecting the superior (SVC) and inferior (IVC)
            vena cava. For more information, see the "Notes" section.

            The default method (``top=None``) might not work for some anatomical structures.
            In such cases, you can define the start and end points by providing a list of
            coordinates like this: ``[[x1, y1, z1], [x2, y2, z2]]``. These two nodes should
            be located on the SVC and IVC rings, approximately at the 12 o'clock position.

            You can also add an intermediate point to enforce the geodesic path, like this:
            ``[[x1, y1, z1], [x3, y3, z3], [x2, y2, z2]]``.

        Returns
        -------
        pv.UnstructuredGrid
            Right atrium with fiber coordinates system in this format: ``e_l``, ``e_t`` and ``e_n``.

        Notes
        -----
        The method is described in `Modeling cardiac muscle fibers in ventricular and atrial
        electrophysiology simulations <https://doi.org/10.1016/j.cma.2020.113468>`_
        """
        LOGGER.info("Computing right atrium fiber...")
        export_directory = os.path.join(self.root_directory, "ra_fiber")

        target = self.run_laplace_problem(
            export_directory, "ra_fiber", raa=np.array(appendage), top=top
        )

        ra_pv = compute_ra_fiber_cs(
            export_directory, self.settings.atrial_fibers, endo_surface=None
        )
        ra_pv.save(os.path.join(export_directory, "ra_fiber.vtu"))
        LOGGER.info("Generating fibers is done.")

        # arrays that save ID map to full model
        ra_pv["cell_ids"] = target["cell_ids"]
        ra_pv["point_ids"] = target["point_ids"]

        LOGGER.info("Assigning fibers to full model...")

        # cell IDs in full model mesh
        ids = ra_pv["cell_ids"]
        self.model.mesh.cell_data["fiber"][ids] = ra_pv["e_l"]
        self.model.mesh.cell_data["sheet"][ids] = ra_pv["e_t"]

        return ra_pv

    def compute_left_atrial_fiber(
        self,
        appendage: list[float] = None,
    ) -> pv.UnstructuredGrid:
        """Compute left atrium fiber with the LDRBD method.

        Parameters
        ----------
        appendage : list[float], default: None
            Coordinates of the appendage. If no value is specified,
            the cap named ``appendage`` is used.

        Returns
        -------
        pv.UnstructuredGrid
            Left atrium with fiber coordinates system in this format: ``e_l``, ``e_t`` and ``e_n``.

        Notes
        -----
        The method is described in `Modeling cardiac muscle fibers in ventricular and atrial
        electrophysiology simulations <https://doi.org/10.1016/j.cma.2020.113468>`.
        """
        LOGGER.info("Computing left atrium fiber...")
        export_directory = os.path.join(self.root_directory, "la_fiber")

        target = self.run_laplace_problem(export_directory, "la_fiber", laa=appendage)

        la_pv = compute_la_fiber_cs(
            export_directory, self.settings.atrial_fibers, endo_surface=None
        )
        la_pv.save(os.path.join(export_directory, "la_fiber.vtu"))
        LOGGER.info("Generating fibers is done.")

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
        Run the Laplace-Dirichlet (thermal) problem in LS-DYNA.

        Parameters
        ----------
        export_directory: str
            LS-DYNA directory
        type: str
            Simulation type.
        kwargs : dict
            Landmarks to create the nodeset. Keys can be ``laa``, ``raa``, and ``top``'.

        Returns
        -------
        UnstructuredGrid
            UnstructuredGrid with array to map data back to the full mesh.

        """
        for k, v in kwargs.items():
            if k not in ["laa", "raa", "top"]:
                LOGGER.error(f"kwarg with {k} can not be identified.")
                raise KeyError(f"kwarg with {k} can not be identified.")

        kwargs = {k: v for k, v in kwargs.items() if v is not None}

        dyna_writer = writers.LaplaceWriter(copy.deepcopy(self.model), type, **kwargs)

        dyna_writer.update()
        dyna_writer.export(export_directory)

        input_file = os.path.join(export_directory, "main.k")
        self._run_dyna(path_to_input=input_file, options="case")

        LOGGER.info("Solving Laplace-Dirichlet problem is done.")

        return dyna_writer.target

    def _run_dyna(self, path_to_input: pathlib, options: str = ""):
        """Run LS-DYNA with the specified input file and options.

        Parameters
        ----------
        path_to_input : Path
            Path to the LS-DYNA simulation file.
        options : str, default: ""
            Additional options to pass to the command line.

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


class EPSimulator(BaseSimulator):
    """EP (electrophysiology) simulator."""

    def __init__(
        self,
        model: models.FullHeart | models.FourChamber | models.BiVentricle | models.LeftVentricle,
        dyna_settings: DynaSettings,
        simulation_directory: pathlib = "",
    ) -> None:
        """Initialize the EP simulator."""
        super().__init__(model, dyna_settings, simulation_directory)

        return

    def simulate(self, folder_name="main-ep", extra_k_files: list[str] | None = None):
        """Launch the EP simulation.

        Parameters
        ----------
        folder_name : str, default: ``'main-ep'``
            Simulation folder name.
        extra_k_files : list[str], default: None
            User-defined k files.
        """
        directory = os.path.join(self.root_directory, folder_name)
        self._write_main_simulation_files(folder_name, extra_k_files=extra_k_files)

        LOGGER.info("Launching main EP simulation...")

        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)

        LOGGER.info("Simulation completed successfully.")

        return

    def _simulate_conduction(self, folder_name="main-ep-onlybeams"):
        """Launch the main EP simulation."""
        directory = os.path.join(self.root_directory, folder_name)
        self._write_main_conduction_simulation_files(folder_name)

        LOGGER.info("Launching main EP simulation...")

        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)

        LOGGER.info("Simulation completed successfully.")

        return

    def compute_purkinje(self):
        """Compute the Purkinje network."""
        directory = os.path.join(self.root_directory, "purkinjegeneration")

        self._write_purkinje_files(directory)

        LOGGER.info("Computing the Purkinje network...")

        # self.settings.save(os.path.join(directory, "simulation_settings.yml"))

        LOGGER.debug("Compute Purkinje network on one CPU.")
        orig_num_cpus = self.dyna_settings.num_cpus
        self.dyna_settings.num_cpus = 1

        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)
        LOGGER.info("Simulation completed successfully.")

        self.dyna_settings.num_cpus = orig_num_cpus
        LOGGER.debug(f"Set number of CPUs back to {orig_num_cpus}.")

        LOGGER.info("Assign the Purkinje network to the model...")

        left_purkinje = ConductionPath.create_from_k_file(
            ConductionPathType.LEFT_PURKINJE,
            k_file=os.path.join(directory, "purkinjeNetwork_001.k"),
            id=1,
            base_mesh=self.model.left_ventricle.endocardium,
            model=self.model,
        )

        if isinstance(self.model, models.LeftVentricle):
            self.model.assign_conduction_paths([left_purkinje])
            return left_purkinje
        else:
            right_purkinje = ConductionPath.create_from_k_file(
                ConductionPathType.RIGHT_PURKINJE,
                k_file=os.path.join(directory, "purkinjeNetwork_002.k"),
                id=2,
                base_mesh=self.model.right_ventricle.endocardium,
                model=self.model,
            )
            self.model.assign_conduction_paths([left_purkinje, right_purkinje])
            return [left_purkinje, right_purkinje]

    def compute_conduction_system(self):
        """Compute the conduction system."""
        if isinstance(self.model, models.FourChamber):
            # TODO: refinement is not correctly used
            # beam_length = self.settings.purkinje.edgelen.m

            beam_list = HeartModelUtils.define_full_conduction_system(
                self.model, os.path.join(self.root_directory, "purkinjegeneration")
            )
            self.model.assign_conduction_paths(beam_list)
        else:
            LOGGER.info("Computation is only implemented for other than FourChamber models.")

        return beam_list

    def _write_main_simulation_files(self, folder_name, extra_k_files: list[str] | None = None):
        """Write LS-DYNA files that are used to start the main EP simulation."""
        export_directory = os.path.join(self.root_directory, folder_name)

        model = copy.deepcopy(self.model)
        dyna_writer = writers.ElectrophysiologyDynaWriter(model, self.settings)
        dyna_writer.update()
        dyna_writer.export(export_directory, user_k=extra_k_files)

        return

    def _write_main_conduction_simulation_files(self, folder_name):
        """Write LS-DYNA files that are used to start the main EP simulation."""
        export_directory = os.path.join(self.root_directory, folder_name)

        model = copy.deepcopy(self.model)
        dyna_writer = writers.ElectrophysiologyBeamsDynaWriter(model, self.settings)
        dyna_writer.update()
        dyna_writer.export(export_directory)

        return

    def _write_purkinje_files(
        self,
        export_directory,
    ) -> pathlib:
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
        model: models.FullHeart | models.FourChamber | models.BiVentricle | models.LeftVentricle,
        dyna_settings: DynaSettings,
        simulation_directory: pathlib = "",
        initial_stress: bool = True,
    ) -> None:
        super().__init__(model, dyna_settings, simulation_directory)

        self.initial_stress = initial_stress
        """If stress-free computation is taken into consideration."""
        self._dynain_name = None
        """LS-DYNA initial state file name from zeropressure."""
        return

    def simulate(
        self,
        folder_name: str = "main-mechanics",
        zerop_folder: str | None = None,
        auto_post: bool = True,
        extra_k_files: list[str] | None = None,
    ):
        """Launch the main mechanical simulation.

        Parameters
        ----------
        folder_name : str, default: ``'main-mechanics'``
            Simulation folder name.
        zerop_folder : str | None, default: None
            Folder containing stress-free simulation.
            If ``None``, the ``zeropressure`` folder under the root directory is used.
        auto_post : bool, default: True
            Whether to run postprocessing scripts.
        extra_k_files : list[str], default: None
            User-defined k files.
        """
        if "apico-basal" not in self.model.mesh.point_data.keys():
            LOGGER.warning(
                "Array named ``apico-basal`` cannot be found. Computing"
                "universal coordinate system (UVC) first."
            )
            self.compute_uhc()

        directory = os.path.join(self.root_directory, folder_name)
        os.makedirs(directory, exist_ok=True)

        if self.initial_stress:
            dynain_file = self._find_dynain_file(zerop_folder)
            self._dynain_name = "dynain.lsda"
            shutil.copy(dynain_file, os.path.join(directory, self._dynain_name))

        self._write_main_simulation_files(folder_name=folder_name, extra_k_files=extra_k_files)

        LOGGER.info("Launching main simulation...")

        input_file = os.path.join(directory, "main.k")
        self._run_dyna(input_file)

        LOGGER.info("done.")

        if auto_post:
            mech_post(pathlib.Path(directory), self.model)

        return

    def _find_dynain_file(self, zerop_folder) -> str:
        """Find the ``dynain.lsda`` file of the last iteration."""
        if zerop_folder is None:
            zerop_folder = os.path.join(self.root_directory, "zeropressure")

        dynain_files = glob.glob(os.path.join(zerop_folder, "iter*.dynain.lsda"))
        # force natural ordering since iteration numbers are not padded with zeros.
        dynain_files = natsort.natsorted(dynain_files)

        if len(dynain_files) == 0:
            error_message = f"Files 'iter*.dynain.lsda` not found in {zerop_folder}."
            LOGGER.error(error_message)
            raise FileNotFoundError(error_message)

        elif len(dynain_files) == 1:
            error_message = (
                f"Only 1 'iter*.dynain.lsda' file is found in {zerop_folder}. Expected at least 2."
            )

            LOGGER.error(error_message)
            raise IndexError(error_message)

        else:
            dynain_file = dynain_files[-1]
            LOGGER.info(f"Using {dynain_file} for initial stress.")

        return dynain_file

    def compute_stress_free_configuration(
        self,
        folder_name="zeropressure",
        overwrite: bool = True,
        extra_k_files: list[str] | None = None,
    ) -> tuple[dict, np.ndarray, np.ndarray]:
        """Compute the stress-free configuration of the model.

        Parameters
        ----------
        folder_name : str, default: ``'zeropressure'``
            Simulation folder name.
        overwrite : bool, default: True
            Whether to run simulation and overwrite files.
        extra_k_files : list[str], default: None
            User-defined k files.

        Returns
        -------
        tuple[dict, np.ndarray, np.ndarray]
            Dictionary with convergence information,
            stress free configuration, and
            (re)computed end-of-diastole configuration.
        """
        directory = os.path.join(self.root_directory, folder_name)

        if not os.path.isdir(directory) or overwrite or len(os.listdir(directory)) == 0:
            os.makedirs(directory, exist_ok=True)

            self._write_stress_free_configuration_files(folder_name, extra_k_files=extra_k_files)
            self.settings.save(pathlib.Path(directory) / "simulation_settings.yml")

            LOGGER.info("Computing stress-free configuration...")
            self._run_dyna(os.path.join(directory, "main.k"), options="case")
            LOGGER.info("Simulation is done.")
        else:
            LOGGER.info(f"Reusing existing results in {directory}.")

        report, stress_free_coord, guess_ed_coord = zerop_post(directory, self.model)

        # replace node coordinates by computed ED geometry
        LOGGER.info("Updating nodes...")

        self.model.mesh.points = guess_ed_coord

        # Update surfaces of all parts
        # TODO: move it into Part or Model
        for part in self.model.parts:
            for key, value in part.__dict__.items():
                if isinstance(value, SurfaceMesh):
                    part.__setattr__(key, self.model.mesh.get_surface(value.id))

        return report, stress_free_coord, guess_ed_coord

    def _write_main_simulation_files(
        self,
        folder_name,
        extra_k_files: list[str] | None = None,
    ):
        """Write LS-DYNA files that are used to start the main simulation."""
        export_directory = os.path.join(self.root_directory, folder_name)

        dyna_writer = writers.MechanicsDynaWriter(
            self.model,
            self.settings,
        )
        dyna_writer.update(dynain_name=self._dynain_name)
        dyna_writer.export(export_directory, user_k=extra_k_files)

        return

    def _write_stress_free_configuration_files(
        self, folder_name, extra_k_files: list[str] | None = None
    ):
        """Write LS-DYNA files to compute the stress-free configuration."""
        export_directory = os.path.join(self.root_directory, folder_name)

        model = copy.deepcopy(self.model)
        # Isolation part need to be created in Zerop because main will use its dynain.lsda
        if isinstance(model, models.FourChamber) and isinstance(self, EPMechanicsSimulator):
            model._create_atrioventricular_isolation()

        dyna_writer = writers.ZeroPressureMechanicsDynaWriter(model, self.settings)
        dyna_writer.update()
        dyna_writer.export(export_directory, user_k=extra_k_files)

        return


class EPMechanicsSimulator(EPSimulator, MechanicsSimulator):
    """Coupled EP-mechanics simulator with computed electrophysiology."""

    def __init__(
        self,
        model: models.FullHeart | models.FourChamber | models.BiVentricle | models.LeftVentricle,
        dyna_settings: DynaSettings,
        simulation_directory: pathlib = "",
    ) -> None:
        MechanicsSimulator.__init__(self, model, dyna_settings, simulation_directory)

        return

    def simulate(
        self,
        folder_name: str = "ep_meca",
        zerop_folder: str | None = None,
        auto_post: bool = True,
        extra_k_files: list[str] | None = None,
    ):
        """Launch the main electro-mechanical simulation.

        Parameters
        ----------
        folder_name : str, default: ``'main-mechanics'``
            Simulation folder name.
        zerop_folder : str | None, default: None
            Folder containing stress-free simulation.
            If ``None``, the ``zeropressure`` folder under the root directory is used.
        auto_post : bool, default: True
            Whether to run postprocessing scripts.
        extra_k_files : list[str], default: None
            User-defined k files.
        """
        # MechanicalSimulator handle dynain file from zerop
        MechanicsSimulator.simulate(
            self,
            folder_name=folder_name,
            zerop_folder=zerop_folder,
            auto_post=auto_post,
            extra_k_files=extra_k_files,
        )

        return

    def _write_main_simulation_files(
        self,
        folder_name,
        extra_k_files: list[str] | None = None,
    ):
        """Write LS-DYNA files that are used to start the main simulation."""
        export_directory = os.path.join(self.root_directory, folder_name)

        dyna_writer = writers.ElectroMechanicsDynaWriter(self.model, self.settings)
        dyna_writer.update(dynain_name=self._dynain_name)
        dyna_writer.export(export_directory, user_k=extra_k_files)

        return


def _kill_all_ansyscl():
    """Kill all Ansys license clients."""
    try:
        for p in psutil.process_iter():
            if "ansyscl" in p.name():
                p.kill()
    except Exception as e:
        LOGGER.warning(f"Failed to kill all ansyscl's: {e}")


def run_lsdyna(
    path_to_input: pathlib,
    settings: DynaSettings = None,
    simulation_directory: pathlib = None,
):
    """Standalone function for running LS-DYNA.

    Parameters
    ----------
    path_to_input : Path
        Input file for LS-DYNA.
    settings : DynaSettings, default: None
        LS-DYNA settings, such as path to the executable file, executable type,
        and platform.
    simulation_directory : Path, default: None
        Directory for the simulation.

    """
    if not settings:
        LOGGER.info("Using default LS-DYNA settings.")
        raise ValueError("Settings must be provided.")

    commands = settings.get_commands(path_to_input)

    os.chdir(os.path.dirname(path_to_input))

    #! Kill all Ansys license clients prior to running LS-DYNA
    #! this to avoid issues with orphan license clients of versions
    #! lower than the one needed by LS-DYNA.
    if _KILL_ANSYSCL_PRIOR_TO_RUN:
        _kill_all_ansyscl()

    LOGGER.info(f"Starting LS-DYNA with {commands}...")

    mess = []
    with subprocess.Popen(commands, stdout=subprocess.PIPE, text=True) as p:
        for line in p.stdout:
            LOGGER.debug(line.rstrip())
            mess.append(line)

    os.chdir(simulation_directory)

    if "N o r m a l    t e r m i n a t i o n" not in "".join(mess):
        if "numNodePurkinje" not in "".join(mess):
            LOGGER.error("LS-DYNA did not terminate properly.")
            raise LSDYNATerminationError(mess)

    LOGGER.info("LS-DYNA simulation successful...")

    return
