import subprocess
import shutil
from unittest import case
import numpy as np
import pandas as pd
import os
import time
import json
from pathlib import Path
from tqdm import tqdm  # for progress bar
from ansys.heart.writer.dynawriter import *

# class purkinje:
#     """define purkinje network parameters"""

#     def __init__(self, origin_coordinates: np.ndarray[shape=(1,3)]) -> None:
#         self.origin = origin_coordinates


class Simulator:
    """perform pre-simulation steps prior to final simulation, steps include fiber generation,
     zeo-pressure geometry and purkinje network generation"""

    def __init__(self, model: HeartModel, lsdynapath: str = "",) -> None:

        """Initializes Simulator by loading a HeartModel

        Parameters
        ----------
        model : HeartModel
            HeartModel object which contains the necessary
            information for the wowrkflow, such as nodes and elements.
        """

        self.model = model
        """Heart model"""
        self.lsdynapath = lsdynapath
        """path of the lsdyna executable"""

    def build(
        self,
        fibers: bool = 0,
        purkinje: bool = 0,
        zeropressure: bool = 0,
        ep: bool = 0,
        mechanics: bool = 0,
    ):
        # TODO: add a "simulation" folder that contains all the precomputed features
        # TODO: make sure the necessary node, segment, part ids consistant in the different steps
        # TODO: make sure the necessary node, segment, part ids consistant in the different steps
        # TODO add getters and setters for fiber angles, purkinje properties, simulation times and
        # other parameters to expose to the user
        
        if fibers:
        # Generate Fibers
            self.write_fibers()

        if zeropressure:
            
            shutil.copy2(
                os.path.join("fibergeneration", "element_solid_ortho.k"),
                os.path.join("mechanics", "solid_elements.k"),
            )
        # TODO copy element_solid_ortho.k to "simulation" folder

        # Generate ZeroP
        self.write_zeropressureconfiguration()
        # TODO handle guess files and nodes.k

        # Generate Purkinje
        self.write_purkinje()
        # TODO rename purkinje files and move them to the Sim. folder

        # TODO continue in the same spirit with zeroP, EP, mechanics, and Purkinje (in the future, add blood, and fluid)

        if (ep) and not (mechanics):
            LOGGER.debug("Pue EP model")

        if not (ep) and (mechanics):
            LOGGER.debug("Pure Mechanical model")
        if ep and mechanics:
            # TODO add coupling stuff and ignore default Ca2+ mechanics active stress/replaced by EP simulation
            LOGGER.debug("Coupled EP + Mechanics problem")

        return

    def write_fibers(
        self,
        workdir: str,
        alpha_endocardium: float = -60,
        alpha_eepicardium: float = 60,
        beta_endocardium: float = 25,
        beta_epicardium: float = -65,
    ):
        dyna_writer = FiberGenerationDynaWriter(self.model)
        dyna_writer.update()
        dyna_writer.export(workdir)

        return

    def write_zeropressureconfiguration(self, workdir: str = ""):
        dyna_writer = ZeroPressureMechanicsDynaWriter(self.model)
        dyna_writer.update()
        dyna_writer.export(workdir)
        return

    def generatePurkinje(
        self,
        workdir: str,
        pointstx: float = 0,  # TODO instanciate this
        pointsty: float = 0,  # TODO instanciate this
        pointstz: float = 0,  # TODO instanciate this
        inodeid: int = 0,  # TODO instanciate this
        iedgeid: int = 0,  # TODO instanciate this
        edgelen: float = 2,  # TODO instanciate this
        ngen: float = 50,
        nbrinit: int = 8,
        nsplit: int = 2,
    ):
        dyna_writer = PurkinjeGenerationDynaWriter(self.model)
        dyna_writer.update()
        dyna_writer.export(workdir)
        return

    def run_lsdyna(sim_file: str, memory: int, lsdynapath: str, NCPU: int = 1, options: str = ""):
        """
        Parameters
        ----------
        sim_file: input file for lsdyna
        NCPU: number of CPUs
        memory: memory usage (ex. 300M for large models)
        options: lsdyna additional run options

        Returns
        -------

        """
        commands = [
            "mpirun",
            "-np",
            str(NCPU),
            lsdynapath,
            "i=" + sim_file,
            "memory=",
            str(memory),
            options,
        ]

        p = subprocess.run(commands, stdout=subprocess.PIPE)
        return
