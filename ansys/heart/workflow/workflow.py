import subprocess
import numpy as np
import pandas as pd
import os
import time
import json
from pathlib import Path
from tqdm import tqdm  # for progress bar
from ansys.heart.writer.dynawriter import *

class purkinje:
    """define purkinje network parameters"""

    def __init__(self, origin_coordinates: np.ndarray[shape=(1,3)]) -> None:
        self.origin = origin_coordinates


class Workflow:
    """perform pre-simulation steps prior to final simulation, steps include fiber generation,
     zeo-pressure geometry and purkinje network generation"""

    def __init__(
        self, 
        model: HeartModel,
        lsdynapath: str="",
        ) -> None:

        """Initializes workflow by loading a HeartModel

        Parameters
        ----------
        model : HeartModel
            HeartModel object which contains the necessary
            information for the wowrkflow, such as nodes and elements.
        """

        self.model = model
        self.lsdynapath = lsdynapath

    def build(
        self,
        Fibers:bool=0,
        Purkinje:bool=0,
        ZeroPressure:bool=0,
        EP:bool=0,
        Mechanics:bool=0,  
    ):
        # TODO: add a "simulation" folder that contains all the precomputed features
        # TODO: make sure the necessary node, segment, part ids consistant in the different steps
        # TODO: make sure the necessary node, segment, part ids consistant in the different steps
        # TODO add getters and setters for fiber angles, purkinje properties, simulation times and other parameters to expose to the user
        if Fibers:
            LOGGER.debug("Generating fibers")
            self.generateFibers()
            # TODO copy element_solid_ortho.k to "simulation" folder
        # TODO continue in the same spirit with zeroP, EP, mechanics, and Purkinje (in the future, add blood, and fluid)



        if EP and Mechanics:
            # TODO add coupling stuff and ignore default Ca2+ mechanics active stress/replaced by EP simulation
            LOGGER.debug("Add Electromechanics coupling")

        return

    def generateZeroPressureConfiguration(self,export_directory:str=""):
        dyna_writer = ZeroPressureMechanicsDynaWriter(self.model)
        dyna_writer.update()
        if not export_directory:
            export_directory = os.path.join(self.model.info.working_directory,"ZeroPressure".lower())
        dyna_writer.export(export_directory)
        sim_file=os.path.join(export_directory,"main.k")
        self.run_lsdyna(sim_file,self.lsdynapath)
        return
    def generateFibers(
        self,
        alpha_endocardium: float=-60,
        alpha_eepicardium: float=60,
        beta_endocardium: float=25,
        beta_epicardium: float=-65,
        ):
        dyna_writer = FiberGenerationDynaWriter(self.model)
        return

    def generatePurkinje(
        self,
        pointstx: float,
        pointsty: float,
        pointstz= float,
        inodeid: int,
        iedgeid: int,
        edgelen: float =2,
        ngen: float =50,
        nbrinit: int =8,
        nsplit: int =2, 
        ):
        dyna_writer = PurkinjeGenerationDynaWriter(self.model)
        return

    def run_lsdyna(sim_file: str, lsdynapath: str, NCPU: int=1, memory: int, options: str = ""):
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
        commands = ["mpirun", "-np", str(NCPU), lsdynapath, "i=" + sim_file, "memory=", str(memory),options]

        p = subprocess.run(commands, stdout=subprocess.PIPE)
        return
