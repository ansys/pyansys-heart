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
        workdir: str,
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
        workdir = os.path.join(self.model.info.path_to_model,"simulation")
        simulationdynawriter = BaseDynaWriter(self.model)
        if fibers:
            path_fibers = os.path.join(self.model.info.path_to_model,"fiber_generation")
            self.write_fibers(path_fibers)
            # run dyna
        if zeropressure:
            path_zeropressure = os.path.join(self.model.info.path_to_model,"zeropressure_generation")
            self.write_zeropressureconfiguration(path_zeropressure)
            if fibers:
                shutil.copy2(
                    os.path.join(path_fibers, "element_solid_ortho.k"),
                    os.path.join(path_zeropressure, "solid_elements.k"),
                )

            # run dyna
            nodes_stressfree = self.get_stressfreenodes(path_zeropressure)
            shutil.copy(os.path.join(path_zeropressure, "nodes.k"),os.path.join(path_zeropressure, "nodes_endofdiastole.k"))
            shutil.copy(os.path.join(path_zeropressure, nodes_stressfree),os.path.join(path_zeropressure, "nodes.k"))

        if purkinje:
            path_purkinje = os.path.join(self.model.info.path_to_model,"purkinje_generation")
            self.write_purkinje(path_purkinje)
            # if exist left ventricle: 
            # run dyna mainLeft
            shutil.copy2(os.path.join(path_purkinje, "purkinjenetwork.k"),os.path.join(path_purkinje, "purkinjenetworkLEFT.k"))
            simulationdynawriter.model.add_part("Left Purkinje")
            # if exist right ventricle: 
            # run dyna mainRight
            shutil.copy2(os.path.join(path_purkinje, "purkinjenetwork.k"),os.path.join(path_purkinje, "purkinjenetworkRIGHT.k"))
            simulationdynawriter.model.add_part("Right Purkinje")
            if zeropressure:
                shutil.copy2(os.path.join(path_zeropressure, "nodes.k"),os.path.join(path_purkinje, "nodes.k"))
            
        # if (ep) and not (mechanics):
        simulationdynawriter.include_files()
        # if not (ep) and (mechanics):
            # write mechanics
        # if ep and mechanics:
            # TODO add coupling stuff and ignore default Ca2+ mechanics active stress/replaced by EP simulation


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

    def write_purkinje(
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



    def get_stressfreenodes(workdir: str):
        """
        Find the result file after zerop
        Returns
        -------
        """
        # TODO check if converged
        guess_files = []
        for file in os.listdir(workdir):
            if file[-5:] == "guess":
                guess_files.append(file)

        return guess_files[-1]

