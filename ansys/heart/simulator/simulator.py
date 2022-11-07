"""Simulator module."""
import os
import shutil

from ansys.heart.preprocessor.models import HeartModel
import ansys.heart.writer.dynawriter as writers


class Simulator:
    """
    Perform pre-simulation steps.

    Some extra info here

    """

    def __init__(self, model: HeartModel, lsdynapath: str = "") -> None:
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
        fibers: bool = 0,
        purkinje: bool = 0,
        zeropressure: bool = 0,
        ep: bool = 0,
        mechanics: bool = 0,
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
