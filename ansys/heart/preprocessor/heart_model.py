from distutils.command.clean import clean
import os
from re import S
from tkinter import SOLID

from multiprocessing.sharedctypes import Value

from typing import List

# NOTE: do specific import!
from ansys.heart.preprocessor.mesh_module import *
from ansys.heart.preprocessor.vtk_module import *

from ansys.heart.preprocessor.heart_mesh import HeartMesh
from ansys.heart.preprocessor.model_information import ModelInformation

# import logger
from ansys.heart.custom_logging import logger


class HeartBoundaryConditions:
    """Class containing relevant functions for the boundary conditions"""

    def __init__(self):
        """Initializes object"""
        return


class HeartModel:
    """Class which represents the heart model of a specific type"""

    # some properties:

    def __init__(self, model_info: ModelInformation):
        """Initializes the heart model by making use of the specified
        model information object

        Parameters
        ----------
        model_info : ModelInformation
            Class which specifies relevant information
        """

        self.info = model_info

        self._mesh = HeartMesh(model_info)

        return

    def load_model(self, filename: str):
        """Loads the model from a model info file
        """
        # exposed to user
        logger.info("Loading heart model from: %s " % filename)
        model_info = self._mesh.load_mesh(filename)
        self.info = model_info
        return

    def dump_model(
        self, filename: str = None, clean_working_directory: bool = False
    ):
        """Dumps model information for future use. Exports
        the simulation-ready mesh in vtk format and the list of
        cavities which are part of the HeartMesh object
        """
        # exposed to user
        if clean_working_directory:
            import glob as glob

            files = glob.glob(os.path.join(self.info.working_directory, "*"))

            if len(files) > 0:
                logger.debug(
                    "Files detected: cleaning all files from directory"
                )
            for f in files:
                os.remove(f)

        if not filename:
            filename = os.path.join(
                self.info.working_directory, "model_info.json"
            )

        # export the mesh
        self._mesh.export_mesh(filename)

        return

    def extract_simulation_mesh(self):
        """Extracts the simulation mesh based on 
        the model information provided
        """
        # exposed to user

        # perform actions on mesh object
        self._mesh.load_raw_mesh()
        self._mesh.extract_parts()
        self._mesh.add_cavities()

        do_remesh = self.info

        if do_remesh:
            mesh_size = self.info.mesh_size
            logger.info("Remeshing uniformly with mesh size: %f " % mesh_size)
            self._mesh.remesh_volume(
                mesher="fluent", mesh_size=mesh_size
            )  # need to make this dynamic.
            self.info.is_remeshed = True
        else:
            raise ValueError(
                "Using original input data for simulation " "not yet supported"
            )

        # add remeshed volume to mesh object:
        self._mesh.set_volume_mesh_vtk(self.info.path_mesh)

        # map data of original mesh to remeshed volume
        # includes (discrete) mapping of the tags
        self._mesh.map_data_to_remeshed_volume()

        # extract surface from remeshed volume mesh
        self._mesh.get_surface_from_volume()

        # closes the cavities
        self._mesh.close_cavities()

        # extract endo/epicardium surfaces
        self._mesh.extract_endocardium_epicardium()

        # validate cavities
        self._mesh._validate_cavities()

        # extract volumetric region of septum 
        if self.info.model_type in ["BiVentricle", "FourChamber"]:
            self._mesh._extract_septum()

        return


if __name__ == "__main__":

    print("Protected")
