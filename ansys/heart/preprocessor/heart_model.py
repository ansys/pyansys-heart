import os

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
        """Loads all relevant model information for the json file

        Parameters
        ----------
        filename : str
            Path to model information in json format
        """

        # exposed to user
        logger.info("Loading heart model from: %s " % filename)
        model_info = self._mesh.load_mesh(filename)
        self.info = model_info
        return

    def dump_model(self, filename: str = "model_info.json", clean_working_directory: bool = False):
        """Dumps model information for future use. Exports simulation mesh in .vtk format

        Parameters
        ----------
        filename : str, optional
            Path to model info json, by default "model_info.json"
        clean_working_directory : bool, optional
            Flag indicating whether to clean the working directory of any
            temporary files, by default False
        """

        # exposed to user
        if clean_working_directory:
            import glob as glob

            files = glob.glob(os.path.join(self.info.working_directory, "*"))

            if len(files) > 0:
                logger.debug("Files detected: cleaning all files from directory")
            for f in files:
                os.remove(f)

        if not filename:
            filename = os.path.join(self.info.working_directory, "model_info.json")

        self.get_model_characteristics()

        # export the mesh
        self._mesh.export_mesh(filename)

        return

    def extract_simulation_mesh(self, remesh: bool = True):
        """Extracts the simulation mesh based on the model information provided

        Parameters
        ----------
        do_remesh : bool, optional
            Flag indicating whether to perform a remeshing step with the
            specified cell size, by default True
        """
        # exposed to user

        # perform actions on mesh object
        self._mesh.load_raw_mesh()
        self._mesh.extract_parts()
        self._mesh.add_cavities()

        if remesh:
            mesh_size = self.info.mesh_size
            logger.info("Remeshing uniformly with mesh size: %f " % mesh_size)
            self._mesh.remesh_volume(
                mesher="fluent", mesh_size=mesh_size
            )  # need to make this dynamic.
            self.info.is_remeshed = True

            # map data of original mesh to remeshed volume
            # includes (discrete) mapping of the tags
            self._mesh.map_data_to_remeshed_volume()
        else:
            myocardium_volume_vtk = self._mesh._extract_myocardium()
            self._mesh._vtk_volume = myocardium_volume_vtk
            self.info.is_remeshed = False
            self.info._mesh_size = None
            # raise ValueError("Using original input data for simulation " "not yet supported")
            logger.warning("Extracting mesh without remeshing not thoroughly tested yet")

        # extract surface from remeshed volume mesh
        self._mesh.get_surface_from_volume()

        # closes the cavities
        self._mesh.close_cavities()

        # extract endo/epicardium surfaces
        self._mesh.extract_endocardium_epicardium()

        # extract apex
        self._mesh.extract_apical_points()

        # validate cavities
        self._mesh._validate_cavities()

        # extract volumetric region of septum
        if self.info.model_type in ["BiVentricle", "FourChamber"]:
            self._mesh._extract_septum()

        # create element sets for myocardium
        self._mesh._create_myocardium_element_sets()

        return

    def get_model_characteristics(self, write_to_file: bool = True) -> dict:
        """Creates a dictionary of model characteristics

        Parameters
        ----------
        write_to_file : bool, optional
            Flag indicating whether to write to file , by default True

        Returns
        -------
        dict
            Dictionary of model characteristics
        """
        characteristics = {"model": {}, "mesh": {"cavity": []}}

        characteristics["model"]["type"] = self.info.model_type
        characteristics["mesh"]["num-cavities"] = len(self._mesh._cavities)
        num_caps = 0
        cavity_names = []
        for cavity in self._mesh._cavities:
            cavity_names.append(cavity.name)
            cap_names = []
            for cap in cavity.closing_caps:
                cap_names.append(cap.name)
                num_caps += 1
            cavity.volume
            cavity_info = {
                "name": cavity.name,
                "cap-names": cap_names,
                "volume": cavity.volume
                # "caps"      : cavity.closing_caps
            }
            characteristics["mesh"]["cavity"].append(cavity_info)

        characteristics["mesh"]["num-caps"] = num_caps
        characteristics["mesh"]["cavity-names"] = cavity_names

        # write to json
        import json

        json_path = os.path.join(self.info.working_directory, "model_characteristics.json")
        with open(json_path, "w") as outfile:
            json.dump(characteristics, indent=4, fp=outfile)

        return characteristics


if __name__ == "__main__":

    print("Protected")
