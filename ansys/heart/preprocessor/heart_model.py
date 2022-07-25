import os

from ansys.heart.preprocessor.heart_mesh import HeartMesh
from ansys.heart.preprocessor.model_information import ModelInformation

# import logger
from ansys.heart.custom_logging import LOGGER
from ansys.heart.preprocessor.vtk_module import (
    read_vtk_polydata_file,
    vtk_remove_arrays,
    write_vtkdata_to_vtkfile,
)


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
        LOGGER.info("Loading heart model from: %s " % filename)
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
                LOGGER.debug("Files detected: cleaning all files from directory")
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
        self._mesh.get_cavity_cap_intersections()

        if remesh:
            mesh_size = self.info.mesh_size
            LOGGER.info("Remeshing uniformly with mesh size: %f " % mesh_size)
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
            LOGGER.warning("Extracting mesh without remeshing not thoroughly tested yet")

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

    def extract_simulation_mesh_from_simplified_geometry(self):
        """Extracts a simulation mesh from a modified/simplified geometry
        from Strocchi or Cristobal et al
        """
        from ansys.heart.preprocessor.vtk_module import (
            read_vtk_polydata_file,
        )
        import numpy as np
        import copy

        # node-tag mapping:
        cavity_tag_map = {"Left ventricle": 1, "Right ventricle": 2}
        node_tag_map = {
            "Left ventricle": {
                "epicardium": 0,
                "mitral-valve-edge": 1,
                "aortic-valve-edge": 2,
                "endocardium": 3,
            },
            "Right ventricle": {
                "epicardium": 4,
                "pulmonary-valve-edge": 5,
                "tricuspid-valve-edge": 6,
                "interventricular-edge": 7,
                "endocardium": 8,
            },
        }

        # read surface mesh
        self._mesh._vtk_volume_raw = read_vtk_polydata_file(self.info.path_original_mesh)

        self._mesh._vtk_volume_temp = self._mesh._vtk_volume_raw

        # convert to PolyData
        self._mesh._vtk_surface = (
            self._mesh._vtk_volume_raw
        )  # = convert_to_polydata(self._mesh._vtk_volume_raw)

        # add cavities
        self._mesh.add_cavities()

        self._mesh.get_cavity_cap_intersections_simplified(node_tag_map)

        self._mesh.mesh_volume_from_simplified(node_tag_map, self.info.mesh_size)

        # get node sets from segment sets
        for cavity in self._mesh._cavities:
            # remove any defined node-sets
            LOGGER.debug("Removing %d nodesets" % len(cavity.node_sets))
            cavity.node_sets = []
            # NOTE: this duplicates the nodes on the intersection of the
            for segset in cavity.segment_sets:
                nodeset = copy.deepcopy(segset)
                nodeset["set"] = np.unique(segset["set"])
                cavity.node_sets.append(nodeset)

        self._mesh._validate_node_sets()
        self._mesh._validate_segment_sets()

        # extract the surface from the remeshed volume mesh
        self._mesh.get_surface_from_volume()

        # closes the cavities
        self._mesh.close_cavities()

        # extract apical points
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
