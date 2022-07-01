import os
from ansys.heart.preprocessor.fluenthdf5_module import fluenthdf5_to_vtk

from ansys.heart.preprocessor.heart_mesh import HeartMesh
from ansys.heart.preprocessor.mesh_module import mesh_by_fluentmeshing
from ansys.heart.preprocessor.model_information import ModelInformation

# import logger
from ansys.heart.custom_logging import logger
from ansys.heart.preprocessor.vtk_module import (
    get_tri_info_from_polydata,
    read_vtk_unstructuredgrid_file,
    threshold_vtk_data,
    threshold_vtk_data_integers,
    vtk_surface_to_stl,
    write_vtkdata_to_vtkfile,
)
from vtk.numpy_interface import dataset_adapter as dsa  # this is an improved numpy integration


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

    def extract_simulation_mesh_from_simplified_geometry(self):
        """Extracts a simulation mesh from a modified/simplified geometry"""
        from ansys.heart.preprocessor.vtk_module import (
            read_vtk_unstructuredgrid_file,
            vtk_surface_to_stl,
            create_vtk_surface_triangles,
        )
        from ansys.heart.preprocessor.mesh_module import add_solid_name_to_stl
        from ansys.heart.preprocessor.fluenthdf5_module import fluenthdf5_to_vtk
        from vtk.numpy_interface import (
            dataset_adapter as dsa,
        )  # this is an improved numpy integration
        import numpy as np
        import vtk

        # node-tag mapping:
        node_tag_map = {
            "Left ventricle": {
                "epicardium": 0,
                "mitral-valve-edge": 1,
                "aortic-valve-edge": 2,
                "endocardium": 3,
            },
            "Right ventricle": {
                "epicardium": 0,
                "pulmonary-valve-edge": 1,
                "tricuspid-valve-edge": 2,
                "interventricular-edge": 3,
                "endocardium": 4,
            },
        }

        # read surface mesh
        self._mesh._vtk_surface = read_vtk_unstructuredgrid_file(self.info.path_original_mesh)

        # add list of cavities
        self._mesh._cavities = self._mesh._add_list_of_cavities()

        # convert to PolyData
        geom = vtk.vtkGeometryFilter()
        geom.SetInputData(self._mesh._vtk_surface)
        geom.Update()
        self._mesh._vtk_surface = geom.GetOutput()

        filename = os.path.join(self.info.working_directory, "p05_1.stl")
        vtk_surface_to_stl(self._mesh._vtk_surface, filename)

        # write separate stl per part and per endo/epicardial surface
        points, tris, cell_data, point_data = get_tri_info_from_polydata(self._mesh._vtk_surface)

        for cavity in self._mesh._cavities:
            tag = cavity.vtk_ids[0]
            part_mask = cell_data["tags"] == tag

            for surface_name in ["endocardium", "epicardium"]:
                # only select faces where one of the nodes is part of the endocardium
                node_tag = node_tag_map[cavity.name][surface_name]

                name_of_stl = "_".join(cavity.name.lower().split() + [surface_name])
                stl_path = os.path.join(
                    self.info.working_directory,
                    "part_{0}.stl".format(name_of_stl),
                )

                point_ids = np.where(point_data["node-tags"] == node_tag)[0]
                surface_mask = np.any(np.isin(tris, point_ids), axis=1)
                tris_to_write = tris[np.all(np.vstack([part_mask, surface_mask]), axis=0)]

                # filt_surface = threshold_vtk_data(vtk_surface, tag, tag, data_name="tags")[0]
                vtk_surface = create_vtk_surface_triangles(points, tris_to_write)

                # write separate files for endo/epicardium

                # write one stl surface per part

                vtk_surface_to_stl(vtk_surface, stl_path)

                name_of_part = "-".join(cavity.name.lower().split() + [surface_name])
                add_solid_name_to_stl(stl_path, "{0}".format(name_of_part), "binary")

        # create volume mesh
        mesh_output = os.path.join(self.info.working_directory, "fluent_volume_mesh.msh.h5")
        mesh_by_fluentmeshing(
            stl_path, mesh_output, mesh_size=self.info.mesh_size, journal_type="simplified_geometry"
        )

        tets, face_zones, nodes = fluenthdf5_to_vtk(
            mesh_output, mesh_output.replace(".msh.h5", ".vtk")
        )

        # assign face zones to cavity segment sets
        for cavity in self._mesh._cavities:
            # append segment-sets:
            cavity.segment_sets.append()

            # append node-sets

        # find edge loops that make up the valve

        # find intersection between face zones

        # use vtk surface to tag volume elements

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
