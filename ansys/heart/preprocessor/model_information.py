import os
import json

from ansys.heart.preprocessor.global_parameters import (
    VALID_MODELS,
    CAVITY_DEFINITIONS,
)

# import logger
from ansys.heart.custom_logging import logger


class ModelInformation:
    """Class that stores Model's meta data such as:
    working directory, original source file, parameters used, etc.
    Also performs some checks whether the 
    provided arguments are valid and as expected.

    Note
    ------
    This class may be subject to change
    """

    def __init__(
        self,
        model_type: str = "BiVentricle",
        database_name: str = "Strocchi2020",
        path_original_mesh: str = "",
        working_directory: str = "",
        is_remeshed: bool = False,
        path_mesh: str = "",
        path_to_cap_file: str = "",
    ):

        # required
        self.model_type = model_type
        self.database_name = database_name
        self.path_original_mesh = path_original_mesh
        self.working_directory = working_directory

        # optional
        self.is_remeshed = is_remeshed
        self.mesh_size = None
        self.path_mesh = path_mesh
        self.path_to_cap_file = path_to_cap_file

        # path to pickle object that stores all other data
        self.path_to_pickle = ""

    # list of properties
    @property
    def model_type(self) -> str:
        return self._model_type

    @model_type.setter
    def model_type(self, value: str):
        valid_models = list(VALID_MODELS.keys())
        if value not in valid_models:
            raise ValueError(
                "{0} not a valid model. Please specify one of the following: {1}".format(
                    value, valid_models
                )
            )
        self._model_type = value
        self.vtk_labels_to_use = VALID_MODELS[self._model_type]["LabelsToUse"]
        logger.info("Model Type: %s" % value)
        logger.info("Labels to use:")
        for label in self.vtk_labels_to_use:
            logger.info("\t%s" % label)
        return

    @property
    def database_name(self) -> str:
        return self._database_name

    @database_name.setter
    def database_name(self, value: str):
        valid_databases = ["Strocchi2020", "Cristobal2021"]
        if value not in valid_databases:
            raise ValueError(
                "Database {0} not recoqnized: may not be able "
                "to map labels to tag ids. Valid database names include: {1}".format(
                    value, valid_databases
                )
            )
        self._database_name = value
        # update labels to extract based on the name of the database
        self._get_list_of_tags_to_use(self._database_name)
        return

    @property
    def vtk_labels_to_use(self) -> dict:
        """This property is inherently linked to the model_type and hence 
        set by @model_type.setter
        """
        return self._vtk_labels_to_use

    @vtk_labels_to_use.setter
    def vtk_labels_to_use(self, value: dict):
        self._vtk_labels_to_use = value

    @property
    def mesh_size(self) -> float:
        """This specifies the mesh size used for remeshing
        the volume
        """
        return self._mesh_size

    @mesh_size.setter
    def mesh_size(self, value: float):
        # could add some checks on the mesh size
        if not value:
            self._mesh_size = 2.0
        else:
            self._mesh_size = value

    def _get_list_of_tags_to_use(self, database_name: str):
        """Relates the labels to tag ids. This is specific to the 
        database the vtk belongs to. That is: from the 24 hearts of Strocchi2020
        or to the database of 20 healthy hearts of Cristobal2021"""

        file_dir = os.path.dirname(__file__)
        f = open(os.path.join(file_dir, "database_label_mapping.json"))
        label_mapping = json.load(f)
        label_mapping = label_mapping[database_name]

        # get list of ids to extract
        vtk_labels_to_use = {}
        for label in self.vtk_labels_to_use:
            vtk_labels_to_use[label] = label_mapping[label]

        # this overwrites the original labels - and is now a dictionary
        self.vtk_labels_to_use = vtk_labels_to_use

        return

    def save_to_file(self, filename: str):
        """Saves model information to file
        """
        logger.info("Writing model info to file: %s" % filename)
        with open(filename, "w") as file:
            file.write(json.dumps(self.__dict__, indent=4))

    def load_from_file(self, filename: str):
        """Loads model information from file
        """
        logger.info("Reading model info from file: %s" % filename)
        with open(filename, "r") as file:
            tmp = json.load(file)

        self.model_type = tmp["_model_type"]
        self.database_name = tmp["_database_name"]

        self.is_remeshed = tmp["is_remeshed"]
        self.mesh_size = tmp["_mesh_size"]
        self.path_original_mesh = tmp["path_original_mesh"]
        self.working_directory = tmp["working_directory"]

        self.path_mesh = tmp["path_mesh"]
        self.path_to_cap_file = tmp["path_to_cap_file"]
        self.path_to_pickle = tmp["path_to_pickle"]

    def clean_working_directory(self):
        """Cleans working directory from any file
        """
        import glob as glob

        files = glob.glob(os.path.join(self.working_directory, "*"))
        if len(files) > 0:
            logger.debug("Files detected: cleaning all files from directory")
        for f in files:
            os.remove(f)

        return
