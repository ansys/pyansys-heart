"""Module containing classes for the various heart models which
can be created with the preprocessor
"""
import json
import pathlib
import typing

from ansys.heart.preprocessor.model_definitions import HEART_PARTS, MODELS, TAG_MAP

# NOTE:
# Parts can contain features. Such as surfaces:

# """
# Models can have attributes such as mesh and parts
# """
class Mesh:
    """Mesh class: contains nodal coordinates and element definitions"""

    def __init__(self) -> None:
        self.elements = None
        """(Volume) elements of the mesh"""
        self.nodes = None
        """Nodes of the mesh"""
        self.cell_data = None
        """Data per mesh cell/element"""
        self.point_data = None
        """Data per mesh point"""
        pass

    def read_mesh_file(filename: pathlib.Path, file_format: str):
        """Reads mesh file"""
        return


class Feature:
    """Feature class"""

    def __init__(self, name: str = None) -> None:
        self.name = name
        """Name of feature"""
        self.type = None
        """Type of feature"""
        pass


class Cavity(Feature):
    """Cavity class"""

    def __init__(self) -> None:
        super().__init__()
        self.type = "cavity"


class Edges(Feature):
    """Edges class"""

    def __init__(self) -> None:
        super().__init__()
        self.type = "edges"


class Surface(Feature):
    """Surface class"""

    def __init__(self, name: str = None) -> None:
        super().__init__(name)

        self.faces = None
        """Faces of surface"""
        self.nodes = None
        """Global node ids of surface"""
        self.type = "surface"


class Cap(Feature):
    """Cap class"""

    def __init__(self, name: str = None) -> None:
        super().__init__(name)
        self.triangles = None
        """Triangulation of cap"""
        self.normal = None
        """Normal of cap"""
        self.centroid = None
        """Centroid of cap"""
        self.type = "cap"


class Part:
    """Part class"""

    @property
    def _features(self) -> typing.List[Feature]:
        """Returns list of part features"""
        features = []
        for key, value in self.__dict__.items():
            attribute = getattr(self, key)
            if type(attribute) == Feature:
                features.append(attribute)
        return features

    @property
    def surfaces(self) -> typing.List[Surface]:
        surfaces = []
        for key, value in self.__dict__.items():
            if type(value) == Surface:
                surfaces.append(value)
        return surfaces

    @property
    def surface_names(self) -> typing.List[str]:
        surface_names = []
        for key, value in self.__dict__.items():
            if type(value) == Surface:
                surface_names.append(value.name)
        return surface_names

    def __init__(self, name: str = None, part_type: str = None) -> None:
        self.name = name
        """Name of the part"""
        self.id = None
        """Part ID"""
        self.part_type = part_type
        """Type of the part"""
        self.tag_labels = None
        """VTK tag labels used in this part"""
        self.tag_ids = None
        """VTK tag ids used in this part"""
        self._add_surfaces()

    def _add_surfaces(self):
        """Adds surfaces to the part"""

        if self.part_type in ["ventricle", "atrium"]:
            self.endocardium = Surface("{0} endocardium".format(self.name))
            """Endocardium"""
            self.epicardium = Surface("{0} epicardium".format(self.name))
            """Epicardium"""
            if self.part_type == "ventricle":
                self.septum = Surface("{0} septum".format(self.name))
                """Septum surface"""
        elif self.part_type in ["artery"]:
            self.wall = Surface("{0} wall".format(self.name))
            """Wall"""
        return

    def _add_myocardium_part(self):
        self.myocardium = Part(name="myocardium", part_type="myocardium")
        return

    def _add_septum_part(self):
        self.septum = Part(name="septum", part_type="septum")
        return


class ModelInfo:
    """Contains model information"""

    def __init__(
        self, database: str, work_directory: pathlib.Path, path_to_case: pathlib.Path
    ) -> None:
        self.database = database
        self.workdir = work_directory
        self.path_to_original_mesh = path_to_case
        pass

    @property
    def database(self) -> str:
        return self._database

    @database.setter
    def database(self, value: str):
        valid_databases = ["Strocchi2020", "Cristobal2021"]
        if value not in valid_databases:
            raise ValueError(
                "{0} not a valid database name. Please specify one of the following database names: {1}".format(
                    value, valid_databases
                )
            )
        self._database = value


class HeartModel:
    """Parent class for heart models"""

    @property
    def parts(self) -> typing.List[Part]:
        """Returns list of parts"""
        parts = []
        for key, value in self.__dict__.items():
            attribute = getattr(self, key)
            if type(attribute) == Part:
                parts.append(attribute)
        return parts

    @property
    def part_names(self) -> typing.List[str]:
        """Returns list of part names"""
        part_names = []
        for part in self.parts:
            part_names.append(part.name)
        return part_names

    def __init__(self, info: ModelInfo) -> None:
        self.info = info
        """Model meta information"""
        self.mesh = Mesh()
        """Model mesh"""

        self._add_subparts()
        """Adds any subparts"""
        self._add_labels_to_parts()
        """Adds appropiate vtk labels to the parts"""
        pass

    def _add_part(self, part_name: str):
        """Dynamically adds a part as an attribute to the object"""
        setattr(self, "_".join(part_name.lower().split()), Part(name=part_name))
        return

    def _add_labels_to_parts(self):
        """Uses model definitions to add corresponding vtk labels to the part"""
        for part in self.parts:
            if part.name == "Septum":
                continue
            part.tag_labels = HEART_PARTS[part.name]["VTKLabels"]
            # add tag ids
            part.tag_ids = []
            for tag_label in part.tag_labels:
                part.tag_ids.append(TAG_MAP[self.info.database][tag_label])
        return

    def _add_subparts(self):
        """Adds subparts to parts of type ventricle"""
        for part in self.parts:
            if part.part_type in ["ventricle"]:
                part._add_myocardium_part()
                if "Left ventricle" in part.name:
                    part._add_septum_part()

        return

    def print_parts(self):
        """Prints the involved part and vtk label ids"""
        for part in self.parts:
            print("{0} : {1}".format(part.name, part.tag_ids))
        return

    def read_input_mesh(filename: pathlib.Path = None):
        """Reads the input mesh"""


class LeftVentricle(HeartModel):
    """Model of just the left ventricle"""

    def __init__(self, info: ModelInfo) -> None:

        self.left_ventricle = Part(name="Left ventricle", part_type="ventricle")
        """Left ventricle part"""
        # remove septum - not used in left ventricle only model
        del self.left_ventricle.septum

        super().__init__(info)
        pass


class BiVentricle(HeartModel):
    """Model of the left and right ventricle"""

    def __init__(self, info: ModelInfo) -> None:

        self.left_ventricle = Part(name="Left ventricle", part_type="ventricle")
        """Left ventricle part"""
        self.right_ventricle = Part(name="Right ventricle", part_type="ventricle")
        """Right ventricle part"""
        self.septum = Part(name="Septum", part_type="septum")
        """Septum"""

        super().__init__(info)
        pass


class FourChamber(HeartModel):
    """Model of the left/right ventricle and left/right atrium"""

    def __init__(self, info: ModelInfo) -> None:

        self.left_ventricle = Part(name="Left ventricle", part_type="ventricle")
        """Left ventricle part"""
        self.right_ventricle = Part(name="Right ventricle", part_type="ventricle")
        """Right ventricle part"""
        self.septum = Part(name="Septum", part_type="septum")
        """Septum"""

        self.left_atrium = Part(name="Left atrium", part_type="atrium")
        """Left atrium part"""
        self.right_atrium = Part(name="Right atrium", part_type="atrium")
        """Right atrium part"""

        super().__init__(info)

        pass


class FullHeart(HeartModel):
    """Model of the left/right ventricle,  left/right atrium, aorta
    and pulmonary artery
    """

    def __init__(self, info: ModelInfo) -> None:

        self.left_ventricle = Part(name="Left ventricle", part_type="ventricle")
        """Left ventricle part"""
        self.right_ventricle = Part(name="Right ventricle", part_type="ventricle")
        """Right ventricle part"""
        self.septum = Part(name="Septum", part_type="septum")
        """Septum"""
        self.left_atrium = Part(name="Left atrium", part_type="atrium")
        """Left atrium part"""
        self.right_atrium = Part(name="Right atrium", part_type="atrium")
        """Right atrium part"""
        self.aorta = Part(name="Aorta", part_type="artery")
        """Aorta part"""
        self.pulmonary_artery = Part(name="Pulmonary artery", part_type="artery")
        """Pulmonary artery part"""

        super().__init__(info)

        pass


if __name__ == "__main__":
    info = ModelInfo(database="Strocchi2020", work_directory="tmp", path_to_case="test.case")

    model = LeftVentricle(info)
    print("LeftVentricle:")
    model.print_parts()

    model = BiVentricle(info)
    print("BiVentricle:")
    model.print_parts()

    model = FourChamber(info)
    print("FourChamber:")
    model.print_parts()

    model = FullHeart(info)
    print("FullHeart:")
    model.print_parts()

    print(model.part_names)
