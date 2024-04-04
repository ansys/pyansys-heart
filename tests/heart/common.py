"""Some common functions to test parts."""
import os

global heart_version

heart_version = os.getenv("ANSYS_HEART_MODEL_VERSION")
if not heart_version:
    heart_version = "v0.1"

if heart_version == "v0.2":
    import ansys.heart.preprocessor.models.v0_2.models as models

elif heart_version == "v0.1":
    import ansys.heart.preprocessor.models.v0_1.models as models

import numpy as np


def compare_part_names(model: models.HeartModel, ref_stats: dict):
    """Test if parts match that of the reference model.

    Note
    ----
    1. tests consistency of part names
    """

    for part_name in model.part_names:
        assert part_name in list(ref_stats["parts"].keys()), (
            "Part name: %s does not exist in reference model" % part_name
        )

    pass


def compare_part_element_ids(model: models.HeartModel, reference_model: models.HeartModel):
    """Test if parts match that of the reference model.

    Note
    ----
    1. tests element ids defined in all parts
    """
    # assert isinstance(reference_model, models.HeartModel), "Expecting model of type HeartModel"

    for part in model.parts:
        ref_part = next(
            ref_part for ref_part in reference_model.parts if part.name == ref_part.name
        )
        assert np.array_equal(part.element_ids, ref_part.element_ids), (
            "%s element-ids do not match with reference model" % part.name
        )

    pass


def compare_surface_names(model: models.HeartModel, ref_stats: dict):
    """Test if surfaces of the parts match with the reference model.

    Note
    ----
    1. tests consistency of surface names
    """

    for part in model.parts:
        # check if surface is in the list of reference surface names
        try:
            ref_surface_names = list(ref_stats["parts"][part.name]["surfaces"].keys())
        except KeyError:
            continue

        for surface_name in ref_surface_names:
            assert surface_name in part.surface_names, (
                "%s does not exist in model but exists in reference model" % surface_name
            )

    pass


def compare_generated_mesh(model: models.HeartModel, ref_stats: dict):
    """Compares the number of tetrahedrons, faces, etc in the model.

    Note
    ----
    Test conditions:
        Difference num tetrahedrons < 100
        Difference num faces < 10
        Difference num faces of valve boundaries < 5
    """
    import os

    if os.name == "nt":
        allowed_difference1 = 0
        allowed_difference2 = 0
        allowed_difference3 = 0
    else:
        allowed_difference1 = 500
        allowed_difference2 = 75
        allowed_difference3 = 10

    # Compare parts.
    for part in model.parts:
        try:
            ref_num_tetra = ref_stats["parts"][part.name]["ntetra"]
        except KeyError:
            continue
        difference = abs(part.element_ids.shape[0] - ref_num_tetra)
        assert difference <= allowed_difference1, (
            "{0}: Difference between reference and generated model is {1} exceeds {2}"
        ).format(part.name, difference, allowed_difference1)

        for surface in part.surfaces:
            # surf_name = "-".join(surface.name.lower().split())
            try:
                ref_nfaces = ref_stats["parts"][part.name]["surfaces"][surface.name]["nfaces"]
            except KeyError:
                print(surface.name + "not found")
                continue
            difference = abs(surface.n_faces - ref_nfaces)
            assert difference <= allowed_difference2, (
                "Boundary: {0} Difference between reference and generated model"
                " is {1} exceeds {2}"
            ).format(part.name, difference, allowed_difference2)

    # Compare other boundaries in Mesh
    for boundary in model.mesh.boundaries:
        ref_nfaces = ref_stats["mesh"]["boundaries"][boundary.name]["nfaces"]

        if "valve" in boundary.name:
            allowed_difference = allowed_difference3
        else:
            allowed_difference = allowed_difference2
        difference = abs(boundary.n_faces - ref_nfaces)
        assert difference <= allowed_difference, (
            "Boundary: {0} Difference between reference and generated model is {1} exceeds {2}"
        ).format(part.name, difference, allowed_difference2)


def compare_cavity_volume(model: models.HeartModel, ref_volumes: dict):
    """Test volume of cavities.

    Note
    ----
    1. Volume of cavity
    """
    # assert isinstance(model, models.HeartModel), "Expecting model of type HeartModel"
    for part in model.parts:
        if not part.cavity:
            continue

        ref_volume = ref_volumes["cavity_volumes"][part.name]
        assert abs(part.cavity.surface.volume - ref_volume) < 1e-2 * ref_volume, (
            "Difference in cavity volume of model %s exceeds 1 percent" % part.name
        )

    pass


def compare_part_names_2(model: models.HeartModel, ref_stats: dict):
    """Test if parts match that of the reference model.

    Notes
    -----
    1. tests consistency of part names
    """
    for part_name in model.part_names:
        assert part_name in list(ref_stats["PARTS"].keys()), (
            "Part name: %s does not exist in reference model" % part_name
        )

    pass


def compare_part_element_ids_2(model: models.HeartModel, reference_model: models.HeartModel):
    """Test if parts match that of the reference model.

    Notes
    -----
    1. tests element ids defined in all parts
    """
    # assert isinstance(reference_model, models.HeartModel), "Expecting model of type HeartModel"

    for part in model.parts:
        ref_part = next(
            ref_part for ref_part in reference_model.parts if part.name == ref_part.name
        )
        assert np.array_equal(part.element_ids, ref_part.element_ids), (
            "%s element-ids do not match with reference model" % part.name
        )

    pass


def compare_surface_names_2(model: models.HeartModel, ref_stats: dict):
    """Test if surfaces of the parts match with the reference model.

    Notes
    -----
    1. tests consistency of surface names
    """
    for part in model.parts:
        # check if surface is in the list of reference surface names
        try:
            ref_surface_names = list(ref_stats["PARTS"][part.name]["SURFACES"].keys())
        except KeyError:
            continue

        for surface_name in ref_surface_names:
            assert surface_name in part.surface_names, (
                "%s does not exist in model but exists in reference model" % surface_name
            )

    pass


def compare_generated_mesh_2(model: models.HeartModel, ref_stats: dict):
    """Compares the number of tetrahedrons, faces, etc in the model.

    Notes
    -----
    Test conditions:
        Difference num tetrahedrons < 100
        Difference num faces < 10
        Difference num faces of valve boundaries < 5
    """
    import os

    if os.name == "nt":
        allowed_difference1 = 0
        allowed_difference2 = 0
        allowed_difference3 = 0
    else:
        allowed_difference1 = 500
        allowed_difference2 = 75
        allowed_difference3 = 10

    # Compare parts.
    for part in model.parts:
        try:
            ref_num_tetra = ref_stats["PARTS"][part.name]["num_tets"]
        except KeyError:
            continue
        difference = abs(part.element_ids.shape[0] - ref_num_tetra)
        assert difference <= allowed_difference1, (
            "{0}: Difference between reference and generated model is {1} exceeds {2}"
        ).format(part.name, difference, allowed_difference1)

        for surface in part.surfaces:
            # surf_name = "-".join(surface.name.lower().split())
            try:
                ref_nfaces = ref_stats["PARTS"][part.name]["SURFACES"][surface.name]["num_faces"]
            except KeyError:
                print(surface.name + "not found")
                continue
            difference = abs(surface.n_faces - ref_nfaces)
            assert difference <= allowed_difference2, (
                "Boundary: {0} Difference between reference and generated model"
                " is {1} exceeds {2}"
            ).format(part.name, difference, allowed_difference2)

    # # Compare other boundaries in Mesh
    # for boundary in model.mesh.boundaries:
    #     ref_nfaces = ref_stats["mesh"]["boundaries"][boundary.name]["nfaces"]

    #     if "valve" in boundary.name:
    #         allowed_difference = allowed_difference3
    #     else:
    #         allowed_difference = allowed_difference2
    #     difference = abs(boundary.n_faces - ref_nfaces)
    #     assert difference <= allowed_difference, (
    #         "Boundary: {0} Difference between reference and generated model is {1} exceeds {2}"
    #     ).format(part.name, difference, allowed_difference2)


def compare_cavity_volume_2(model: models.HeartModel, ref_volumes: dict):
    """Test volume of cavities.

    Notes
    -----
    1. Volume of cavity
    """
    # assert isinstance(model, models.HeartModel), "Expecting model of type HeartModel"
    for part in model.parts:
        if not part.cavity:
            continue

        ref_volume = ref_volumes["CAVITIES"][part.name]["volume"]
        assert abs(part.cavity.surface.volume - ref_volume) < 1e-2 * ref_volume, (
            "Difference in cavity volume of model %s exceeds 1 percent" % part.name
        )

    pass
