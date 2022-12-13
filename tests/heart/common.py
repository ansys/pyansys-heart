"""Some common functions to test parts."""
import ansys.heart.preprocessor.models as models
import numpy as np


def compare_part_names(model: models.HeartModel, reference_model: models.HeartModel):
    """Test if parts match that of the reference model.

    Note
    ----
    1. tests consistency part names
    """
    assert isinstance(reference_model, models.HeartModel), "Expecting model of type HeartModel"

    for part_name in model.part_names:
        assert part_name in reference_model.part_names, (
            "Part name: %s does not exist in reference model" % part_name
        )

    pass


def compare_part_element_ids(model: models.HeartModel, reference_model: models.HeartModel):
    """Test if parts match that of the reference model.

    Note
    ----
    1. tests element ids defined in all parts
    """
    assert isinstance(reference_model, models.HeartModel), "Expecting model of type HeartModel"

    for part in model.parts:
        ref_part = next(
            ref_part for ref_part in reference_model.parts if part.name == ref_part.name
        )
        assert np.array_equal(part.element_ids, ref_part.element_ids), (
            "%s element-ids do not match with reference model" % part.name
        )

    pass


def compare_surface_names(model: models.HeartModel, reference_model: models.HeartModel):
    """Test if surfaces of the parts match with the reference model.

    Note
    ----
    1. tests consistency of surface names
    """

    for part in model.parts:
        ref_part = next(
            ref_part for ref_part in reference_model.parts if part.name == ref_part.name
        )

        # check surface names in reference part
        for surface_name in ref_part.surface_names:
            assert surface_name in part.surface_names, (
                "%s does not exist in model but exists in reference model" % surface_name
            )

    pass


def compare_surface_faces(model: models.HeartModel, reference_model: models.HeartModel):
    """Test if surfaces of the parts match with the reference model.

    Note
    ----
    1. tests surface topology (face id and ordering) of surfaces
    """

    for part in model.parts:
        ref_part = next(
            ref_part for ref_part in reference_model.parts if part.name == ref_part.name
        )

        # check surface contents in reference part
        for ref_surface in ref_part.surfaces:
            surface = next(surface for surface in part.surfaces if surface.name == ref_surface.name)

            assert np.array_equal(surface.faces, ref_surface.faces), (
                "Surface topology of reference surface %s is different" % ref_surface.name
            )
    pass


def compare_caps_nodeids(model: models.HeartModel, reference_model: models.HeartModel):
    """Test caps and cap definitions."""
    assert isinstance(reference_model, models.HeartModel), "Expecting model of type HeartModel"
    for part in model.parts:
        ref_part = next(
            ref_part for ref_part in reference_model.parts if part.name == ref_part.name
        )

        for cap in part.caps:
            ref_cap = next(ref_cap for ref_cap in ref_part.caps if ref_cap.name == cap.name)
            assert np.array_equal(ref_cap.node_ids, cap.node_ids), (
                "One or more nodes in %s do not match" % ref_cap.name
            )

    pass


def compare_caps_num_nodeids(model: models.HeartModel, reference_model: models.HeartModel):
    """Test caps and cap definitions."""
    assert isinstance(reference_model, models.HeartModel), "Expecting model of type HeartModel"
    for part in model.parts:
        ref_part = next(
            ref_part for ref_part in reference_model.parts if part.name == ref_part.name
        )

        for cap in part.caps:
            ref_cap = next(ref_cap for ref_cap in ref_part.caps if ref_cap.name == cap.name)
            assert len(cap.node_ids > 50), "%s does not contain more than 50 nodes" % ref_cap.name

    pass


def compare_cavity_topology(model: models.HeartModel, reference_model: models.HeartModel):
    """Test topology of cavities.

    Note
    ----
    1. Topology of cavity surface
    """
    assert isinstance(reference_model, models.HeartModel), "Expecting model of type HeartModel"
    for part in model.parts:
        ref_part = next(
            ref_part for ref_part in reference_model.parts if part.name == ref_part.name
        )
        if not part.cavity:
            continue

        ref_num_faces = ref_part.cavity.surface.faces.shape[0]
        num_faces = part.cavity.surface.faces.shape[0]
        assert ref_num_faces == num_faces, (
            "Cavity of part {0} has {1} faces in the reference model, "
            "but model has {2} faces".format(ref_part.name, ref_num_faces, num_faces)
        )

        assert np.array_equal(part.cavity.surface.faces, ref_part.cavity.surface.faces), (
            "Cavity faces of part %s do not match" % ref_part.name
        )

    pass


def compare_cavity_volume(model: models.HeartModel, reference_model: models.HeartModel):
    """Test volume of cavities.

    Note
    ----
    1. Volume of cavity
    """
    assert isinstance(reference_model, models.HeartModel), "Expecting model of type HeartModel"
    for part in model.parts:
        ref_part = next(
            ref_part for ref_part in reference_model.parts if part.name == ref_part.name
        )
        if not part.cavity:
            continue

        assert abs(part.cavity.volume - ref_part.cavity.volume) < 1e-2 * ref_part.cavity.volume, (
            "Difference in cavity volume of model %s exceeds 1%" % ref_part.name
        )

    pass
