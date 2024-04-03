# populate a summary dictionary:
import os 
from ansys.heart.custom_logging import LOGGER
heart_version = os.environ["ANSYS_HEART_MODEL_VERSION"]

if not heart_version or heart_version=="v0.1":
    from ansys.heart.preprocessor.models.v0_1.models import HeartModel        
elif heart_version=="v0.2":
    from ansys.heart.preprocessor.models.v0_2.models import HeartModel        

def model_summary(model: HeartModel) -> dict:
    """Generates a dictionary with model information

    Parameters
    ----------
    model : HeartModel
        HeartModel for which to generate the summary dictionary

    Returns
    -------
    dict
        Dictionary with organized model information.
    """

    sum_dict = {}
    sum_dict["GENERAL"] = {}
    
    try:
        sum_dict["GENERAL"]["total_num_tets"] = model.mesh.tetrahedrons.shape[0]
        sum_dict["GENERAL"]["total_num_nodes"] = model.mesh.nodes.shape[0]
    except TypeError:
        LOGGER.info("Failed to format General model information.")

    
    sum_dict["PARTS"] = {}
    sum_dict["CAVITIES"] = {}
    for ii, part in enumerate(model.parts):
        sum_dict["PARTS"][part.name] = {}       
        sum_dict["PARTS"][part.name]["num_tets"] = len(part.element_ids)
        
        sum_dict["PARTS"][part.name]["SURFACES"] = {}
        sum_dict["PARTS"][part.name]["CAPS"] = {}
        
        for surface in part.surfaces:
            sum_dict["PARTS"][part.name]["SURFACES"][surface.name] = {}
            sum_dict["PARTS"][part.name]["SURFACES"][surface.name]["num_faces"] = surface.triangles.shape[0]
        for cap in part.caps:
            sum_dict["PARTS"][part.name]["CAPS"][cap.name] = {}
            sum_dict["PARTS"][part.name]["CAPS"][cap.name]["num_nodes"] = len(cap.node_ids)
            
    for cavity in model.cavities:
        sum_dict["CAVITIES"][cavity.name] = {}
        sum_dict["CAVITIES"][cavity.name]["volume"] = cavity.surface.volume

    return sum_dict
