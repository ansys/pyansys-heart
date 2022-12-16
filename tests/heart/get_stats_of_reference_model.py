"""Use this file to extract stats of a reference model.

Example
-------
>>> python -m get_stats_of_reference_model heart_model.pickle
"""
import os
import sys

import ansys.heart.preprocessor.models as models

if __name__ == "__main__":
    """
    Extracts stats of a reference model in json format.
    """

    try:
        path_to_reference_model = sys.argv[1]
    except:
        # default.
        path_to_reference_model = (
            "D:\\development\\pyheart-lib\\pyheart-lib\\tests\\heart"
            "\\assets\\reference_models\\strocchi2020\\01\\BiVentricle\\heart_model.pickle"
        )
    if not os.path.isfile(path_to_reference_model):
        print("Model not found. Please check path.")
        exit()

    print("Generating reference stats of model {}".format(path_to_reference_model))
    # Load model.
    model = models.HeartModel.load_model(path_to_reference_model)
    if isinstance(model, models.BiVentricle):
        model_type = "BiVentricle"
    elif isinstance(model, models.FullHeart):
        model_type = "FullHeart"
    else:
        print("Model type not supported.")
        exit()

    path_to_stats = "stats_reference_model_{:}.json".format(model_type)
    print("Storing info in: {}".format(path_to_stats))

    # store relevant stats in a dictionary.
    ref_stats = {}
    ref_stats["cavity_volumes"] = {
        model.left_ventricle.cavity.name: model.left_ventricle.cavity.volume,
        model.right_ventricle.cavity.name: model.right_ventricle.cavity.volume,
    }

    ref_stats["parts"] = dict.fromkeys([p.name for p in model.parts], {})
    for part in model.parts:
        ref_stats["parts"][part.name] = {"ntetra": 0}
        ref_stats["parts"][part.name]["ntetra"] = part.element_ids.shape[0]
        ref_stats["parts"][part.name]["surfaces"] = {}
        for surface in part.surfaces:
            ref_stats["parts"][part.name]["surfaces"][surface.name] = {
                "nfaces": surface.faces.shape[0]
            }

    ref_stats["mesh"] = {"boundaries": {}}
    for b in model.mesh.boundaries:
        ref_stats["mesh"]["boundaries"][b.name] = {"nfaces": b.faces.shape[0]}

    import json

    with open(path_to_stats, "w") as outfile:
        outfile.write(json.dumps(ref_stats, indent=4))
