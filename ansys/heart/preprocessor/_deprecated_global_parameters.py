"""Contains global parameters and constants
"""

# cavity definitions
CAVITY_DEFINITIONS = {
    "LeftVentricle": ["Left ventricle myocardium", "Aortic valve plane", "Mitral valve plane"],
    "RightVentricle": [
        "Right ventricle myocardium",
        "Tricuspid valve plane",
        "Pulmonary valve plane",
    ],
    "LeftAtrium": [
        "Left atrium myocardium",
        "Mitral valve plane",
        "Left atrium appendage inlet",
        "Left superior pulmonary vein inlet",
        "Left inferior pulmonary vein inlet",
        "Right inferior pulmonary vein inlet",
        "Right superior pulmonary vein inlet",
    ],
    "RightAtrium": [
        "Right atrium myocardium",
        "Tricuspid valve plane",
        "Superior vena cava inlet",
        "Inferior vena cava inlet",
    ],
}

# dictionary of valid models
VALID_MODELS = {
    "LeftVentricle": {
        "LabelsToUse": ["Left ventricle myocardium", "Aortic valve plane", "Mitral valve plane"],
        "CavityDefinition": [CAVITY_DEFINITIONS["LeftVentricle"]],
    },
    "BiVentricle": {
        "LabelsToUse": [
            "Left ventricle myocardium",
            "Aortic valve plane",
            "Mitral valve plane",
            "Right ventricle myocardium",
            "Tricuspid valve plane",
            "Pulmonary valve plane",
        ],
        "CavityDefinition": [
            CAVITY_DEFINITIONS["LeftVentricle"],
            CAVITY_DEFINITIONS["RightVentricle"],
        ],
    },
    "FourChamber": {
        "LabelsToUse": [
            "Left ventricle myocardium",
            "Aortic valve plane",
            "Mitral valve plane",
            "Right ventricle myocardium",
            "Tricuspid valve plane",
            "Pulmonary valve plane",
            "Left atrium myocardium",
            "Left atrium appendage inlet",
            "Left superior pulmonary vein inlet",
            "Left inferior pulmonary vein inlet",
            "Right inferior pulmonary vein inlet",
            "Right superior pulmonary vein inlet",
            "Right atrium myocardium",
            "Superior vena cava inlet",
            "Inferior vena cava inlet",
        ],
        "CavityDefinition": [
            CAVITY_DEFINITIONS["LeftVentricle"],
            CAVITY_DEFINITIONS["RightVentricle"],
            CAVITY_DEFINITIONS["LeftAtrium"],
            CAVITY_DEFINITIONS["RightAtrium"],
        ],
    },
    # IMPROVED MODELS:
    "LeftVentricleImproved": {
        "LabelsToUse": ["Left ventricle myocardium", "Aortic valve plane", "Mitral valve plane"],
        "CavityDefinition": [CAVITY_DEFINITIONS["LeftVentricle"]],
    },
    "BiVentricleImproved": {
        "LabelsToUse": [
            "Left ventricle myocardium",
            "Aortic valve plane",
            "Mitral valve plane",
            "Right ventricle myocardium",
            "Tricuspid valve plane",
            "Pulmonary valve plane",
        ],
        "CavityDefinition": [
            CAVITY_DEFINITIONS["LeftVentricle"],
            CAVITY_DEFINITIONS["RightVentricle"],
        ],
    },
    "FourChamberImproved": {
        "LabelsToUse": [
            "Left ventricle myocardium",
            "Aortic valve plane",
            "Mitral valve plane",
            "Right ventricle myocardium",
            "Tricuspid valve plane",
            "Pulmonary valve plane",
            "Left atrium myocardium",
            "Left atrium appendage inlet",
            "Left superior pulmonary vein inlet",
            "Left inferior pulmonary vein inlet",
            "Right inferior pulmonary vein inlet",
            "Right superior pulmonary vein inlet",
            "Right atrium myocardium",
            "Superior vena cava inlet",
            "Inferior vena cava inlet",
            "Left atrial appendage border",
            "Left superior pulmonary vein border",
            "Left inferior pulmonary vein border",
            "Right inferior pulmonary vein border",
            "Right superior pulmonary vein border",
            "Superior vena cava border",
            "Inferior vena cava border",
        ],
        "CavityDefinition": [
            CAVITY_DEFINITIONS["LeftVentricle"],
            CAVITY_DEFINITIONS["RightVentricle"],
            CAVITY_DEFINITIONS["LeftAtrium"],
            CAVITY_DEFINITIONS["RightAtrium"],
        ],
    },
    "FullHeartImproved": {
        "LabelsToUse": [
            "Left ventricle myocardium",
            "Aortic valve plane",
            "Mitral valve plane",
            "Right ventricle myocardium",
            "Tricuspid valve plane",
            "Pulmonary valve plane",
            "Left atrium myocardium",
            "Left atrium appendage inlet",
            "Left superior pulmonary vein inlet",
            "Left inferior pulmonary vein inlet",
            "Right inferior pulmonary vein inlet",
            "Right superior pulmonary vein inlet",
            "Right atrium myocardium",
            "Superior vena cava inlet",
            "Inferior vena cava inlet",
            "Left atrial appendage border",
            "Left superior pulmonary vein border",
            "Left inferior pulmonary vein border",
            "Right inferior pulmonary vein border",
            "Right superior pulmonary vein border",
            "Superior vena cava border",
            "Inferior vena cava border",
            "Aorta wall",
            "Pulmonary artery wall",
        ],
        "CavityDefinition": [
            CAVITY_DEFINITIONS["LeftVentricle"],
            CAVITY_DEFINITIONS["RightVentricle"],
            CAVITY_DEFINITIONS["LeftAtrium"],
            CAVITY_DEFINITIONS["RightAtrium"],
        ],
    },
}

# dictionary of improved valid models
CAVITY_DEFINITIONS_IMPROVED = {
    "LeftVentricle": ["Left ventricle myocardium", "Aortic valve plane", "Mitral valve plane"],
    "RightVentricle": [
        "Right ventricle myocardium",
        "Tricuspid valve plane",
        "Pulmonary valve plane",
    ],
    "LeftAtrium": [
        "Left atrium myocardium",
        "Mitral valve plane",
        "Left atrium appendage inlet",
        "Left superior pulmonary vein inlet",
        "Left inferior pulmonary vein inlet",
        "Right inferior pulmonary vein inlet",
        "Right superior pulmonary vein inlet",
    ],
    "RightAtrium": [
        "Right atrium myocardium",
        "Tricuspid valve plane",
        "Superior vena cava inlet",
        "Inferior vena cava inlet",
    ],
}

VALID_MODELS_IMPROVED = {}

if __name__ == "__main__":
    print("Protected")
