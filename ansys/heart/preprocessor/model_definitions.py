"""Model definitions."""

# Dictionary of heart parts
HEART_PARTS = {
    "Left ventricle": {
        "VTKLabels": [
            "Left ventricle myocardium",
            "Aortic valve plane",
            "Mitral valve plane",
        ]
    },
    "Right ventricle": {
        "VTKLabels": [
            "Right ventricle myocardium",
            "Tricuspid valve plane",
            "Pulmonary valve plane",
        ]
    },
    "Left atrium": {
        "VTKLabels": [
            "Left atrium myocardium",
            "Left atrium appendage inlet",
            "Left superior pulmonary vein inlet",
            "Left inferior pulmonary vein inlet",
            "Right inferior pulmonary vein inlet",
            "Right superior pulmonary vein inlet",
            "Left atrial appendage border",
            "Left superior pulmonary vein border",
            "Left inferior pulmonary vein border",
            "Right inferior pulmonary vein border",
            "Right superior pulmonary vein border",
        ]
    },
    "Right atrium": {
        "VTKLabels": [
            "Right atrium myocardium",
            "Superior vena cava inlet",
            "Inferior vena cava inlet",
            "Superior vena cava border",
            "Inferior vena cava border",
        ]
    },
    "Aorta": {"VTKLabels": ["Aorta wall"]},
    "Pulmonary artery": {"VTKLabels": ["Pulmonary artery wall"]},
}
# List of valid models
MODELS = {
    "LeftVentricle": {"Parts": ["Left ventricle"]},
    "BiVentricle": {"Parts": ["Left ventricle", "Right ventricle"]},
    "FourChamber": {"Parts": ["Left ventricle", "Right ventricle", "Left atrium", "Right atrium"]},
    "FullHeart": {
        "Parts": [
            "Left ventricle",
            "Right ventricle",
            "Left atrium",
            "Right atrium",
            "Aorta",
            "Pulmonary artery",
        ]
    },
}
# map for mapping labels to ID
LABELS_TO_ID = {
    "Strocchi2020": {
        "Left ventricle myocardium": 1,
        "Right ventricle myocardium": 2,
        "Left atrium myocardium": 3,
        "Right atrium myocardium": 4,
        "Aorta wall": 5,
        "Pulmonary artery wall": 6,
        "Left atrial appendage border": 7,
        "Left superior pulmonary vein border": 8,
        "Left inferior pulmonary vein border": 9,
        "Right inferior pulmonary vein border": 10,
        "Right superior pulmonary vein border": 11,
        "Superior vena cava border": 12,
        "Inferior vena cava border": 13,
        "Mitral valve plane": 14,
        "Tricuspid valve plane": 15,
        "Aortic valve plane": 16,
        "Pulmonary valve plane": 17,
        "Left atrium appendage inlet": 18,
        "Left superior pulmonary vein inlet": 19,
        "Left inferior pulmonary vein inlet": 20,
        "Right inferior pulmonary vein inlet": 21,
        "Right superior pulmonary vein inlet": 22,
        "Superior vena cava inlet": 23,
        "Inferior vena cava inlet": 24,
    },
    "Rodero2021": {
        "Left ventricle myocardium": 1,
        "Right ventricle myocardium": 2,
        "Left atrium myocardium": 3,
        "Right atrium myocardium": 4,
        "Aorta wall": 5,
        "Pulmonary artery wall": 6,
        "Left atrial appendage border": 18,
        "Left superior pulmonary vein border": 21,
        "Left inferior pulmonary vein border": 20,
        "Right inferior pulmonary vein border": 19,
        "Right superior pulmonary vein border": 22,
        "Superior vena cava border": 23,
        "Inferior vena cava border": 24,
        "Mitral valve plane": 7,
        "Tricuspid valve plane": 8,
        "Aortic valve plane": 9,
        "Pulmonary valve plane": 10,
        "Left atrium appendage inlet": 11,
        "Left superior pulmonary vein inlet": 12,
        "Left inferior pulmonary vein inlet": 13,
        "Right inferior pulmonary vein inlet": 14,
        "Right superior pulmonary vein inlet": 15,
        "Superior vena cava inlet": 16,
        "Inferior vena cava inlet": 17,
    },
    "LabeledSurface": {
        "Left ventricle myocardium": 100,
        "Right ventricle myocardium": 101,
        "Left ventricle endocardium": 1,
        "Left ventricle epicardium": 2,
        "Mitral valve plane": 3,
        "Aortic valve plane": 4,
        "Right ventricle endocardium": 5,
        "Right ventricle epicardium": 6,
        "Pulmonary valve plane": 7,
        "Tricuspid valve plane": 8,
    },
}
