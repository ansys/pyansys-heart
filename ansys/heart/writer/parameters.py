# Define the time unit
time_unit = "s"

# Scale parameters by given unit system
if time_unit == "s":
    scale = 1.0
elif time_unit == "ms":
    scale = 1000.0
else:
    print("Error in unit")

# define simulation parameters
parameters = {
    "Unit": {
        "Time": "s" if time_unit == "s" else "ms",
        "Consistent System": "kPa, mm, mN, s, kg" if time_unit == "s" else "MPa, mm, N, ms, g",
    },
    "Time": {
        "End Time": 5.0 * scale,
        "dtmin": 0.001 * scale,
        "dtmax": 0.01 * scale,
        "dt_d3plot": 0.05 * scale,
        "dt_icvout": 0.001 * scale,
    },
    "Material": {
        "Myocardium": {
            "Isotropic": {"rho": 1e-6 * scale, "k1": 2.36 / scale, "k2": 1.75},
            "Anisotropic": {"k1": 0.49 / scale, "k2": 9.01},
            "Active": {"Tmax": 125 / scale, "ca2ionm": 4.35, "Prefill": 1.0},
        },
        "Atrium": {
            "rho": 1e-6 * scale,
            "itype": -1,
            "mu1": 34.9 / scale,
            "alpha1": 2,
            "Comment": "Shoule be equivalent with MAT_077_H",
        },
    },
    "Cap": {"Density": 1.0e-6 * scale, "nu": 0.499, "c10": 1000 / scale, "Thickness": 5.0},
    "Global damping": 500.0 / scale,
    "Boundary Condition": {
        "Valve Spring": {"BV": 5.0 / scale, "4C": 20.0 / scale},
        "Normal Scale factor": 0.5,
        "Radial Scale factor": 1.0,
    },
    "Pericardium": {
        "Penalty function": [0.1, 25],
        "Spring Stiffness": 50.0 / scale,
        "Spring Type": "apex-mitral-direction",
    },
    "ED pressure": {"Left Ventricle": 2.0 / scale, "Right Ventricle": 0.53333 / scale},
    "Circulation System": {
        "Name": "ConstantPreloadWindkesselAfterload",
        "Left Ventricle": {
            "Constant": {
                "Rv": 5.0e-6,
                "Ra": 1.0e-5,
                "Rp": 1.2e-4,
                "Ca": 2.5e4 * scale,
                "Pven": None,  # define this later
            },
            "Initial Value": {"part_init": 8 / scale},
        },
        "Right Ventricle": {
            "Constant": {
                "Rv": 2.5e-6,
                "Ra": 0.35e-5,
                "Rp": 0.15e-4,
                "Ca": 11.25e4 * scale,
                "Pven": None,  # define this later
            },
            "Initial Value": {"part_init": 2 / scale},
        },
    },
}

parameters["Circulation System"]["Left Ventricle"]["Constant"]["Pven"] = parameters["ED pressure"][
    "Left Ventricle"
]
parameters["Circulation System"]["Right Ventricle"]["Constant"]["Pven"] = parameters["ED pressure"][
    "Right Ventricle"
]

if __name__ == "__main__":
    print(parameters)
