import numpy as np
import meshio
import vtk
import json


def get_cavity_volume(cavity):
    ids, a = np.unique(cavity, return_inverse=True)
    coords = x_m[ids]
    connectivity = a.reshape(cavity.shape)
    meshio.write_points_cells(
        "cavity.stl",
        coords,
        [("triangle", connectivity)],
    )

    return compute_stl_volume("cavity.stl")


def compute_stl_volume(f_stl):
    reader = vtk.vtkSTLReader()  # noqa
    reader.SetFileName(f_stl)
    reader.Update()

    mass_property = vtk.vtkMassProperties()  # noqa
    mass_property.SetInputData(reader.GetOutput())
    mass_property.Update()
    volume = mass_property.GetVolume()  # mm3

    return float(volume)


if __name__ == "__main__":
    """
    Update the unstressed cavity volume in the system model settings.
    
    How to use:
    Copy and run this script under the main simulation directory
    """

    nodes_file = "nodes.k"
    # BV or 4C
    lv_cavity_file = r"cavity_left_ventricle.segment"
    rv_cavity_file = r"cavity_right_ventricle.segment"

    #
    data = []
    with open(nodes_file) as f:
        for line in f.readlines():
            if line[0] != "*" and line[0] != "$":
                data.append(line)
    x_m = np.genfromtxt(data, delimiter=[8, 16, 16, 16])[:, 1:4]

    lv_cavity = np.loadtxt(lv_cavity_file,delimiter=',',dtype=int)
    rv_cavity = np.loadtxt(rv_cavity_file,delimiter=',',dtype=int)
    lv_volume = get_cavity_volume(lv_cavity)
    rv_volume = get_cavity_volume(rv_cavity)

    json_file = r"system_model_settings.json"
    with open(json_file) as f:
        sys_dct = json.load(f)
    with open(json_file + "_old", "w") as f:
        f.write(json.dumps(sys_dct, indent=2, separators=(",", ": ")))

    sys_dct["SystemModelInitialValues"]["UnstressedVolumes"]["lv"] = lv_volume
    sys_dct["SystemModelInitialValues"]["UnstressedVolumes"]["rv"] = rv_volume

    with open(json_file, "w") as f:
        f.write(json.dumps(sys_dct, indent=2, separators=(",", ": ")))
