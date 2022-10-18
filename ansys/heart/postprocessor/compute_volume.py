"""Module for computing the volume."""
import json

import meshio
import numpy as np


class ClosedSurface:
    """Closed surface."""

    def __init__(self, coord, surface_mesh):
        """
        Init a closed surface.

        Parameters
        ----------
        coord: coordinate can have more nodes than surface mesh.
        surface_mesh: connectivity table of triangles.
        """
        ids, a = np.unique(surface_mesh, return_inverse=True)
        self.coord = coord[ids]
        self.connect = a.reshape(surface_mesh.shape)

        self.x = coord[surface_mesh][:, :, 0]
        self.y = coord[surface_mesh][:, :, 1]
        self.z = coord[surface_mesh][:, :, 2]

    def get_volume(self):
        """Compute volume.

        Note
        ----
        From numpy-stl
        http://www.geometrictools.com/Documentation/PolyhedralMassProperties.pdf.

        """

        def subexpression(x):
            w0, w1, w2 = x[:, 0], x[:, 1], x[:, 2]
            temp0 = w0 + w1
            f1 = temp0 + w2
            temp1 = w0 * w0
            temp2 = temp1 + w1 * temp0
            f2 = temp2 + w2 * f1
            f3 = w0 * temp1 + w1 * temp2 + w2 * f2
            g0 = f2 + w0 * (f1 + w0)
            g1 = f2 + w1 * (f1 + w1)
            g2 = f2 + w2 * (f1 + w2)
            return f1, f2, f3, g0, g1, g2

        x0, x1, x2 = self.x[:, 0], self.x[:, 1], self.x[:, 2]
        y0, y1, y2 = self.y[:, 0], self.y[:, 1], self.y[:, 2]
        z0, z1, z2 = self.z[:, 0], self.z[:, 1], self.z[:, 2]
        a1, b1, c1 = x1 - x0, y1 - y0, z1 - z0
        a2, b2, c2 = x2 - x0, y2 - y0, z2 - z0
        d0, d1, d2 = b1 * c2 - b2 * c1, a2 * c1 - a1 * c2, a1 * b2 - a2 * b1

        f1x, f2x, f3x, g0x, g1x, g2x = subexpression(self.x)
        f1y, f2y, f3y, g0y, g1y, g2y = subexpression(self.y)
        f1z, f2z, f3z, g0z, g1z, g2z = subexpression(self.z)

        intg = np.zeros(10)
        intg[0] = sum(d0 * f1x)
        intg[1:4] = sum(d0 * f2x), sum(d1 * f2y), sum(d2 * f2z)
        intg[4:7] = sum(d0 * f3x), sum(d1 * f3y), sum(d2 * f3z)
        intg[7] = sum(d0 * (y0 * g0x + y1 * g1x + y2 * g2x))
        intg[8] = sum(d1 * (z0 * g0y + z1 * g1y + z2 * g2y))
        intg[9] = sum(d2 * (x0 * g0z + x1 * g1z + x2 * g2z))
        intg /= np.array([6, 24, 24, 24, 60, 60, 60, 120, 120, 120])
        volume = intg[0]
        cog = intg[1:4] / volume
        cogsq = cog**2
        inertia = np.zeros((3, 3))
        inertia[0, 0] = intg[5] + intg[6] - volume * (cogsq[1] + cogsq[2])
        inertia[1, 1] = intg[4] + intg[6] - volume * (cogsq[2] + cogsq[0])
        inertia[2, 2] = intg[4] + intg[5] - volume * (cogsq[0] + cogsq[1])
        inertia[0, 1] = inertia[1, 0] = -(intg[7] - volume * cog[0] * cog[1])
        inertia[1, 2] = inertia[2, 1] = -(intg[8] - volume * cog[1] * cog[2])
        inertia[0, 2] = inertia[2, 0] = -(intg[9] - volume * cog[2] * cog[0])
        # return volume, cog, inertia
        return abs(volume)

    def write_stl(self, name):
        """
        Write closed surface in a stl file.

        Parameters
        ----------
        name: file name
        -------

        """
        meshio.write_points_cells(
            name + ".stl",
            self.coord,
            [("triangle", self.connect)],
        )


def update_system_json(nodes_file):
    """
    Update the unstressed cavity volume in the system model settings.

    Note
    ----
    How to use: run this script under the main simulation directory

    """
    # nodes_file = "iter4.guess"
    # BV or 4C
    lv_cavity_file = r"cavity_left_ventricle.segment"
    rv_cavity_file = r"cavity_right_ventricle.segment"

    #
    data = []
    with open(nodes_file) as f:
        for line in f.readlines():
            if line[0] != "*" and line[0] != "$":
                data.append(line)
    x_m = np.genfromtxt(data, delimiter=[8, 16, 16, 16])
    x_m = x_m[x_m[:, 0].argsort()][:, 1:]

    lv_cavity = np.loadtxt(lv_cavity_file, delimiter=",", dtype=int)
    rv_cavity = np.loadtxt(rv_cavity_file, delimiter=",", dtype=int)
    lv_volume = ClosedSurface.get_volume(x_m, lv_cavity)
    rv_volume = ClosedSurface(x_m, rv_cavity)
    print(lv_volume)
    print(rv_volume)

    json_file = r"system_model_settings.json"
    with open(json_file) as f:
        sys_dct = json.load(f)
    with open(json_file + "_old", "w") as f:
        f.write(json.dumps(sys_dct, indent=2, separators=(",", ": ")))

    sys_dct["SystemModelInitialValues"]["UnstressedVolumes"]["lv"] = lv_volume
    sys_dct["SystemModelInitialValues"]["UnstressedVolumes"]["rv"] = rv_volume

    with open(json_file, "w") as f:
        f.write(json.dumps(sys_dct, indent=2, separators=(",", ": ")))

    return


if __name__ == "__main__":
    update_system_json("nodes.k")
