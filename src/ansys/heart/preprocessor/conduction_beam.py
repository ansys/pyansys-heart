# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Module containing class for creating conduxtion system."""

import networkx as nx
import numpy as np
import pyvista as pv

from ansys.heart.core import LOG as LOGGER
from ansys.heart.preprocessor.mesh.objects import BeamMesh, Point, SurfaceMesh
from ansys.heart.preprocessor.models import FourChamber


def _create_line(point_start: np.array, point_end: np.array, beam_length: float):
    """Create points in a line defined by a start point and an end point.

    Parameters
    ----------
    point_start : np.array
        Start Point.
    point_end : np.array
        End point.
    beam_length : float
        Beam length.

    Returns
    -------
    points:
        List of created points.
    """
    line_vector = point_end - point_start
    line_length = np.linalg.norm(line_vector)
    n_points = int(np.round(line_length / beam_length)) + 1
    points = np.zeros([n_points, 3])
    # beams = np.zeros([n_points - 1, 2])
    points = np.linspace(point_start, point_end, n_points)
    # beams[:, 0] = np.linspace(0, n_points - 2, n_points - 1, dtype=int)
    # beams[:, 1] = np.linspace(0, n_points - 2, n_points - 1, dtype=int) + 1
    return points


def _refine_line(nodes: np.array, beam_length: float):
    new_nodes = [nodes[0, :]]
    for beam_id in range(len(nodes) - 1):
        point_start = nodes[beam_id, :]
        point_end = nodes[beam_id + 1, :]
        points = _create_line(point_start, point_end, beam_length=beam_length)
        new_nodes = np.vstack((new_nodes, points[1:, :]))
    return new_nodes


class ConductionSystem:
    """Methods to generate conduction system."""

    def __init__(self, model: FourChamber):
        self.m = model

    def compute_SA_node(self, target_coord=None) -> Point:
        """
        Compute SinoAtrial node.

        SinoAtrial node is defined on the endocardium of the right atrium and
        between sup vena cava and inf vena cave.
        """
        if target_coord == None:
            for cap in self.m.right_atrium.caps:
                if "superior" in cap.name:
                    sup_vcava_centroid = cap.centroid
                elif "inferior" in cap.name:
                    inf_vcava_centroid = cap.centroid

            # define SinoAtrial node:
            target_coord = sup_vcava_centroid - (inf_vcava_centroid - sup_vcava_centroid) / 2

        right_atrium_endo = self.m.mesh.get_surface(self.m.right_atrium.endocardium.id)

        target_id = pv.PolyData(
            self.m.mesh.nodes[right_atrium_endo.global_node_ids, :]
        ).find_closest_point(target_coord)

        SA_node_id = right_atrium_endo.global_node_ids[target_id]

        SA_point = Point(name="SA_node", xyz=self.m.mesh.nodes[SA_node_id, :], node_id=SA_node_id)
        # TODO
        self.m.right_atrium.points.append(SA_point)

        return SA_point

    def compute_AV_node(self, target_coord=None) -> Point:
        """
        Compute Atrio-Ventricular node.

        AtrioVentricular node is on right artrium endocardium surface and closest to septum.

        Returns
        -------
        Point
            returns the AV node.
        """
        right_atrium_endo = self.m.mesh.get_surface(self.m.right_atrium.endocardium.id)

        if target_coord is None:
            for surface in self.m.right_ventricle.surfaces:
                if "endocardium" in surface.name and "septum" in surface.name:
                    right_septum = self.m.mesh.get_surface(surface.id)
            # define AtrioVentricular as the closest point to septum
            target_id = pv.PolyData(
                self.m.mesh.points[right_atrium_endo.global_node_ids, :]
            ).find_closest_point(right_septum.center)

        else:
            target_id = pv.PolyData(
                self.m.mesh.points[right_atrium_endo.global_node_ids, :]
            ).find_closest_point(target_coord)

        # assign a point
        av_id = right_atrium_endo.global_node_ids[target_id]
        AV_point = Point(name="AV_node", xyz=self.m.mesh.nodes[av_id, :], node_id=av_id)

        self.m.right_atrium.points.append(AV_point)

        return AV_point

    def compute_av_conduction(self, beam_length: float = 1.5) -> BeamMesh:
        """Compute Atrio-Ventricular conduction by means of beams following a geodesic path."""
        right_atrium_endo = self.m.mesh.get_surface(self.m.right_atrium.endocardium.id)

        try:
            SA_id = self.m.right_atrium.get_point("SA_node").node_id
        except AttributeError:
            LOGGER.info("SA node is not defined, creating with default option.")
            SA_id = self.m.compute_SA_node().node_id

        try:
            AV_id = self.m.right_atrium.get_point("AV_node").node_id
        except AttributeError:
            LOGGER.info("AV node is not defined, creating with default option.")
            AV_id = self.m.compute_AV_node().node_id

        #! get local SA/AV ids.
        SA_id_local = np.argwhere(right_atrium_endo.global_node_ids == SA_id).flatten()[0]
        AV_id_local = np.argwhere(right_atrium_endo.global_node_ids == AV_id).flatten()[0]
        path_SAN_AVN = right_atrium_endo.geodesic(SA_id_local, AV_id_local)
        beam_nodes = path_SAN_AVN.points

        beam_nodes = _refine_line(beam_nodes, beam_length=beam_length)[1:, :]

        # duplicate nodes inside the line, connect only SA node (the first) with 3D
        point_ids = np.linspace(0, len(beam_nodes) - 1, len(beam_nodes), dtype=int)
        point_ids = np.insert(point_ids, 0, SA_id)
        # build connectivity table
        edges = np.vstack((point_ids[:-1], point_ids[1:])).T

        mask = np.ones(edges.shape, dtype=bool)
        mask[0, 0] = False  # SA point at solid

        beam_net = self.m.add_beam_net(beam_nodes, edges, mask, pid=0, name="SAN_to_AVN")

        return beam_net

    def _get_hisbundle_bifurcation(self):
        """
        Define start points of the bundle of His.

        End point: create a point inside of septum part and close to AV node .
        """
        AV_node = self.m.right_atrium.get_point("AV_node")

        septum_point_ids = np.unique(np.ravel(self.m.mesh.tetrahedrons[self.m.septum.element_ids]))

        # remove nodes on surface, to make sure His bundle nodes are inside of septum
        septum_point_ids = np.setdiff1d(
            septum_point_ids, self.m.left_ventricle.endocardium.node_ids
        )
        septum_point_ids = np.setdiff1d(septum_point_ids, self.m.right_ventricle.septum.node_ids)

        septum_pointcloud = pv.PolyData(self.m.mesh.nodes[septum_point_ids, :])

        # Define start point: closest to artria
        pointcloud_id = septum_pointcloud.find_closest_point(AV_node.xyz)

        pointcloud_id = septum_pointcloud.find_closest_point(AV_node.xyz)

        bifurcation_id = septum_point_ids[pointcloud_id]
        bifurcation_coord = self.m.mesh.nodes[bifurcation_id, :]

        return bifurcation_coord

    def compute_His_conduction(self, beam_length: float = 1.5):
        """Compute His bundle conduction."""
        bifurcation_coord = self._get_hisbundle_bifurcation()

        # path start from AV point, to septum start point, then to septum end point
        av_id = None
        for beam in self.m.beam_network:
            if beam.name == "SAN_to_AVN":
                av_id = beam.edges[-1, -1]
                break

        if av_id is None:
            LOGGER.error(
                "Unable to find the last node of SAN_to_AVN branch, you need to define it."
            )
            exit()

        # create nodes from start to end
        # new_nodes = np.array(
        #     [
        #         self.m.right_atrium.get_point("AV_node").xyz,
        #         bifurcation_coord,
        #     ]
        # )
        sgmt_top, nodes = self.find_path(
            self.m.mesh,
            self.m.right_atrium.get_point("AV_node").xyz,
            bifurcation_coord,
        )
        new_nodes = self.m.mesh.points[nodes]
        new_nodes = _refine_line(new_nodes, beam_length=beam_length)
        new_nodes = new_nodes[1:, :]

        point_ids = np.concatenate(
            (
                np.array([av_id]),
                +np.linspace(0, len(new_nodes) - 1, len(new_nodes), dtype=int),
            )
        )
        bifurcation_id = point_ids[-1]
        # create beam
        edges = np.vstack((point_ids[:-1], point_ids[1:])).T

        (
            position_id_His_end_left,
            his_end_left_coord,
            new_nodes,
            edges,
            sgmt_left,
        ) = self._create_His_side(
            side="left",
            new_nodes=new_nodes,
            edges=edges,
            beam_length=beam_length,
            bifurcation_coord=bifurcation_coord,
            bifurcation_id=bifurcation_id,
        )
        (
            position_id_His_end_right,
            his_end_right_coord,
            new_nodes,
            edges,
            sgmt_right,
        ) = self._create_His_side(
            side="right",
            new_nodes=new_nodes,
            edges=edges,
            beam_length=beam_length,
            bifurcation_coord=bifurcation_coord,
            bifurcation_id=bifurcation_id,
        )
        # finally
        mask = np.ones(edges.shape, dtype=bool)
        mask[0, 0] = False  # AV point at previous, not offset in creation

        beam_net = self.m.add_beam_net(new_nodes, edges, mask, pid=0, name="His")
        beam_net.beam_nodes_mask[0, 0] = True  # offset in writer

        #

        surf = SurfaceMesh(
            name="his_bundle_segment",
            triangles=np.vstack((sgmt_top, sgmt_left, sgmt_right)),
            nodes=self.m.mesh.points,
        )
        #!? why?
        self.m.mesh.boundaries.append(surf)
        
        #! add surface to central mesh object for future use.
        surface_id = int(np.max(self.m.mesh.surface_ids)+1)
        self.m.mesh.add_surface(surf.clean(), surface_id)
        self.m.mesh._surface_id_to_name[surface_id] = "his_bundle_segment"
        self.m.mesh = self.m.mesh.clean()

        return Point(
            xyz=his_end_left_coord,
            node_id=beam_net.edges[position_id_His_end_left[0], position_id_His_end_left[1]],
        ), Point(
            xyz=his_end_right_coord,
            node_id=beam_net.edges[position_id_His_end_right[0], position_id_His_end_right[1]],
        )

    @staticmethod
    def find_path(
        mesh: pv.UnstructuredGrid, start: np.ndarray, end: np.ndarray, return_segment=True
    ):
        """Find shortest path between two nodes.

        Notes
        -----
        Unlike geodesic, this method searches a path inside of a 3D mesh.

        Parameters
        ----------
        mesh : pv.UnstructuredGrid
            Must be with tetra cells.
        start : np.ndarray
            Start point coordinates.
        end : np.ndarray
            End point coordinates
        return_segment : bool, optional
            Return a segment set (list of triangles) on which the path relies, by default True
        """
        #! mesh can now have multiple element types: TETRA, TRIANGLE, etc.
        mesh = mesh.extract_cells_by_type(pv.CellType.TETRA)

        def _mesh_to_nx_graph(mesh):
            # convert tetra mesh to graph
            G = nx.Graph()
            # Add nodes
            for i, point in enumerate(mesh.points):
                G.add_node(i, pos=tuple(point))
            # Assume all cells are tetra
            cells = np.array(mesh.cells).reshape(-1, 5)[:, 1:]
            # Add edges
            for cell in cells:
                G.add_edge(cell[0], cell[1])
                G.add_edge(cell[1], cell[2])
                G.add_edge(cell[2], cell[0])
                G.add_edge(cell[0], cell[3])
                G.add_edge(cell[1], cell[3])
                G.add_edge(cell[2], cell[3])
            return G

        # do the search in a small region for efficiency
        center = 0.5 * (start + end)
        radius = 3 * np.linalg.norm(start - center)
        sphere = pv.Sphere(center=center, radius=radius)

        # extract region
        cell_center = mesh.cell_centers()
        ids = np.where(cell_center.select_enclosed_points(sphere)["SelectedPoints"])[0]
        mesh.point_data["temp_ids"] = np.linspace(0, mesh.n_points - 1, mesh.n_points, dtype=int)
        sub_mesh = mesh.extract_cells(ids)

        # search shortes path across cells
        source_id = sub_mesh.find_closest_point(start)
        target_id = sub_mesh.find_closest_point(end)
        graph = _mesh_to_nx_graph(sub_mesh)
        # ids in submesh
        path = nx.shortest_path(graph, source=source_id, target=target_id)
        # ids in mesh
        path2 = sub_mesh["temp_ids"][path]

        if return_segment:
            tetras = sub_mesh.cells.reshape(-1, 5)[:, 1:]
            triangles = np.vstack(
                (
                    tetras[:, [0, 1, 2]],
                    tetras[:, [0, 1, 3]],
                    tetras[:, [0, 2, 3]],
                    tetras[:, [1, 2, 3]],
                )
            )
            segment = []
            for i, j in zip(path[0:-1], path[1:]):
                for tri in triangles:
                    if i in tri and j in tri:
                        segment.append(tri)
                        break
            segment = np.array(segment)
            segment2 = sub_mesh["temp_ids"][segment]

            return segment2, path2
        else:
            return path2

    def _create_His_side(
        self, side: str, new_nodes, edges, beam_length, bifurcation_coord, bifurcation_id
    ):
        """Create His side after bifucation."""
        if side.lower() == "left":
            endo = self.m.left_ventricle.endocardium
        elif side.lower() == "right":
            endo = self.m.right_ventricle.septum

        n = 20  # avoid too close to bifurcation point
        temp_id = pv.PolyData(self.m.mesh.points[endo.node_ids, :]).find_closest_point(
            bifurcation_coord, n=n
        )[n - 1]

        his_end_id = endo.node_ids[temp_id]
        his_end_coord = self.m.mesh.points[his_end_id, :]

        # side_his = np.array([bifurcation_coord, his_end_coord])
        sgmt, nodes = self.find_path(self.m.mesh, bifurcation_coord, his_end_coord)
        side_his = self.m.mesh.points[nodes]

        side_his = _refine_line(side_his, beam_length=beam_length)
        new_nodes = np.vstack((new_nodes, side_his[1:, :]))

        side_his_point_ids = np.concatenate(
            (
                np.array([bifurcation_id]),
                edges[-1, -1]
                + 1
                + np.linspace(0, len(side_his[1:, :]) - 1, len(side_his[1:, :]), dtype=int),
            )
        )

        edges = np.vstack(
            (edges, np.column_stack((side_his_point_ids[:-1], side_his_point_ids[1:])))
        )
        position_id_His_end = np.argwhere(edges == side_his_point_ids[-1])[0]
        return (position_id_His_end, his_end_coord, new_nodes, edges, sgmt)

    def compute_left_right_bundle(self, start_coord, start_id, side: str, beam_length: float = 1.5):
        """Bundle brunch."""
        if side == "Left":
            ventricle = self.m.left_ventricle
            endo_surface = self.m.mesh.get_surface(self.m.left_ventricle.endocardium.id)
        elif side == "Right":
            ventricle = self.m.right_ventricle
            surface_ids = [ventricle.endocardium.id, ventricle.septum.id]
            endo_surface = self.m.mesh.get_surface(surface_ids)

        #! this will give local ids.
        bundle_branch = endo_surface.geodesic(
            endo_surface.find_closest_point(start_coord),
            endo_surface.find_closest_point(self.m.mesh.points[ventricle.apex_points[0].node_id]),
        )

        new_nodes = bundle_branch.points
        new_nodes = _refine_line(new_nodes, beam_length=beam_length)
        # exclude first and last (apex) node which belongs to purkinje beam
        new_nodes = new_nodes[1:-1, :]
        point_ids = np.linspace(0, len(new_nodes) - 1, len(new_nodes), dtype=int)
        point_ids = np.insert(point_ids, 0, start_id)
        apex = ventricle.apex_points[0].node_id
        for network in self.m.beam_network:
            if network.name == side + "-purkinje":
                apex = network.edges[0, 0]
        point_ids = np.append(point_ids, apex)

        edges = np.vstack((point_ids[:-1], point_ids[1:])).T

        mask = np.ones(edges.shape, dtype=bool)
        mask[0, 0] = False  # His end point of previous, no offset at creation
        mask[-1, -1] = False  # Apex point, no offset

        beam_net = self.m.add_beam_net(new_nodes, edges, mask, pid=0, name=side + " bundle branch")
        # used in dynawriter, to write beam connectivity, need offset since it's a beam node
        beam_net.beam_nodes_mask[0, 0] = True
        beam_net.beam_nodes_mask[-1, -1] = True

        return beam_net

    @staticmethod
    def _get_closest_point_id(surface: pv.PolyData, point):
        # Surface contains all mesh node, find_closest_point could be wrong
        cell_id = surface.find_closest_cell(point)
        return surface.get_cell(cell_id).point_ids[0]

    def compute_Bachman_bundle(self, start_coord, end_coord, beam_length: float = 1.5):
        """Compute Bachman bundle conduction system."""
        la_epi = self.m.left_atrium.epicardium
        ra_epi = self.m.right_atrium.epicardium

        start_id = self._get_closest_point_id(ra_epi, start_coord)
        end_id = self._get_closest_point_id(la_epi, end_coord)

        epi = la_epi.merge(ra_epi)
        path = epi.geodesic(start_id, end_id)

        #
        beam_nodes = _refine_line(path.points, beam_length=beam_length)[1:-1, :]
        point_ids = np.linspace(0, len(beam_nodes) - 1, len(beam_nodes), dtype=int)
        point_ids = np.insert(point_ids, 0, start_id)
        point_ids = np.append(point_ids, end_id)

        # build connectivity table
        edges = np.vstack((point_ids[:-1], point_ids[1:])).T

        mask = np.ones(edges.shape, dtype=bool)
        mask[0, 0] = False  # Start point at solid
        mask[-1, 1] = False  # End point at solid

        tri = np.vstack((la_epi.triangles, ra_epi.triangles))
        surface = SurfaceMesh(name="Bachman segment", triangles=tri, nodes=self.m.mesh.nodes)
        self.m.mesh.boundaries.append(surface)
        
        #! add surface to central mesh. mesh.boundaries is deprecated.
        surface_id = int(np.max(self.m.mesh.surface_ids)+1)
        self.m.mesh.add_surface(surface.clean(), surface_id)
        self.m.mesh._surface_id_to_name[surface_id] = "Bachman segment"
        self.m.mesh.clean()
        
        beam_net = self.m.add_beam_net(beam_nodes, edges, mask, pid=0, name="Bachman bundle")

        return beam_net


if __name__ == "__main__":
    model: FourChamber = FourChamber.load_model(
        r"D:\ansysdev\pyansys-heart\downloads\Strocchi2020\01\FourChamber\heart_model.pickle"
    )

    test = ConductionSystem(model)
    sa_point = test.compute_SA_node()
    test.compute_AV_node()

    test.compute_av_conduction(
        # midpoints=[[-74, 90, 388], [70, 111, 372]]
    )

    # a,b =model.compute_His_conduction()
    a, b = test.compute_His_conduction()
    print(a.xyz, a.node_id)
    print(b.xyz, b.node_id)

    test.compute_left_right_bundle(a.xyz, a.node_id, side="Left")
    test.compute_left_right_bundle(b.xyz, b.node_id, side="Right")
    test.compute_Bachman_bundle(start_coord=sa_point.xyz, end_coord=np.array([-34, 163, 413]))

    model.plot_purkinje()
