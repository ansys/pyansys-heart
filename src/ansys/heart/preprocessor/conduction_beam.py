# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
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
from ansys.heart.core.models import FourChamber
from ansys.heart.core.objects import BeamMesh, CapType, Point, SurfaceMesh


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

    def compute_sa_node(self, target_coord=None) -> Point:
        """
        Compute SinoAtrial node.

        SinoAtrial node is defined on the endocardium of the right atrium and
        between sup vena cava and inf vena cave.
        """
        if target_coord is None:
            sup_vcava_centroid = next(
                cap.centroid
                for cap in self.m.right_atrium.caps
                if cap.type == CapType.SUPERIOR_VENA_CAVA
            )
            inf_vcava_centroid = next(
                cap.centroid
                for cap in self.m.right_atrium.caps
                if cap.type == CapType.INFERIOR_VENA_CAVA
            )

            # define SinoAtrial node:
            target_coord = sup_vcava_centroid - (inf_vcava_centroid - sup_vcava_centroid) / 2

        right_atrium_endo = self.m.mesh.get_surface(self.m.right_atrium.endocardium.id)

        target_id = pv.PolyData(
            self.m.mesh.points[right_atrium_endo.global_node_ids_triangles, :]
        ).find_closest_point(target_coord)

        sino_atrial_node_id = right_atrium_endo.global_node_ids_triangles[target_id]

        sino_atrial_point = Point(
            name="SA_node",
            xyz=self.m.mesh.points[sino_atrial_node_id, :],
            node_id=sino_atrial_node_id,
        )
        self.m.right_atrium.points.append(sino_atrial_point)

        return sino_atrial_point

    def compute_av_node(self, target_coord=None) -> Point:
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
                self.m.mesh.points[right_atrium_endo.global_node_ids_triangles, :]
            ).find_closest_point(right_septum.center)

        else:
            target_id = pv.PolyData(
                self.m.mesh.points[right_atrium_endo.global_node_ids_triangles, :]
            ).find_closest_point(target_coord)

        # assign a point
        av_id = right_atrium_endo.global_node_ids_triangles[target_id]
        atrioventricular_point = Point(
            name="AV_node", xyz=self.m.mesh.points[av_id, :], node_id=av_id
        )

        self.m.right_atrium.points.append(atrioventricular_point)

        return atrioventricular_point

    def compute_av_conduction(self, beam_length: float = 1.5) -> BeamMesh:
        """Compute Atrio-Ventricular conduction by means of beams following a geodesic path."""
        right_atrium_endo = self.m.mesh.get_surface(self.m.right_atrium.endocardium.id)

        try:
            sino_atrial_id = self.m.right_atrium.get_point("SA_node").node_id
        except AttributeError:
            LOGGER.info("SA node is not defined, creating with default option.")
            sino_atrial_id = self.m.compute_sa_node().node_id

        try:
            atrio_ventricular_id = self.m.right_atrium.get_point("AV_node").node_id
        except AttributeError:
            LOGGER.info("AV node is not defined, creating with default option.")
            atrio_ventricular_id = self.m.compute_av_node().node_id

        #! get local SA/AV ids.
        sino_atrial_id_local = np.argwhere(
            right_atrium_endo.global_node_ids_triangles == sino_atrial_id
        ).flatten()[0]
        atrio_ventricular_id_local = np.argwhere(
            right_atrium_endo.global_node_ids_triangles == atrio_ventricular_id
        ).flatten()[0]
        path_sinoatrial_atrioventricular = right_atrium_endo.geodesic(
            sino_atrial_id_local, atrio_ventricular_id_local
        )
        beam_nodes = path_sinoatrial_atrioventricular.points

        beam_nodes = _refine_line(beam_nodes, beam_length=beam_length)[1:, :]

        # duplicate nodes inside the line, connect only SA node (the first) with 3D
        point_ids = np.linspace(0, len(beam_nodes) - 1, len(beam_nodes), dtype=int)
        point_ids = np.insert(point_ids, 0, sino_atrial_id)
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
        atrio_ventricular_node = self.m.right_atrium.get_point("AV_node")

        septum_point_ids = np.unique(np.ravel(self.m.mesh.tetrahedrons[self.m.septum.element_ids]))

        # remove nodes on surface, to make sure His bundle nodes are inside of septum
        septum_point_ids = np.setdiff1d(
            septum_point_ids,
            self.m.mesh.get_surface(self.m.left_ventricle.endocardium.id).global_node_ids_triangles,
        )
        septum_point_ids = np.setdiff1d(
            septum_point_ids,
            self.m.mesh.get_surface(self.m.right_ventricle.septum.id).global_node_ids_triangles,
        )

        septum_pointcloud = pv.PolyData(self.m.mesh.points[septum_point_ids, :])

        # Define start point: closest to artria
        pointcloud_id = septum_pointcloud.find_closest_point(atrio_ventricular_node.xyz)

        pointcloud_id = septum_pointcloud.find_closest_point(atrio_ventricular_node.xyz)

        bifurcation_id = septum_point_ids[pointcloud_id]
        bifurcation_coord = self.m.mesh.points[bifurcation_id, :]

        return bifurcation_coord

    def compute_his_conduction(self, beam_length: float = 1.5) -> tuple[BeamMesh, Point, Point]:
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
            position_id_his_end_left,
            his_end_left_coord,
            new_nodes,
            edges,
            sgmt_left,
        ) = self._create_his_side(
            side="left",
            new_nodes=new_nodes,
            edges=edges,
            beam_length=beam_length,
            bifurcation_coord=bifurcation_coord,
            bifurcation_id=bifurcation_id,
        )
        (
            position_id_his_end_right,
            his_end_right_coord,
            new_nodes,
            edges,
            sgmt_right,
        ) = self._create_his_side(
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

        #! add surface to central mesh object for future use.
        surface_id = int(np.max(self.m.mesh.surface_ids) + 1)
        self.m.mesh.add_surface(surf.clean(), surface_id, name="his_bundle_segment")
        self.m.mesh = self.m.mesh.clean()

        return (
            beam_net,
            Point(
                xyz=his_end_left_coord,
                node_id=beam_net.edges[position_id_his_end_left[0], position_id_his_end_left[1]],
            ),
            Point(
                xyz=his_end_right_coord,
                node_id=beam_net.edges[position_id_his_end_right[0], position_id_his_end_right[1]],
            ),
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
            graph = nx.Graph()
            # Add nodes
            for i, point in enumerate(mesh.points):
                graph.add_node(i, pos=tuple(point))
            # Assume all cells are tetra
            cells = np.array(mesh.cells).reshape(-1, 5)[:, 1:]
            # Add edges
            for cell in cells:
                graph.add_edge(cell[0], cell[1])
                graph.add_edge(cell[1], cell[2])
                graph.add_edge(cell[2], cell[0])
                graph.add_edge(cell[0], cell[3])
                graph.add_edge(cell[1], cell[3])
                graph.add_edge(cell[2], cell[3])
            return graph

        # do the search in a small region for efficiency
        center = 0.5 * (start + end)
        radius = 3 * np.linalg.norm(start - center)
        sphere = pv.Sphere(center=center, radius=radius)

        # extract region
        cell_center = mesh.cell_centers()
        ids = np.where(cell_center.select_enclosed_points(sphere)["SelectedPoints"])[0]
        sub_mesh = mesh.extract_cells(ids)

        # search shortes path across cells
        source_id = sub_mesh.find_closest_point(start)
        target_id = sub_mesh.find_closest_point(end)
        graph = _mesh_to_nx_graph(sub_mesh)
        # ids in submesh
        path = nx.shortest_path(graph, source=source_id, target=target_id)
        # ids in mesh
        path2 = sub_mesh["_global-point-ids"][path]

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
            segment2 = sub_mesh["_global-point-ids"][segment]

            return segment2, path2
        else:
            return path2

    def _create_his_side(
        self, side: str, new_nodes, edges, beam_length, bifurcation_coord, bifurcation_id
    ):
        """Create His side after bifucation."""
        if side.lower() == "left":
            endo = self.m.mesh.get_surface(self.m.left_ventricle.endocardium.id)
        elif side.lower() == "right":
            endo = self.m.mesh.get_surface(self.m.right_ventricle.septum.id)

        n = 20  # avoid too close to bifurcation point
        temp_id = pv.PolyData(
            self.m.mesh.points[endo.global_node_ids_triangles, :]
        ).find_closest_point(bifurcation_coord, n=n)[n - 1]

        his_end_id = endo.global_node_ids_triangles[temp_id]
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
        position_id_his_end = np.argwhere(edges == side_his_point_ids[-1])[0]
        return (position_id_his_end, his_end_coord, new_nodes, edges, sgmt)

    def compute_left_right_bundle(self, start_coord, start_id, side: str, beam_length: float = 1.5):
        """Bundle branch."""
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

    def compute_bachman_bundle(self, start_coord, end_coord, beam_length: float = 1.5) -> BeamMesh:
        """Compute Bachman bundle conduction system."""
        la_epi = self.m.mesh.get_surface(self.m.left_atrium.epicardium.id)
        ra_epi = self.m.mesh.get_surface(self.m.right_atrium.epicardium.id)
        epi = la_epi.merge(ra_epi)

        start_id = epi.find_closest_point(start_coord)
        end_id = epi.find_closest_point(end_coord)
        path = epi.geodesic(start_id, end_id)

        #
        beam_nodes = _refine_line(path.points, beam_length=beam_length)[1:-1, :]
        point_ids = np.linspace(0, len(beam_nodes) - 1, len(beam_nodes), dtype=int)

        glob_start_id = epi.point_data["_global-point-ids"][start_id]
        glob_end_id = epi.point_data["_global-point-ids"][end_id]
        point_ids = np.insert(point_ids, 0, glob_start_id)
        point_ids = np.append(point_ids, glob_end_id)

        # build connectivity table
        edges = np.vstack((point_ids[:-1], point_ids[1:])).T

        mask = np.ones(edges.shape, dtype=bool)
        mask[0, 0] = False  # Start point at solid
        mask[-1, 1] = False  # End point at solid

        #! add surface to central mesh.
        surface_id = int(np.max(self.m.mesh.surface_ids) + 1)
        self.m.mesh.add_surface(epi, surface_id, name="Bachman segment")
        self.m.mesh = self.m.mesh.clean()

        beam_net = self.m.add_beam_net(beam_nodes, edges, mask, pid=0, name="Bachman bundle")

        return beam_net
