"""Module containing class for creating conduxtion system."""
import os

heart_version = os.getenv("ANSYS_HEART_MODEL_VERSION")
if not heart_version:
    heart_version = "v0.1"

if heart_version == "v0.2":
    from ansys.heart.preprocessor.models.v0_2.models import FourChamber
elif heart_version == "v0.1":
    from ansys.heart.preprocessor.models.v0_1.models import FourChamber

import logging

from ansys.heart.preprocessor.mesh.objects import BeamMesh, Point, _create_line
import pyvista as pv

LOGGER = logging.getLogger("pyheart_global.preprocessor")
import numpy as np


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

        SinoAtrial node is defined on endocardium surface and
        between sup vena cave and inf vena cave.
        """
        if target_coord == None:
            for cap in self.m.right_atrium.caps:
                if "superior" in cap.name:
                    sup_vcava_centroid = cap.centroid
                elif "inferior" in cap.name:
                    inf_vcava_centroid = cap.centroid

            # define SinoAtrial node:
            target_coord = sup_vcava_centroid - (inf_vcava_centroid - sup_vcava_centroid) / 2

        right_atrium_endo = self.m.right_atrium.endocardium

        target_id = pv.PolyData(
            self.m.mesh.nodes[right_atrium_endo.node_ids, :]
        ).find_closest_point(target_coord)

        SA_node_id = right_atrium_endo.node_ids[target_id]

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
        right_atrium_endo = self.m.right_atrium.endocardium

        if target_coord is None:
            for surface in self.m.right_ventricle.surfaces:
                if "endocardium" in surface.name and "septum" in surface.name:
                    right_septum = surface
            # define AtrioVentricular as the closest point to septum
            target_id = pv.PolyData(
                self.m.mesh.points[right_atrium_endo.node_ids, :]
            ).find_closest_point(right_septum.center)

        else:
            target_id = pv.PolyData(
                self.m.mesh.points[right_atrium_endo.node_ids, :]
            ).find_closest_point(target_coord)

        # assign a point
        av_id = right_atrium_endo.node_ids[target_id]
        AV_point = Point(name="AV_node", xyz=self.m.mesh.nodes[av_id, :], node_id=av_id)

        self.m.right_atrium.points.append(AV_point)

        return AV_point

    def compute_av_conduction(self, beam_length: float = 1.5) -> BeamMesh:
        """Compute Atrio-Ventricular conduction by means of beams following a geodesic path."""
        right_atrium_endo = self.m.right_atrium.endocardium

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

        path_SAN_AVN = right_atrium_endo.geodesic(SA_id, AV_id)
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
        new_nodes = np.array(
            [
                self.m.right_atrium.get_point("AV_node").xyz,
                bifurcation_coord,
            ]
        )

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

        position_id_His_end_left, his_end_left_coord, new_nodes, edges = self.m._create_His_side(
            side="left",
            new_nodes=new_nodes,
            edges=edges,
            beam_length=beam_length,
            bifurcation_coord=bifurcation_coord,
            bifurcation_id=bifurcation_id,
        )
        position_id_His_end_right, his_end_right_coord, new_nodes, edges = self._create_His_side(
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

        return Point(
            xyz=his_end_left_coord,
            node_id=beam_net.edges[position_id_His_end_left[0], position_id_His_end_left[1]],
        ), Point(
            xyz=his_end_right_coord,
            node_id=beam_net.edges[position_id_His_end_right[0], position_id_His_end_right[1]],
        )

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

        side_his = np.array([bifurcation_coord, his_end_coord])
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
        return (
            position_id_His_end,
            his_end_coord,
            new_nodes,
            edges,
        )

    def compute_left_right_bundle(self, start_coord, start_id, side: str, beam_length: float = 1.5):
        """Bundle brunch."""
        if side == "Left":
            ventricle = self.m.left_ventricle
            endo_surface = self.m.left_ventricle.endocardium
        elif side == "Right":
            ventricle = self.m.right_ventricle
            face = np.hstack(
                (self.m.right_ventricle.endocardium.faces, self.m.right_ventricle.septum.faces)
            )
            endo_surface = pv.PolyData(self.m.mesh.points, face)

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

    def compute_Bachman_bundle(self):
        """Compute Bachman bundle conduction system."""
        raise NotImplementedError
        return


if __name__ == "__main__":
    model: FourChamber = FourChamber.load_model(
        r"D:\ansysdev\pyansys-heart\downloads\Strocchi2020\01\FourChamber\heart_model.pickle"
    )

    test = ConductionSystem(model)
    test.compute_SA_node([-48, 136, 393])
    test.compute_AV_node()

    test.compute_av_conduction()

    # a,b =model.compute_His_conduction()
    a, b = test.compute_His_conduction()
    print(a.xyz, a.node_id)
    print(b.xyz, b.node_id)

    test.compute_left_right_bundle(a.xyz, a.node_id, side="Left")
    test.compute_left_right_bundle(b.xyz, b.node_id, side="Right")
    model.plot_purkinje()
