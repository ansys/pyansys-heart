"""
Module that contains classes relevant for the mesh.

Such as a Mesh object, Part object, Features, etc.

"""
import copy
import logging
import pathlib
from typing import List, Optional, Tuple, Union

LOGGER = logging.getLogger("pyheart_global.preprocessor")
import ansys.heart.preprocessor.mesh.connectivity as connect
import ansys.heart.preprocessor.mesh.geodisc as geodisc
import ansys.heart.preprocessor.mesh.vtkmethods as vtkmethods
import numpy as np

try:
    import pyvista as pv
except ImportError:
    LOGGER.warning("Importing pyvista failed. Install with: pip install pyvista")


class Feature:
    """Feature class."""

    def __init__(self, name: str = None) -> None:
        self.name = name
        """Name of feature."""
        self.type = None
        """Type of feature."""
        self.nsid: int = None
        """Node set id associated with feature."""
        self.pid: int = None
        """Part id associated with the feature."""

        pass


class BoundaryEdges(Feature):
    """Edges class."""

    def __init__(self, edges: np.ndarray = None) -> None:
        super().__init__()
        self.type = "edges"
        """Feature type."""

        self.node_ids = None
        """List of edges."""
        self.groups = None
        """Grouped edges based on connectivity."""


class EdgeGroup:
    """Edge group class, contains info on connected edges."""

    def __init__(self, edges: np.ndarray = None, type: str = None) -> None:
        self.edges = edges
        """Edges in edge group."""
        self.type = type
        """Type of edge group: 'open' or 'closed'."""
        pass


class SurfaceMesh(pv.PolyData, Feature):
    """Surface class."""

    @property
    def nodes(self):
        """Node coordinates."""
        return np.array(self.points)

    @nodes.setter
    def nodes(self, array: np.ndarray):
        if isinstance(array, type(None)):
            return
        try:
            self.points = array
        except:
            LOGGER.warning("Failed to set points.")
            return

    @property
    def triangles(self):
        """Triangular faces of the surface num_faces x 3."""
        faces = np.reshape(self.faces, (self.n_cells, 3 + 1))[:, 1:]
        return faces

    @triangles.setter
    def triangles(self, value: np.ndarray):
        # sets faces of PolyData
        try:
            num_faces = value.shape[0]
            faces = np.hstack([np.full((num_faces, 1), 3, dtype=np.int8), value])
            self.faces = faces
        except:
            return

    def __init__(
        self,
        name: str = None,
        triangles: np.ndarray = None,
        nodes: np.ndarray = None,
        id: int = None,
    ) -> None:
        super().__init__(self)
        Feature.__init__(self, name)

        self.type = "surface"
        """Surface type."""
        self.boundary_edges: np.ndarray = np.empty((0, 2), dtype=int)
        """Boundary edges."""
        self.edge_groups: List[EdgeGroup] = []
        """Edge groups."""
        self.id: int = id
        """ID of surface."""
        self.nsid: int = None
        """ID of corresponding set of nodes."""
        # self.cell_data: dict = {}
        # """Data associated with each face/cell of surface."""
        # self.point_data: dict = {}
        """Data associated with each point on surface."""
        self._pv_polydata: pv.PolyData = pv.PolyData()
        """Pyvista representation of the surface mesh."""

        self.triangles = triangles
        """Triangular faces of the surface num_faces x 3."""
        self.nodes = nodes
        """Node coordinates."""

    @property
    def node_ids(self) -> np.ndarray:
        """Global node ids - sorted by earliest occurrence."""
        _, idx = np.unique(self.triangles.flatten(), return_index=True)
        node_ids = self.triangles.flatten()[np.sort(idx)]
        return node_ids

    @property
    def boundary_nodes(self) -> np.ndarray:
        """Global node ids of nodes on the boundary of the mesh (if any)."""
        _, idx = np.unique(self.boundary_edges.flatten(), return_index=True)
        node_ids = self.boundary_edges.flatten()[np.sort(idx)]
        return node_ids

    def compute_centroid(self) -> np.ndarray:
        """Compute the centroid of the surface."""
        return np.mean(self.nodes[np.unique(self.triangles), :], axis=0)

    def compute_bounding_box(self) -> Tuple[np.ndarray, float]:
        """Compute the bounding box of the surface."""
        bounding_box = np.reshape(self.clean().bounds, (3, 2)).T
        volume = np.prod(np.diff(bounding_box, axis=0))
        return bounding_box, volume

    def get_boundary_edges(self, append_triangles=None) -> List[EdgeGroup]:
        """
        Get boundary edges (if any) of the surface and groups them by connectivity.

        Parameters
        ----------
        append_triangles: optional
            special fix for right ventricle endocardium surface since it needs one part
            from spetum.
        """
        write_vtk = False

        if append_triangles is not None:
            self.boundary_edges = connect.get_free_edges(
                np.vstack((self.triangles, append_triangles))
            )
        else:
            self.boundary_edges = connect.get_free_edges(self.triangles)

        edge_groups, group_types = connect.edge_connectivity(
            self.boundary_edges, return_type=True, sort_closed=True
        )

        for ii, edge_group in enumerate(edge_groups):
            group = EdgeGroup(edges=edge_group, type=group_types[ii])
            self.edge_groups.append(group)

            if write_vtk:
                tris = np.vstack([edge_group[:, 0], edge_group.T]).T
                vtk_surf = vtkmethods.create_vtk_surface_triangles(self.nodes, tris)
                vtkmethods.write_vtkdata_to_vtkfile(
                    vtk_surf, "edges_{0}_{1}.vtk".format(ii, self.name)
                )

        return self.edge_groups

    def separate_connected_regions(self):
        """Use vtk to get connected regions and separate into different surfaces."""
        region_ids = vtkmethods.get_connected_regions(self.nodes, self.triangles)

        return region_ids

    def smooth_boundary_edges(self, window_size: int = 3) -> np.ndarray:
        """Smooth the boundary edges if they are closed."""
        if window_size % 2 != 1:
            raise ValueError("Please specify window size to be an uneven number")

        # self.write_feature_edges_to_vtk("unsmoothed")
        modified_nodes = np.empty((0), dtype=int)
        for ii, edge_group in enumerate(self.edge_groups):
            edges = edge_group.edges

            idx = np.unique(edges.flatten(), return_index=True)[1]
            idx = np.sort(idx)  # maintain order
            node_ids = edges.flatten()[idx]
            nodes = copy.deepcopy(self.nodes[node_ids, :])
            # project points onto plane
            nodes = geodisc.project_3d_points(nodes)[0]

            # self.write_feature_edges_to_vtk("projected")

            if edge_group.type == "closed":
                # use a window average method to smooth edges nodes
                num_points_to_add = int((window_size - 1) / 2)
                nodes = np.concatenate(
                    (nodes[-num_points_to_add:], nodes, nodes[0:num_points_to_add])
                )
                offset = num_points_to_add
                for ii, node in enumerate(nodes[:-num_points_to_add]):
                    nodes[ii + offset] = np.mean(nodes[ii : ii + window_size, :], axis=0)
                nodes = nodes[num_points_to_add:-num_points_to_add]

            self.points[node_ids, :] = nodes

            modified_nodes = np.append(modified_nodes, node_ids)

        # self.write_feature_edges_to_vtk("smoothed")

        return modified_nodes

    def write_to_stl(self, filename: pathlib.Path = None) -> None:
        """Write the surface to a vtk file."""
        if not filename:
            filename = "_".join(self.name.lower().split()) + ".stl"
        if filename[-4:] != ".stl":
            filename = filename + ".stl"

        # NOTE: The below should yield the same stls, but somehow fluent meshing
        # produces a slightly different mesh. Should still be valid though
        # cleaned = self.clean()
        # cleaned.save(filename)
        # vtkmethods.add_solid_name_to_stl(filename, self.name, file_type="binary")

        vtk_surface = vtkmethods.create_vtk_surface_triangles(self.nodes, self.triangles)
        vtkmethods.vtk_surface_to_stl(vtk_surface, filename, self.name)
        return

    def write_feature_edges_to_vtk(self, prefix: str = None, per_edge_group: bool = False) -> None:
        """Write the feature edges to a vtk file."""
        edges = np.array((0, 2))
        for ii, edge_group in enumerate(self.edge_groups):
            edges = np.vstack([edges, edge_group.edges])
            if per_edge_group:
                tris = np.vstack([edge_group.edges[:, 0], edge_group.edges.T]).T
                vtk_surf = vtkmethods.create_vtk_surface_triangles(self.nodes, tris)
                filename = "{0}_groupid_{1}_edges_{2}.vtk".format(prefix, ii, self.name)
                vtkmethods.write_vtkdata_to_vtkfile(vtk_surf, filename)
        if not per_edge_group:
            tris = np.vstack([edges[:, 0], edges.T]).T
            vtk_surf = vtkmethods.create_vtk_surface_triangles(self.nodes, tris)
            filename = "{0}_groupid_edges_{1}.vtk".format(prefix, self.name)
            vtkmethods.write_vtkdata_to_vtkfile(vtk_surf, filename)

        return

    def _to_pyvista_object(self) -> pv.PolyData:
        """Convert to pyvista polydata object.

        Returns
        -------
        pv.PolyData
            pyvista PolyData object
        """
        DeprecationWarning("_to_pyvista_object is deprecated.")
        faces = np.hstack([np.ones((self.triangles.shape[0], 1), dtype=int) * 3, self.triangles])
        nodes = self.nodes
        faces = np.reshape(faces, (faces.size))
        polydata = pv.PolyData(nodes, faces)
        if self.cell_data:
            for key, value in self.cell_data.items():
                polydata.cell_data[key] = value

        if self.point_data:
            for key, value in self.point_data.items():
                polydata.point_data[key] = value

        return polydata


class BeamMesh(pv.UnstructuredGrid, Feature):
    """Beam class."""

    all_beam_nodes = []
    # beam nodes array

    @property
    def nodes(self):
        """Node coordinates."""
        return np.array(self.points)

    @nodes.setter
    def nodes(self, array: np.ndarray):
        if isinstance(array, type(None)):
            return
        try:
            self.points = array
        except:
            LOGGER.warning("Failed to set nodes.")
            return

    @property
    def edges(self):
        """Tetrahedrons num_tetra x 4."""
        return self.cells_dict[pv.CellType.LINE]

    @edges.setter
    def edges(self, value: np.ndarray):
        # sets lines of UnstructuredGrid
        try:
            points = self.points
            celltypes = np.full(value.shape[0], pv.CellType.LINE, dtype=np.int8)
            lines = np.hstack([np.full(len(celltypes), 2)[:, None], value])
            super().__init__(lines, celltypes, points)
        except:
            LOGGER.warning("Failed to set lines.")
            return

    def __init__(
        self,
        name: str = None,
        edges: np.ndarray = None,
        nodes: np.ndarray = None,
        beam_nodes_mask: np.ndarray = None,
        pid: int = None,
        nsid: int = -1,
    ) -> None:
        super().__init__(self)
        Feature.__init__(self, name)

        self.edges = edges
        """Beams edges."""

        self.nodes = nodes
        """Node coordinates."""

        self.beam_nodes_mask = beam_nodes_mask
        """True for beam nodes, False for solid nodes."""

        self.pid = pid
        """Part id associated with the network."""

        self.nsid: int = nsid
        """Surface id associated with the network."""

        self._all_beam_nodes: np.ndarray = np.empty((0, 3))
        """Temporary attribute to save all previously created beam nodes."""


class Cavity(Feature):
    """Cavity class."""

    def __init__(self, surface: SurfaceMesh = None, centroid: np.ndarray = None, name=None) -> None:
        super().__init__(name)
        self.type = "cavity"
        """Type."""
        self.surface: SurfaceMesh = surface
        """Surface mesh making up the cavity."""
        self.centroid: np.ndarray = centroid
        """Centroid of the cavity."""

    @property
    def volume(self):
        """Volume of the cavity."""
        return self.surface.volume

    def compute_volume(self) -> float:
        """Compute the volume of the (enclosed) cavity.

        Notes
        -----
        - Writes stl and computes volume from stl
        - Assumes normals are pointing inwards
        """
        DeprecationWarning(
            "compute_volume() is deprecated. Use the volume property of this class instead."
        )
        return self.surface.volume

    def compute_centroid(self):
        """Compute the centroid of the cavity."""
        # self.centroid = np.mean(self.surface.nodes[np.unique(self.surface.triangles), :], axis=0)
        self.centroid = self.surface.center
        return self.centroid


class Cap(Feature):
    """Cap class."""

    def __init__(self, name: str = None, node_ids: Union[List[int], np.ndarray] = []) -> None:
        super().__init__(name)
        self.node_ids = node_ids
        """(Global) node ids of the cap."""
        self.triangles = None
        """Triangulation of cap."""
        self.normal = None
        """Normal of cap."""
        self.centroid = None
        """Centroid of cap."""
        self.centroid_id = None
        """Centroid of cap ID (in case centroid node is created)."""
        self.type = "cap"
        """Type."""
        return

    def tessellate(self, center_point_id=None) -> np.ndarray:
        """
        Form triangles with the node ids.

        Parameters
        ----------
        center_point_id: ID of the center point of cap

        Returns
        -------
        Mesh connectivity of cap (triangles)

        """
        if center_point_id is None:
            tris = []
            for ii, _ in enumerate(self.node_ids[0:-2]):
                # first node is reference node
                tri = [self.node_ids[0], self.node_ids[ii + 1], self.node_ids[ii + 2]]
                tris.append(tri)
            self.triangles = np.array(tris, dtype=int)
        else:
            ref_node = center_point_id[0]
            num_triangles = self.node_ids.shape[0] + 1
            tris = [[ref_node, self.node_ids[0], self.node_ids[1]]]
            for ii, _ in enumerate(self.node_ids[0:-2]):
                tri = [ref_node, self.node_ids[ii + 1], self.node_ids[ii + 2]]
                tris.append(tri)
            tris.append([ref_node, self.node_ids[-1], self.node_ids[0]])
            self.triangles = np.array(tris, dtype=int)

        return self.triangles


class Point(Feature):
    """Point class. Can be used to collect relevant points in the mesh."""

    def __init__(self, name: str = None, xyz: np.ndarray = None, node_id: int = None) -> None:
        super().__init__(name)

        self.xyz: np.ndarray = xyz
        """XYZ Coordinates of point."""
        self.node_id: int = node_id
        """Global node id of point."""


class Mesh(pv.UnstructuredGrid):
    """Mesh class: inherits from pyvista UnstructuredGrid.

    Notes
    -----
    Only tetrahedrons are supported.
    Additional attributes are added on top of the pyvista UnstructuredGrid class
    """

    @property
    def nodes(self):
        """Node coordinates."""
        return np.array(self.points)

    @nodes.setter
    def nodes(self, array: np.ndarray):
        if isinstance(array, type(None)):
            return
        try:
            self.points = array
        except:
            LOGGER.warning("Failed to set nodes.")
            return

    @property
    def tetrahedrons(self):
        """Tetrahedrons num_tetra x 4."""
        return self.cells_dict[pv.CellType.TETRA]

    @tetrahedrons.setter
    def tetrahedrons(self, value: np.ndarray):
        # sets tetrahedrons of UnstructuredGrid
        try:
            points = self.points
            celltypes = np.full(value.shape[0], pv.CellType.TETRA, dtype=np.int8)
            tetra = np.hstack([np.full(len(celltypes), 4)[:, None], value])
            super().__init__(tetra, celltypes, points)
        except:
            LOGGER.warning("Failed to set tetrahedrons.")
            return

    def __init__(self, *args):
        super().__init__(*args)

        # self.tetrahedrons: np.ndarray = None
        # """Tetrahedral volume elements of the mesh."""
        # self.nodes: np.ndarray = None
        # """Nodes of the mesh."""
        # self.cell_data: dict = None
        # """Data per mesh cell/element."""
        # self.point_data: dict = None
        # """Data per mesh point."""
        self.triangles: np.ndarray = None
        """Faces that make up the tetrahedrons."""
        self.face_types: np.ndarray = None
        """Type of face: 1: interior face, 2: boundary face, 3: interface face."""
        self.conn = {"c0": [], "c1": []}
        """Face-tetra connectivity array."""
        # TODO: just store used nodes in interfaces and boundaries
        # and add mapper to map from local to global (volume mesh) node ids
        self.interfaces: List[SurfaceMesh] = []
        """List of surface meshes that make up the interface between different parts."""
        self.boundaries: List[SurfaceMesh] = []
        """List of boundary surface meshes within the part."""
        pass

    @property
    def part_ids(self) -> np.ndarray:
        """Array of part ids indicating to which part the tetrahedron belongs.

        Notes
        -----
        This is derived from the "tags" field in cell data
        """
        try:
            value = self.cell_data["tags"].astype(int)
        except (KeyError, NameError):
            LOGGER.warning("'tags' field not found in self.cell_data")
            value = None
        return value

    @property
    def boundary_names(self) -> List[str]:
        """Iterate over boundaries and returns their names."""
        return [b.name for b in self.boundaries]

    def read_mesh_file(self, filename: pathlib.Path) -> None:
        """Read mesh file."""
        mesh = pv.read(filename)
        # .case gives multiblock
        if isinstance(mesh, pv.MultiBlock):
            mesh: pv.UnstructuredGrid = mesh.GetBlock(0)

        if not isinstance(mesh, pv.UnstructuredGrid):
            LOGGER.warning("Failed to read mesh file. Expecting .vtk unstructured grid or .case")
            return

        self.points = mesh.points
        self.tetrahedrons = mesh.cells_dict[pv.CellType.TETRA]
        for key, value in mesh.cell_data.items():
            self.cell_data[key] = mesh.cell_data[key]
        for key, value in mesh.point_data.items():
            self.point_data[key] = mesh.point_data[key]

        return

    def read_mesh_file_rodero2021(self, filename: pathlib.Path) -> None:
        """Read mesh file - but modifies the fields to match data of Strocchi 2020."""
        mesh = pv.read(filename)
        # .case gives multiblock
        if isinstance(mesh, pv.MultiBlock):
            mesh: pv.UnstructuredGrid = mesh.GetBlock(0)

        if not isinstance(mesh, pv.UnstructuredGrid):
            LOGGER.warning("Failed to read mesh file. Expecting .vtk unstructured grid or .case")
            return

        name_array_mapping = [
            ["tags", "ID", "cell"],
            ["fiber", "fibres", "cell"],
            ["sheet", "sheets", "cell"],
            ["uvc_longitudinal", "Z.dat", "point"],
            ["uvc_rotational", "PHI.dat", "point"],
            ["uvc_transmural", "RHO.dat", "point"],
            ["uvc_intraventricular", "V.dat", "point"],
        ]

        # rename tags in rodero
        for item in name_array_mapping:
            mesh.rename_array(item[1], item[0], item[2])

        self.points = mesh.points
        self.tetrahedrons = mesh.cells_dict[pv.CellType.TETRA]
        for key, value in mesh.cell_data.items():
            self.cell_data[key] = mesh.cell_data[key]
        for key, value in mesh.point_data.items():
            self.point_data[key] = mesh.point_data[key]

        if np.issubdtype(self.cell_data["tags"].dtype, np.integer):
            self.cell_data["tags"] = np.array(self.cell_data["tags"], dtype=float)

        return None

    def write_to_vtk(self, filename: pathlib.Path) -> None:
        """Write mesh to VTK file."""
        self.save(filename)
        return

    def keep_elements_with_value(
        self,
        values: List[int],
        field_name: str,
    ) -> None:
        """Remove elements that satisfy a certain cell value of a specific field."""
        mask = np.isin(self.cell_data[field_name], values)
        self.tetrahedrons = self.tetrahedrons[mask, :]
        for key in self.cell_data.keys():
            self.cell_data[key] = self.cell_data[key][mask]
        return

    def establish_connectivity(self) -> None:
        """Establish the connetivity of the tetrahedrons."""
        self.triangles, self.conn["c0"], self.conn["c1"] = connect.face_tetra_connectivity(
            self.tetrahedrons
        )
        # get the face types
        c0c1_matrix = np.array([self.conn["c0"], self.conn["c1"]]).transpose()
        self.face_types = connect.get_face_type(self.triangles, c0c1_matrix)

        return

    def get_mask_interface_faces(
        self, return_pairs: bool = False
    ) -> Tuple[np.ndarray, Optional[List[int]]]:
        """Get the (interface) faces between two parts."""
        c0 = self.conn["c0"]
        c1 = self.conn["c1"]
        part_ids = self.part_ids
        face_types = self.face_types

        mask_interface_faces = np.all(
            np.vstack([part_ids[c0] != part_ids[c1], face_types == 1]), axis=0
        )

        # get interface pairs
        interface_pairs = np.vstack(
            [part_ids[c0][mask_interface_faces], part_ids[c1][mask_interface_faces]]
        )
        interface_pairs = np.array(interface_pairs, dtype=int)
        interface_pairs_sorted = np.sort(interface_pairs, axis=0)
        # get unique pairs
        part_pairs = np.unique(interface_pairs_sorted, axis=1).transpose()
        part_pairs = part_pairs.tolist()

        # mark interface faces in type array
        self.face_types[mask_interface_faces] = 3

        if return_pairs:
            return mask_interface_faces, part_pairs
        else:
            return mask_interface_faces

    def add_interfaces(
        self,
        pairs: List[List[int]],
        pair_names: List[str],
    ) -> None:
        """Add the interfaces between the parts to the mesh."""
        part_ids = self.part_ids
        c0 = self.conn["c0"]
        c1 = self.conn["c1"]

        # loop over pairs and add to the list of interfaces
        for ii, pair in enumerate(pairs):
            name = pair_names[ii]
            part_mask1 = (
                np.sum(
                    np.array(
                        [
                            part_ids[c0] == pair[0],  # mask part 1
                            part_ids[c1] == pair[0],  # mask part 1
                        ]
                    ),
                    axis=0,
                )
                == 1
            )
            part_mask2 = (
                np.sum(
                    np.array(
                        [
                            part_ids[c0] == pair[1],  # mask part 2
                            part_ids[c1] == pair[1],  # mask part 2
                        ]
                    ),
                    axis=0,
                )
                == 1
            )
            pair_mask = np.all(np.array([part_mask1, part_mask2]), axis=0)

            faces = self.triangles[pair_mask, :]
            # NOTE: Nodes are shallow copied
            self.interfaces.append(SurfaceMesh(name, faces, self.nodes))

    def smooth_interfaces(self) -> None:
        """Smooth the interfaces between the different parts."""
        for interface in self.interfaces:
            interface.get_boundary_edges()
            node_ids_smoothed = interface.smooth_boundary_edges()
            # make sure nodes of the (volume) mesh are updated
            self.points[node_ids_smoothed, :] = interface.nodes[node_ids_smoothed, :]
        return

    def add_boundaries(self, add_part_ids: List[int] = [], boundary_names: List[str] = []) -> None:
        """Add boundary surfaces to the mesh object. One surface per part."""
        part_ids = self.part_ids
        c0 = self.conn["c0"]
        c1 = self.conn["c1"]

        for ii, part_id in enumerate(add_part_ids):
            boundary_mask = np.all(
                np.array([part_ids[c0] == part_id, part_ids[c1] == part_id, self.face_types == 2]),
                axis=0,
            )
            boundary_faces = self.triangles[boundary_mask, :]
            self.boundaries.append(SurfaceMesh(boundary_names[ii], boundary_faces, self.nodes))

        return

    def get_surface_from_name(self, name: str = None):
        """Return a list of surfaces that match the given list of names.

        Notes
        -----
        Returns single surface. When multiple matches are found returns list of surfaces
        """
        surfaces_search = self.boundaries + self.interfaces
        surfaces = [s for s in surfaces_search if s.name == name]
        if len(surfaces) == 0:
            return None
        if len(surfaces) == 1:
            return surfaces[0]
        else:
            return surfaces

    def _to_pyvista_object(self) -> pv.UnstructuredGrid:
        """Convert mesh object into pyvista unstructured grid object.

        Returns
        -------
        pv.UnstructuredGrid
            Pyvista unstructured grid object.
        """
        num_tets = self.tetrahedrons.shape[0]
        cells = np.hstack(
            [np.ones((self.tetrahedrons.shape[0], 1), dtype=int) * 4, self.tetrahedrons]
        ).flatten()
        celltypes = np.ones(num_tets, dtype=int) * pv.CellType.TETRA
        points = self.nodes
        grid = pv.UnstructuredGrid(cells, celltypes, points)
        # add cell and point data
        if self.cell_data:
            for key, value in self.cell_data.items():
                if value.size == value.shape[0]:
                    grid.cell_data.set_scalars(name=key, scalars=value)
                elif len(value.shape) > 1:
                    grid.cell_data.set_vectors(name=key, vectors=value)

        if self.point_data:
            for key, value in self.point_data.items():
                if value.size == value.shape[0]:
                    grid.point_data.set_array(name=key, data=value)
                elif len(value.shape) > 1:
                    grid.point_data.set_vectors(name=key, vectors=value)

        return grid


class Part:
    """Part class."""

    @property
    def _features(self) -> List[Feature]:
        """Return list of part features."""
        features = []
        for key, value in self.__dict__.items():
            attribute = getattr(self, key)
            if isinstance(attribute, Feature):
                features.append(attribute)
        return features

    @property
    def surfaces(self) -> List[SurfaceMesh]:
        """List of surfaces belonging to part."""
        surfaces = []
        for key, value in self.__dict__.items():
            if isinstance(value, SurfaceMesh):
                surfaces.append(value)
        return surfaces

    @property
    def surface_names(self) -> List[str]:
        """List of surface names belonging to part."""
        surface_names = []
        for key, value in self.__dict__.items():
            if isinstance(value, SurfaceMesh):
                surface_names.append(value.name)
        return surface_names

    def get_point(self, pointname: str) -> Point:
        """Get point from part."""
        for point in self.points:
            if point.name == pointname:
                return point
        LOGGER.error("Cannot find point {0:s}.".format(pointname))
        return None

    def __init__(self, name: str = None, part_type: str = None) -> None:
        self.name = name
        """Name of the part."""
        self.pid = None
        """Part ID."""
        self.mid = None
        """Material id associated with part."""
        self.part_type = part_type
        """Type of the part."""
        self.tag_labels = None
        """VTK tag labels used in this part."""
        self.tag_ids = None
        """VTK tag ids used in this part."""
        self.element_ids: np.ndarray = np.empty((0, 4), dtype=int)
        """Array holding element ids that make up this part."""
        self.points: List[Point] = []
        """Points of interest belonging to the part."""
        self.caps: List[Cap] = []
        """List of caps belonging to the part."""
        self.cavity: Cavity = None

        self.has_fiber: bool = False
        """If this part has fiber/sheet data."""
        self.is_active: bool = False
        """If active stress will be established."""

        """Cavity belonging to the part."""
        if self.part_type in ["ventricle"]:
            self.apex_points: List[Point] = []
            """Points on apex."""

        self._add_surfaces()

    def _add_surfaces(self):
        """Add surfaces to the part."""
        if self.part_type in ["ventricle", "atrium"]:
            self.endocardium = SurfaceMesh("{0} endocardium".format(self.name))
            """Endocardium."""
            self.epicardium = SurfaceMesh("{0} epicardium".format(self.name))
            """Epicardium."""
            if self.part_type == "ventricle":
                self.septum = SurfaceMesh("{0} septum".format(self.name))
                """Septum surface."""
        elif self.part_type in ["artery"]:
            self.wall = SurfaceMesh("{0} wall".format(self.name))
            """Wall."""
        return

    def _add_myocardium_part(self):
        self.myocardium = Part(name="myocardium", part_type="myocardium")
        return

    def _add_septum_part(self):
        self.septum = Part(name="septum", part_type="septum")
        return

    # def get_mesh(self, mesh: Mesh = None) -> Mesh:
    #     tets: np.ndarray = mesh.tetrahedrons
    #     tets = tets[self.element_ids :]
    #     partmesh = Mesh()
    #     partmesh.tetrahedrons = tets

    #     partmesh.nodes =
    #     return partmesh
