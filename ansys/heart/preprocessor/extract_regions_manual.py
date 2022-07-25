from matplotlib.pyplot import get
import numpy as np
import pathlib
import copy

# import tqdm as tqdm

from ansys.heart.preprocessor.vtk_module import (
    # read_vtk_unstructuredgrid_file,
    threshold_vtk_data,
    vtk_read_mesh_file,
    get_tetra_info_from_unstructgrid,
    write_vtkdata_to_vtkfile,
)


def get_faces_tetra(tetra):
    """Gets faces that make up the tetra"""
    num_tetra = tetra.shape[0]
    faces = np.zeros((num_tetra, 3, 4), dtype=int)
    masks = np.array(
        [
            [True, True, True, False],
            [True, True, False, True],
            [True, False, True, True],
            [False, True, True, True],
        ]
    )
    for ii, mask in enumerate(masks):
        faces[:, :, ii] = tetra[:, mask]

    return faces


def tetra_to_faces(tetra):
    """Creates list of unique faces from tetrahedrons and returns tetra_face_map"""
    faces = get_faces_tetra(tetra)

    # reshape to (4*NumTetra, 3)
    num_tetra = tetra.shape[0]
    faces_1 = np.reshape(faces.transpose(0, 2, 1), (4 * num_tetra, 3))

    # construct tetra_face_map. Constructs tetrahedron such that it points to face array
    tetra_face_map = np.reshape(np.arange(0, num_tetra * 4, 1), (num_tetra, 4))

    # sort faces to facilitate finding uniques
    faces_1 = np.sort(faces_1, axis=1)

    # find unique faces
    unique_faces, index, inverse, counts = np.unique(
        faces_1, axis=0, return_index=True, return_inverse=True, return_counts=True
    )

    # update tetra_face_map accordingly
    tetra_face_map = inverse[tetra_face_map]

    # find connectivity c0 an c1 (cell indices that are connected to the face)
    tet_ids = np.arange(0, num_tetra, 1)
    face_ids = np.arange(0, len(unique_faces), 1)

    # this giives faces_1 again:
    # unique_faces[tetra_face_map]

    return tetra_face_map, unique_faces


def face_tetra_connectivity(tetra: np.array):
    """Computes the tetra-face connectivity tables"""

    import time as time

    print("Establishing face connectivity...")
    t0 = time.time()

    # num_tetra = tetra.shape[0]
    # faces = np.zeros( (num_tetra, 3, 4), dtype=int )
    # masks = np.array([
    #     [ True, True, True, False],
    #     [ True, True, False, True],
    #     [ True, False, True, True],
    #     [ False, True, True, True]
    #      ] )
    # for ii, mask in enumerate( masks ):
    #     faces[:,:,ii] = tetra[:, mask]
    faces = get_faces_tetra(tetra)
    # reshape faces such that shape is (NumTetra*4, 3) with following structure:
    # i n1 n2 n3
    # 0 tet1_face1
    # 1 tet1_face2
    # 2 tet1_face3
    # 3 tet1_face4
    # 4 tet2_face1
    # 5 tet2_face2
    # 6 tet2_face3
    # 7 tet2_face4
    # ...
    num_tetra = tetra.shape[0]
    faces_1 = np.reshape(faces.transpose(0, 2, 1), (4 * num_tetra, 3))

    # sort faces in order to find duplicates
    faces_sorted = np.sort(faces_1, axis=1)
    np.sort(faces_sorted, axis=0)

    # find duplicate rows
    faces_unique, index, inverse, counts = np.unique(
        faces_sorted, return_index=True, return_inverse=True, return_counts=True, axis=0
    )

    # find duplicate rows in reversed order
    faces_sorted_flip = np.flipud(faces_sorted)
    faces_unique_r, index_r, inverse_r, counts_r = np.unique(
        faces_sorted_flip, return_index=True, return_inverse=True, return_counts=True, axis=0
    )

    # number of cells connected to each face
    num_cells_connected = counts[inverse]

    tetra_ids = np.repeat(np.arange(0, num_tetra, 1), 4)
    tetra_ids_flip = np.flipud(tetra_ids)

    # get connected tetra id for each face (two for interior, one for boundary face)
    c0 = tetra_ids[index][inverse]
    c1 = np.flip(tetra_ids_flip[index_r][inverse_r])

    # print("num cells connected")
    # print("  " + str(num_cells_connected))
    # print("c0" + str(c0))
    # print("c1" + str(c1))

    # NOTE: cell ordering makes sense in the following, but way slower.
    # c0 = tetra_ids
    # c1 = tetra_ids[inverse]
    # c1 = copy.deepcopy(c0)
    # c1[ num_cells_connected == 1 ] = c0[num_cells_connected == 1]
    # # populate c1 with a loop
    # for ii in np.where(num_cells_connected == 2)[0]:
    #     mask = np.all( faces_sorted[ii,:] == faces_sorted, axis = 1)
    #     mask[ii] = False
    #     c1[mask] = c0[ii]

    # # for all interior faces find matching cell
    # # TODO: How to get c1 efficiently without using a loop
    # print("num cells connected")
    # print("  " + str(num_cells_connected))
    # print("c0" + str(c0))
    # print("c1" + str(c1))
    t1 = time.time()
    print("Time elapsed: {:.1f}".format(t1 - t0))

    return faces_1, c0, c1


def get_face_type(faces: np.array, face_cell_connectivity: np.array) -> np.array:
    """Establishes face type. Either boundary faces or interior faces

    Parameters
    ----------
    faces : np.array
        Array with face definitions
    face_cell_connectivity : np.array
        Array describing to which cells each of the faces is connected to

    Returns
    -------
    np.array
        Type of face. Either interior (face_type = 1) or boundary (face_type = 2)
    """

    interior_face_ids = face_cell_connectivity[:, 0] != face_cell_connectivity[:, 1]
    boundary_face_ids = face_cell_connectivity[:, 0] == face_cell_connectivity[:, 1]
    face_types = np.zeros((faces.shape[0]), dtype=int)
    face_types[interior_face_ids] = 1
    face_types[boundary_face_ids] = 2
    num_assigned = np.sum(boundary_face_ids) + np.sum(interior_face_ids)
    assert num_assigned == faces.shape[0], "Not all faces assigned as either interior or boundary"
    return face_types


if __name__ == "__main__":

    # ************** Simple 2 tetrahedron example *******
    nodes = np.array([[1, 0, 0], [-1, 0, 0], [0, 0, 1], [0, -1, 0], [0, 1, 0]], dtype=float)
    tetra = np.array([[0, 1, 2, 3], [0, 1, 2, 4]])
    part_ids = np.array([1, 2])

    # get face-cell connectivity table
    faces, c0, c1 = face_tetra_connectivity(tetra)
    face_types = get_face_type(faces, np.array([c0, c1]).T)

    interior_faces = faces[face_types == 1, :]
    boundary_faces = faces[face_types == 2, :]
    interface_faces = faces[part_ids[c0] != part_ids[c1], :]

    # ************** Real example **********************

    path_mesh_file = pathlib.PurePath(
        pathlib.Path(__file__).parents[3],
        "tests",
        "heart",
        "assets",
        "cases",
        "strocchi2020",
        "01",
        "01.case",
    )
    vtk_volume = vtk_read_mesh_file(path_mesh_file)

    from ansys.heart.preprocessor.vtk_module import (
        threshold_vtk_data,
        vtk_surface_filter,
        vtk_surface_to_stl,
    )
    from ansys.heart.preprocessor.mesh_module import add_solid_name_to_stl
    import os

    # ** create stl of surface for each part:
    tags = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
    for tag in tags:
        volume_reduced_vtk = threshold_vtk_data(vtk_volume, tag, tag, "tags")[0]
        surface_vtk = vtk_surface_filter(volume_reduced_vtk)
        filename = os.path.join("..", "tmp", "part_{:02d}.stl".format(tag))
        vtk_surface_to_stl(surface_vtk, filename)
        add_solid_name_to_stl(filename, "part_{:02d}".format(tag), "binary")

    # ** full heart example:
    points, tetra, cell_data, point_data = get_tetra_info_from_unstructgrid(vtk_volume)
    # get part ids
    part_ids = cell_data["tags"]
    part_ids = np.array(part_ids, dtype=int)
    # get face-tetra connectivity
    faces, c0, c1 = face_tetra_connectivity(tetra)
    # get face types
    face_types = get_face_type(faces, np.array([c0, c1]).T)
    interior_faces = faces[face_types == 1, :]
    boundary_faces = faces[face_types == 2, :]

    # get interface between two parts
    # loop over these interface pairs:
    part_pairs = [
        [1, 14],
        [1, 16],
        [2, 15],
        [2, 17],
        [3, 14],
        [3, 18],
        [3, 19],
        [3, 20],
        [3, 21],
        [3, 22],
        [3, 7],
        [3, 8],
        [3, 9],
        [3, 10],
        [3, 11],
        [4, 15],
        [4, 23],
        [4, 24],
        [4, 12],
        [4, 13],
        [5, 16],
        [6, 17],
    ]
    parts_ids_to_extract = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
    from ansys.heart.preprocessor.vtk_module import (
        create_vtk_surface_triangles,
        write_vtkdata_to_vtkfile,
        get_connected_regions,
    )

    # write interface pairs:
    for pair in part_pairs:
        mask_interface = part_ids[c0] != part_ids[c1]
        mask_part1 = (
            np.sum(np.vstack([part_ids[c0] == pair[0], part_ids[c1] == pair[0]]), axis=0) == 1
        )
        mask_part2 = (
            np.sum(np.vstack([part_ids[c0] == pair[1], part_ids[c1] == pair[1]]), axis=0) == 1
        )

        face_ids_interface = np.where(
            np.all(np.vstack([mask_interface, mask_part1, mask_part2]), axis=0)
        )[0]
        faces_interface = faces[face_ids_interface, :]
        vtk_surface = create_vtk_surface_triangles(points, faces_interface)
        write_vtkdata_to_vtkfile(
            vtk_surface, "..\\tmp\\pair_{:02d}_{:02d}.vtk".format(pair[0], pair[1])
        )
        part_name = "interface_pair_{:02d}_{:02d}".format(*pair)
        filename = "..\\tmp\\{:}.stl".format(part_name)
        vtk_surface_to_stl(vtk_surface, filename)
        add_solid_name_to_stl(filename, part_name, "binary")

    # write boundary faces
    for part_id in parts_ids_to_extract:
        mask = np.vstack([face_types == 2, part_ids[c0] == part_id])
        mask = np.all(mask, axis=0)
        boundary_faces = faces[mask, :]

        # extract regions for parts 1-4 (LV, RV, LA, RA)
        if part_id in [1, 2, 3, 4]:
            region_ids, vtk_surface = get_connected_regions(points, boundary_faces, True)
            write_vtkdata_to_vtkfile(
                vtk_surface, "..\\tmp\\part_with_region_ids_{:02d}.vtk".format(part_id)
            )

            # smalles region endocardium, largest region epicardium. In case of Left ventricle
            # smallest region is septum
            name_map = ["endocardium", "epicardium", "septum"]
            unique_region_ids, counts = np.unique(region_ids, return_counts=True)
            for ii, index in enumerate(np.argsort(counts)):
                unique_region_ids[index]
                vtk_surface1 = create_vtk_surface_triangles(
                    points, boundary_faces[region_ids == unique_region_ids[index]]
                )
                write_vtkdata_to_vtkfile(
                    vtk_surface1,
                    "..\\tmp\\part_{:02d}_region_{:02d}_{:}.vtk".format(
                        part_id, ii, name_map[index]
                    ),
                )
                part_name = "part_{:02d}_{:}".format(part_id, name_map[index])
                filename = "..\\tmp\\{0}.stl".format(part_name)
                vtk_surface_to_stl(vtk_surface1, filename)
                add_solid_name_to_stl(filename, part_name, "binary")
        else:
            vtk_surface = create_vtk_surface_triangles(points, boundary_faces)
            write_vtkdata_to_vtkfile(vtk_surface, "..\\tmp\\part_{:02d}.vtk".format(part_id))
            part_name = "part_{:02d}".format(part_id)
            filename = "..\\tmp\\{0}.stl".format(part_name)
            vtk_surface_to_stl(vtk_surface, filename)
            add_solid_name_to_stl(filename, part_name, "binary")

    # ** bi-ventricular example
    vtk_volume_reduced = threshold_vtk_data(vtk_volume, 1, 2, "tags")[0]

    points, tetra, cell_data, point_data = get_tetra_info_from_unstructgrid(vtk_volume_reduced)
    part_ids = cell_data["tags"]
    part_ids = np.array(part_ids, dtype=int)

    # get face-cell connectivity table
    faces, c0, c1 = face_tetra_connectivity(tetra)
    face_types = get_face_type(faces, np.array([c0, c1]).T)

    interior_faces = faces[face_types == 1, :]
    boundary_faces = faces[face_types == 2, :]

    # get interface between two parts
    interface_faces = faces[part_ids[c0] != part_ids[c1], :]

    # get interface between 1 and 2
    connected_parts = np.array([part_ids[c0], part_ids[c1]]).T

    mask_part1 = np.any(np.isin(connected_parts, 1), axis=1)
    mask_part2 = np.any(np.isin(connected_parts, 2), axis=1)
    mask_interior = face_types == 1

    mask_boundary_part1 = np.all(np.array([mask_part1, face_types == 2]), axis=0)
    mask_boundary_part2 = np.all(np.array([mask_part2, face_types == 2]), axis=0)
    mask_interface_p1p2 = np.all(np.array([mask_interior, mask_part1, mask_part2]), axis=0)

    boundary_faces_p1 = faces[mask_boundary_part1]
    boundary_faces_p2 = faces[mask_boundary_part2]
    interface_p1p2 = np.unique(np.sort(faces[mask_interface_p1p2], axis=1), axis=0)

    import meshio

    cells = [
        ("triangle", boundary_faces_p1),
        ("triangle", boundary_faces_p2),
        ("triangle", interface_p1p2),
    ]
    mesh = meshio.Mesh(points=points, cells=cells)
    mesh.write("tmp/p1_p2.vtk")
    mesh.write("tmp/p1_p2.stl")

    import meshio

    cells = [
        ("triangle", boundary_faces_p1),
    ]
    mesh = meshio.Mesh(points=points, cells=cells)
    mesh.write("tmp/boundary_p1.stl")

    cells = [
        ("triangle", boundary_faces_p2),
    ]
    mesh = meshio.Mesh(points=points, cells=cells)
    mesh.write("tmp/boundary_p2.stl")

    cells = [
        ("triangle", interface_p1p2),
    ]
    mesh = meshio.Mesh(points=points, cells=cells)
    mesh.write("tmp/interface_p1p2.stl")

    # write to vtk file
    write_to_file = True

    interface_faces1 = np.unique(np.sort(interface_faces, axis=1), axis=0)
    if write_to_file:
        import meshio

        cells = [("triangle", boundary_faces)]
        mesh = meshio.Mesh(points, cells)
        mesh.write("boundary_faces1.vtk")

        # write interface faces

        cells = [("triangle", interface_faces1)]
        # cell_data = {"face_type" : [face_type] }
        mesh = meshio.Mesh(
            points,
            cells
            # cell_data = cell_data
        )
        mesh.write("interface_faces.vtk")

    print("Protected")
