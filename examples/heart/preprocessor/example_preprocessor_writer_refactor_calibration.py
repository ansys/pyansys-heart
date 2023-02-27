"""Example to pre-process (morphed) data from Casis.

Note
----
Uses 'protected' methods from the HeartModel class to chain
various operations. Moreover, some non-generic steps are required such
as compute uvc, providing a list of surfaces instead of an entire model, etc

"""
import os
import pathlib
from typing import List

from ansys.heart.preprocessor.mesh.objects import SurfaceMesh
import ansys.heart.preprocessor.mesh.vtkmethods as vtkmethods
import ansys.heart.preprocessor.models as models
import ansys.heart.writer.dynawriter as writers
import numpy as np
from vtk.numpy_interface import dataset_adapter as dsa  # type: ignore # noqa

if __name__ == "__main__":
    """BiVentricle morphed to Casis data."""
    os.environ["PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION"] = "python"
    os.environ["REMOTING_SERVER_ADDRESS"] = "localhost"
    # merge two vtks
    path_to_lv = os.path.join(
        "data",
        "endoepi_LV_tagged.vtk",
    )
    path_to_rv = os.path.join(
        "data",
        "Calib_RV_volume_array.vtk",
    )
    folder_name = pathlib.Path(path_to_lv).parent

    preproc = True
    if preproc:
        vtk_lv = vtkmethods.read_vtk_polydata_file(path_to_lv)
        vtk_rv = vtkmethods.read_vtk_polydata_file(path_to_rv)

        # remove capitalization of cell-tags
        # vtk_rv = vtkmethods.rename_vtk_array(vtk_rv, "Cell-tags", "cell-tags")
        # vtkmethods.write_vtkdata_to_vtkfile(vtk_rv, path_to_rv)

        # Merge data of both files
        points_1, tris_1, cell_data_1, point_data_1 = vtkmethods.get_tri_info_from_polydata(vtk_lv)
        points_2, tris_2, cell_data_2, point_data_2 = vtkmethods.get_tri_info_from_polydata(vtk_rv)

        points = np.vstack([points_1, points_2])
        tris = np.vstack([tris_1, tris_2 + points_1.shape[0]])

        # append cell and point data
        cell_data = {}
        for key in cell_data_1.keys():
            cell_data[key] = np.append(cell_data_1[key], cell_data_2[key])

        # point_data = {}
        # for key in point_data_1.keys():
        #     point_data[key] = np.append(point_data_1[key], point_data_2[key])

        label_to_id = {
            "Left ventricle endocardium": 1,
            "Left ventricle epicardium": 2,
            "Left ventricle myocardium mitral valve": 3,
            "Left ventricle myocardium aortic valve": 4,
            "Right ventricle endocardium": 5,
            "Right ventricle epicardium": 6,
            "Right ventricle myocardium pulmonary valve": 7,
            "Right ventricle myocardium tricuspid valve": 8,
            "Right ventricle septum": 9,
        }
        id_to_label = dict((v, k) for k, v in label_to_id.items())

        label_to_id_part = {
            "Left ventricle myocardium": 1,
            "Right ventricle myocardium": 2,
        }
        id_to_label_part = dict((v, k) for k, v in label_to_id_part.items())

        # Instantiate model
        info = models.ModelInfo(
            "CasisCalibration",
            work_directory=os.path.join(pathlib.Path(path_to_lv).parent, "BiVentricle"),
            path_to_case="",
            mesh_size=1.5,
        )
        info.create_workdir()
        info.clean_workdir()
        model = models.BiVentricle(info)

        # create lists of relevant surfaces for remeshing
        boundary_surfaces = []
        interface_surfaces = []

        # write stl input for volume meshing
        for key, tagid in label_to_id.items():
            surface_name = "-".join(key.lower().split())
            mask = cell_data["cell-tags"] == tagid
            surface = SurfaceMesh(surface_name, tris[mask, :], nodes=points)
            surface.write_to_stl(os.path.join(model.info.workdir, "part_" + surface_name + ".stl"))
            if "valve" in key or "septum" in key:
                interface_surfaces.append(surface)
            else:
                boundary_surfaces.append(surface)

        # add surface enclosing myocardium as surface
        part_surfaces: List[SurfaceMesh] = []
        for key, tagid in label_to_id_part.items():
            surface_name = "-".join(key.lower().split())
            mask = cell_data["tags"] == label_to_id_part[key]
            surface = SurfaceMesh(surface_name, tris[mask, :], nodes=points)
            part_surfaces.append(surface)
            surface.save(os.path.join(model.info.workdir, surface.name + ".vtk"))

        # Add these surfaces to the model
        model.mesh_raw.boundaries = boundary_surfaces
        model.mesh_raw.interfaces = interface_surfaces
        # Start remeshing
        model._remesh(post_mapping=False)

        # split epicardium boundary into 'real' epicardium and septum
        for b in model.mesh.boundaries:
            if "left-ventricle-epicardium" in b.name:
                triangles = b.faces.reshape(-1, 4)[:, 1:]
                region_ids = vtkmethods.get_connected_regions(b.nodes, triangles)
                unique_region_ids, counts = np.unique(region_ids, return_counts=True)
                sort_idx = np.argsort(counts)
                assert len(unique_region_ids) == 2, "Expecting two regions"
                unique_region_ids = unique_region_ids[sort_idx]
                mask_septum = region_ids == unique_region_ids[0]
                mask_epicardium = region_ids == unique_region_ids[1]
                model.mesh.boundaries.append(
                    SurfaceMesh(
                        "right-ventricle-septum", triangles=triangles[mask_septum, :], nodes=b.nodes
                    )
                )
                vertex = np.ones((len(triangles[mask_epicardium]), 1), dtype=int) * 3
                bb = np.hstack((vertex, triangles[mask_epicardium, :]))
                b.faces = bb.ravel()

        # write boundaries for debugging purposes
        for b in model.mesh.boundaries:
            b.save(os.path.join(model.info.workdir, b.name + ".stl"))

        # update lv and rv parts (assign tet ids to these parts). Use `enclosed by` surface
        # method for that
        # temp_vtk_file = os.path.join(model.info.workdir, "volume_mesh_postremesh.vtk")
        # model.mesh.write_to_vtk(temp_vtk_file)
        # vtk_volume = vtkmethods.read_vtk_unstructuredgrid_file(temp_vtk_file)
        # model.mesh.cell_data = {}
        # model.mesh.cell_data["tags"] = np.zeros(model.mesh.tetrahedrons.shape[0], dtype=float)
        # for surface_part in part_surfaces:
        #     # element_ids_inside = vtkmethods.cell_ids_inside_enclosed_surface(
        #     #     model.mesh, surface_part
        #     # )
        #     # print(len(element_ids_inside))
        #     x = model.mesh.cell_centers().select_enclosed_points(surface_part)['SelectedPoints']
        #     element_ids_inside = np.where(x == 1)[0]
        #     if "left" in surface_part.name:
        #         model.mesh.cell_data["tags"][element_ids_inside] = 1
        #     elif "right" in surface_part.name:
        #         model.mesh.cell_data["tags"][element_ids_inside] = 2

        import meshio

        tag1 = np.ones(len(points_1))
        tag2 = np.ones(len(points_2)) * 2
        vertex = np.linspace(0, len(points) - 1, num=len(points), dtype=int).reshape(-1, 1)
        tags = np.hstack((tag1, tag2))
        meshio.write_points_cells("x.vtk", points, [("vertex", vertex)], cell_data={"tags": [tags]})

        xx = vtkmethods.read_vtk_unstructuredgrid_file("x.vtk")
        model.mesh = vtkmethods.vtk_map_discrete_cell_data(xx, model.mesh, "tags")

        # extract septum
        model._extract_septum()
        for part in model.parts:
            if "Left" in part.name:
                part.tag_ids = [1]
            elif "Right" in part.name:
                part.tag_ids = [2]

        # call some additional methods to assign parts, surfaces, caps, etc
        model._assign_elements_to_parts()
        model._assign_surfaces_to_parts()
        model._assign_caps_to_parts()
        model._assign_cavities_to_parts()
        model._extract_apex()

        # add uvc coordinates to the model
        # 0. Compute longitudinal axis with:
        #   average cap centroid
        #   average apical point
        # 1. use longitudal axis to compute rotate to xy plane and compute uvc_l by z-coordinate
        lv_apex = model.left_ventricle.apex_points[1].xyz
        mv_centroid = [c.centroid for p in model.parts for c in p.caps if "mitral" in c.name][0]
        longitudinal_axis = lv_apex - mv_centroid
        # np.savetxt(os.path.join(model.info.workdir, "points.txt"),
        #            np.vstack([lv_apex,mv_centroid]),delimiter=",")
        from ansys.heart.preprocessor.mesh.geodisc import rodrigues_rot

        Prot = rodrigues_rot(model.mesh.nodes - lv_apex, longitudinal_axis, [0, 0, -1])
        Prot[:, 2] = Prot[:, 2] - np.min(Prot, axis=0)[2]
        scaling = Prot[:, 2] / np.max(Prot[:, 2])
        # model.mesh.point_data = {}
        model.mesh.point_data["uvc_longitudinal"] = scaling
        model.mesh.cell_data["fiber"] = np.tile(
            [0.0, 0.0, 1.0], (model.mesh.tetrahedrons.shape[0], 1)
        )
        model.mesh.cell_data["sheet"] = np.tile(
            [0.0, 1.0, 0.0], (model.mesh.tetrahedrons.shape[0], 1)
        )

        # rename cap
        for part in model.parts:
            for cap in part.caps:
                cap.name = cap.name.replace("left-ventricle-myocardium-", "")
                cap.name = cap.name.replace("right-ventricle-myocardium-", "")

        # write mesh
        model.mesh.write_to_vtk(os.path.join(model.info.workdir, "volume_mesh_postremesh.vtk"))

        model.dump_model()
        model.info.clean_workdir(["*.vtk", "*.trn", "*.log", "part_*.stl"])

    # write LS-DYNA files
    # Load model (e.g. when you skip the preprocessor):
    path_to_model = os.path.join(
        pathlib.Path(path_to_lv).parent, "BiVentricle", "heart_model.pickle"
    )
    model = models.HeartModel.load_model(path_to_model)
    model._add_nodal_areas()

    write_lsdyna_files = True
    if write_lsdyna_files:
        for writer in (
            writers.FiberGenerationDynaWriter(model),
            #     writers.ElectrophysiologyDynaWriter(model),
            # writers.MechanicsDynaWriter(model),
            # writers.ZeroPressureMechanicsDynaWriter(model),
            # writers.PurkinjeGenerationDynaWriter(model),
        ):
            exportdir = os.path.join(
                writer.model.info.workdir,
                writer.__class__.__name__.lower().replace("dynawriter", ""),
            )

            writer.model.mesh.write_to_vtk(
                os.path.join(writer.model.info.workdir, "volume_mesh.vtk")
            )
            writer.update()
            writer.export(exportdir)

    pass
