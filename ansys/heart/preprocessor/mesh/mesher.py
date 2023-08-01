"""Module contains methods for interaction with Fluent meshing."""
import glob
import os
from pathlib import Path
import shutil
import subprocess
from typing import List, Union

from ansys.heart.custom_logging import LOGGER
from ansys.heart.preprocessor._load_template import load_template
from ansys.heart.preprocessor.input import _InputModel
import ansys.heart.preprocessor.mesh.fluenthdf5 as hdf5  # noqa: F401
from ansys.heart.preprocessor.mesh.fluenthdf5 import FluentCellZone, FluentMesh
import numpy as np

# from pkg_resources import resource_filename
import pkg_resources
import pyvista as pv

_template_directory = pkg_resources.resource_filename("ansys.heart.preprocessor", "templates")

_fluent_version = "23.1.0"

try:
    import ansys.fluent.core as pyfluent
except ImportError:
    LOGGER.info(
        "Failed to import PyFluent. Considering installing "
        "pyfluent with `pip install ansys-fluent-core`."
    )


def mesh_from_manifold_input_model(
    model: _InputModel,
    workdir: Union[str, Path],
    path_to_output: Union[str, Path],
    mesh_size: float = 2.0,
) -> FluentMesh:
    """Create mesh from good-quality manifold input model.

    Parameters
    ----------
    model : _InputModel
        Input model.
    workdir : Union[str, Path]
        Working directory.
    path_to_output : Union[str, Path]
        Path to the resulting Fluent mesh file.
    mesh_size : float, optional
        Uniform mesh size to use for both wrapping and filling the volume, by default 2.0

    Returns
    -------
    FluentMesh
        The volume mesh with cell and face zones.
    """
    if not isinstance(model, _InputModel):
        raise ValueError(f"Expecting input to be of type {str(_InputModel)}")

    if not os.path.isdir(workdir):
        os.makedirs(workdir)

    min_size = mesh_size
    max_size = mesh_size
    growth_rate = 1.2

    # clean up any stls in the directory
    stls = glob.glob(os.path.join(workdir, "*.stl"))
    for stl in stls:
        os.remove(stl)

    # write all boundaries
    model.write_part_boundaries(workdir)

    import ansys.fluent.core as pyfluent

    session = pyfluent.launch_fluent(
        mode="meshing",
        precision="double",
        processor_count=2,
        start_transcript=True,
        show_gui=True,
        product_version=_fluent_version,
    )

    # import files
    session.tui.file.import_.cad("no " + workdir + " *.stl yes 40 yes mm")
    session.tui.file.start_transcript(os.path.join(workdir, "fluent_meshing.log"))
    session.tui.objects.merge("'(*) heart")
    session.tui.objects.labels.create_label_per_zone("heart '(*)")
    session.tui.diagnostics.face_connectivity.fix_free_faces("objects '(*) merge-nodes yes 1e-3")
    session.tui.diagnostics.face_connectivity.fix_self_intersections(
        "objects '(heart) fix-self-intersection"
    )
    # smooth all zones
    face_zone_names = session.scheme_eval.scheme_eval(
        '(tgapi-util-convert-zone-ids-to-name-strings (get-face-zones-of-filter "*"))'
    )
    for fz in face_zone_names:
        session.tui.boundary.modify.select_zone(fz)
        session.tui.boundary.modify.smooth()

    session.tui.objects.create_intersection_loops("collectively '(*)")
    session.tui.boundary.feature.create_edge_zones("(*) fixed-angle 70 yes")
    # create size field
    session.tui.size_functions.set_global_controls(min_size, max_size, growth_rate)
    session.tui.scoped_sizing.compute("yes")

    # remesh surface
    session.tui.boundary.remesh.remesh_face_zones_conformally("'(*) '(*) 40 20 yes")

    # some diagnostics
    session.tui.diagnostics.face_connectivity.fix_self_intersections(
        "objects '(heart) fix-self-intersection"
    )
    session.tui.diagnostics.face_connectivity.fix_duplicate_faces("objects '(heart)")

    # convert to mesh object
    session.tui.objects.change_object_type("'(heart) mesh y")

    # compute volumes
    session.tui.objects.volumetric_regions.compute("heart", "no")

    # start auto meshing
    session.tui.mesh.tet.controls.cell_sizing("size-field")
    session.tui.mesh.auto_mesh("heart", "yes", "pyramids", "tet", "no")
    session.tui.mesh.modify.auto_node_move("(*)", "(*)", 0.3, 50, 120, "yes", 5)
    session.tui.objects.delete_all_geom()
    session.tui.mesh.zone_names_clean_up()
    # session.tui.mesh.check_mesh()
    # session.tui.mesh.check_quality()
    session.tui.boundary.manage.remove_suffix("(*)")

    session.tui.mesh.prepare_for_solve("yes")

    # write to file
    session.tui.file.write_mesh(path_to_output)
    # session.meshing.tui.file.read_journal(script)
    session.exit()

    mesh = FluentMesh()
    mesh.load_mesh(path_to_output)

    return mesh


def mesh_from_non_manifold_input_model(
    model: _InputModel,
    workdir: Union[str, Path],
    path_to_output: Union[str, Path],
    mesh_size: float = 2.0,
) -> FluentMesh:
    """Generate mesh from non-manifold poor quality input model.

    Parameters
    ----------
    model : _InputModel
        Input model.
    workdir : Union[str, Path]
        Working directory.
    path_to_output : Union[str, Path]
        Path to the resulting Fluent mesh file.
    mesh_size : float, optional
        Uniform mesh size to use for both wrapping and filling the volume, by default 2.0

    Notes
    -----
    Uses Fluent wrapping technology to wrap the individual parts first to create manifold
    parts. Consequently wrap the entire model and use the manifold parts to split the
    wrapped model into the different cell zones.

    Returns
    -------
    FluentMesh
        The volume mesh with cell and face zones.
    """
    if not isinstance(model, _InputModel):
        raise ValueError(f"Expecting input to be of type {str(_InputModel)}")

    if not os.path.isdir(workdir):
        os.makedirs(workdir)

    min_size = mesh_size
    max_size = mesh_size
    growth_rate = 1.2

    # clean up any stls in the directory
    stls = glob.glob(os.path.join(workdir, "*.stl"))
    for stl in stls:
        os.remove(stl)

    # change boundary names to ensure max length is not exceeded
    boundary_name_map_old_to_new = {}
    for ii, b in enumerate(model.boundaries):
        if b.name in boundary_name_map_old_to_new.keys():
            b.name = boundary_name_map_old_to_new[b.name]
        else:
            if "interface" in b.name:
                tmp_name = "interface{:03d}".format(ii)
            else:
                tmp_name = "boundary{:03d}".format(ii)
            boundary_name_map_old_to_new[b.name] = tmp_name
            b.name = tmp_name

    # find interface names
    interface_boundary_names = [
        new_name
        for old_name, new_name in boundary_name_map_old_to_new.items()
        if "interface-" in old_name
    ]

    # write all boundaries
    model.write_part_boundaries(workdir)

    # launch pyfluent
    session = pyfluent.launch_fluent(
        mode="meshing",
        precision="double",
        processor_count=2,
        start_transcript=True,
        show_gui=True,
        product_version=_fluent_version,
    )

    # import stls
    session.tui.file.import_.cad("no " + workdir + " *.stl yes 40 yes mm")
    session.tui.file.start_transcript(os.path.join(workdir, "fluent_meshing.log"))

    # each stl is imported as a separate object. Wrap the different collections of stls to create
    # new surface meshes for each of the parts.
    session.tui.size_functions.set_global_controls(min_size, max_size, growth_rate)
    session.tui.scoped_sizing.compute("yes")

    session.tui.objects.extract_edges("'(*) feature 40")
    visited_parts = []
    for part in model.parts:
        LOGGER.info("Wrapping " + part.name + "...")
        # wrap object.
        session.tui.objects.wrap.wrap(
            "'({0}) collectively {1} shrink-wrap external wrapped hybrid".format(
                " ".join(part.boundary_names), part.name
            )
        )
        # manage boundary names of wrapped surfaces
        face_zone_names = session.scheme_eval.scheme_eval(
            '(tgapi-util-convert-zone-ids-to-name-strings (get-face-zones-of-filter "*:*"))'
        )
        for face_zone_name in face_zone_names:
            # Exclude renaming of boundaries that include name of visited parts
            if any([s + ":" in face_zone_name for s in visited_parts]):
                continue

            old_name = face_zone_name
            new_name = part.name + ":" + old_name.split(":")[0]
            session.tui.boundary.manage.name(old_name + " " + new_name)

        visited_parts += [part.name]

    LOGGER.info("Wrapping model...")
    # create a usable material point which is inside the model that can be used for shrink-wrapping
    material_point = [52, 152, 369.93]
    session.tui.material_point.create_material_point(
        "myocardium {:f} {:f} {:f}".format(*material_point)
    )
    resolution_factor = 0.25

    # wrap entire model in one pass so that we can create a single volume mesh.
    boundaries_to_use_for_wrapping = session.scheme_eval.scheme_eval(
        '(tgapi-util-convert-zone-ids-to-name-strings (get-face-zones-of-filter "boundary*"))'
    )

    session.tui.objects.wrap.wrap(
        "'({0}) collectively {1} shrink-wrap myocardium hybrid {2}".format(
            " ".join(boundaries_to_use_for_wrapping), "model", resolution_factor
        )
    )

    # with external material point
    # session.tui.objects.wrap.wrap(
    #     "'({0}) collectively {1} shrink-wrap external wrapped hybrid".format(
    #         " ".join(boundaries_to_use_for_wrapping), "model"
    #     )
    # )

    # rename boundaries accordingly.
    face_zone_names = session.scheme_eval.scheme_eval(
        '(tgapi-util-convert-zone-ids-to-name-strings (get-face-zones-of-filter "*:*"))'
    )
    for face_zone_name in face_zone_names:
        # Exclude renaming of boundaries that include name of visited parts
        if any([s + ":" in face_zone_name for s in visited_parts]):
            continue

        old_name = face_zone_name
        new_name = "model" + ":" + old_name.split(":")[0]
        session.tui.boundary.manage.name(old_name + " " + new_name)

    # mesh the entire model in one go.
    session.tui.objects.volumetric_regions.compute("model")
    session.tui.mesh.auto_mesh("model")

    # clean up geometry objects
    session.tui.objects.delete_all_geom()

    # write mesh
    session.tui.file.write_mesh(path_to_output)
    session.exit()

    # Update the cell zones such that for each part we have a separate cell zone.
    mesh = FluentMesh()
    mesh.load_mesh(path_to_output)

    num_cells = mesh.cell_zones[0].cells.shape[0]

    # convert to unstructured grid.
    cells = np.hstack([np.ones((num_cells, 1), dtype=int) * 4, mesh.cell_zones[0].cells])
    celltypes = [pv.CellType.TETRA] * num_cells
    grid = pv.UnstructuredGrid(cells.flatten(), celltypes, mesh.nodes)

    # represent cell centroids as point cloud assign part-ids to cells.
    cell_centroids = grid.cell_centers()
    cell_centroids.point_data.set_scalars(name="part-id", scalars=0)

    for part in model.parts:
        LOGGER.warning("Disabled check for manifold surface before computing the enclosed points.")
        cell_centroids = cell_centroids.select_enclosed_points(
            part.combined_boundaries, check_surface=False
        )
        cell_centroids.point_data["part-id"][
            cell_centroids.point_data["SelectedPoints"] == 1
        ] = part.id

    # Use closest-point interpolation to assign part-ids to cell centers that are
    # not enclosed by any of the wrapped parts
    cell_centroids["orig_indices"] = np.arange(cell_centroids.n_points, dtype=np.int32)
    cell_centroids.point_data.remove("SelectedPoints")
    cell_centroids_1 = cell_centroids.remove_cells(
        cell_centroids.point_data["part-id"] != 0, inplace=False
    )
    orig_indices_1 = cell_centroids_1.point_data["orig_indices"]
    cell_centroids_2 = cell_centroids.remove_cells(
        cell_centroids.point_data["part-id"] == 0, inplace=False
    )

    cell_centroids_1.point_data.remove("part-id")
    cell_centroids_1 = cell_centroids_1.interpolate(
        cell_centroids_2, n_points=1, pass_cell_data=False
    )
    cell_centroids.point_data["part-id"][orig_indices_1] = cell_centroids_1.point_data["part-id"]

    # assign part-ids to grid
    grid.cell_data.set_scalars(scalars=cell_centroids.point_data["part-id"], name="part-id")

    # change FluentMesh object accordingly.
    idx_sorted = np.argsort(np.array(grid.cell_data["part-id"], dtype=int))
    partids_sorted = np.sort(np.array(grid.cell_data["part-id"], dtype=int))

    new_mesh = mesh
    new_mesh.cells = new_mesh.cells[idx_sorted]
    new_mesh.cell_zones: List[FluentCellZone] = []

    for part in model.parts:
        cell_zone = FluentCellZone(
            min_id=np.where(partids_sorted == part.id)[0][0],
            max_id=np.where(partids_sorted == part.id)[0][-1],
            name=part.name,
            cid=part.id,
        )
        cell_zone.get_cells(new_mesh.cells)
        new_mesh.cell_zones.append(cell_zone)

    # remove any unused face zones.
    new_mesh.face_zones = [fz for fz in new_mesh.face_zones if "part" not in fz.name.lower()]

    # rename face zones
    for fz in new_mesh.face_zones:
        fz.name = fz.name.replace("model:", "")

    return new_mesh


def mesh_heart_model_by_fluent(
    path_to_stl_directory: str,
    path_to_output: str,
    mesh_size: float = 2.0,
    add_blood_pool: bool = False,
    show_gui: bool = False,
):
    """
    Use Fluent meshing to wrap the surface and create tetrahedral volume mesh.

    Note
    ----
    Optionally extracts the blood pool.
    """
    import ansys.fluent.core as pyfluent

    # make sure we are using absolute path
    path_to_stl_directory = os.path.abspath(path_to_stl_directory)
    path_to_output = os.path.abspath(path_to_output)

    # change directory to directory of stl file
    old_directory = os.getcwd()
    working_directory = path_to_stl_directory
    os.chdir(working_directory)

    if add_blood_pool:
        path_to_blood_pool_script = pkg_resources.resource_filename(
            "ansys.heart.preprocessor", "templates/fluent_meshing_add_blood_mesh_template.jou"
        )
        f = open(path_to_blood_pool_script, "r")
        blood_pool_script = "".join(f.readlines())
        f.close()
    else:
        blood_pool_script = ""

    var_for_template = {
        "work_directory": working_directory,
        "output_path": path_to_output,
        "mesh_size": mesh_size,
        "blood_pool_script": blood_pool_script,
    }

    template = load_template("fluent_meshing_template_improved_2.jou")

    script = os.path.join(path_to_stl_directory, "fluent_meshing.jou")

    with open(script, "w") as f:
        f.write(template.render(var_for_template))

    num_cpus = 2

    # check whether containerized version of Fluent is used
    if os.getenv("PYFLUENT_LAUNCH_CONTAINER"):
        LOGGER.debug("Launching Fluent as container...")
        uses_container = True
    else:
        uses_container = False

    # NOTE: when using containerized version - we need to copy all the files
    # to and from the mounted volume given by pyfluent.EXAMPLES_PATH (default)
    if uses_container:
        mounted_volume = pyfluent.EXAMPLES_PATH
        work_dir_meshing = os.path.join(mounted_volume, "tmp_meshing")
        num_cpus = 1
        show_gui = False
    else:
        work_dir_meshing = os.path.abspath(os.path.join(working_directory, "meshing"))

    if os.path.isdir(work_dir_meshing):
        shutil.rmtree(work_dir_meshing)
    os.mkdir(work_dir_meshing)

    path_to_output_old = path_to_output
    path_to_output = os.path.join(work_dir_meshing, "volume-mesh.msh.h5")

    # copy all necessary files to meshing directory
    files_to_copy = glob.glob("part*.stl") + glob.glob("fluent_meshing.jou")
    for file in files_to_copy:
        shutil.copyfile(file, os.path.join(work_dir_meshing, file))

    LOGGER.debug("Starting meshing in directory: {}".format(work_dir_meshing))
    # start fluent session
    session = pyfluent.launch_fluent(
        mode="meshing",
        precision="double",
        processor_count=num_cpus,
        start_transcript=True,
        show_gui=show_gui,
        product_version=_fluent_version,
    )
    if session.check_health() != "SERVING":
        LOGGER.error("Fluent session failed. Exiting Fluent")
        session.stop_transcript()
        session.exit()
        exit()

    session.start_transcript()

    min_size = mesh_size
    max_size = mesh_size
    growth_rate_wrap = 1.2

    # import files
    session.tui.file.import_.cad("no " + work_dir_meshing + " part_*.stl yes 40 yes mm")
    session.tui.file.start_transcript(work_dir_meshing, "fluent_meshing.log")
    session.tui.objects.merge("'(*) heart")
    session.tui.objects.labels.create_label_per_zone("heart '(*)")
    session.tui.diagnostics.face_connectivity.fix_free_faces("objects '(*) merge-nodes yes 1e-3")
    session.tui.objects.create_intersection_loops("collectively '(*)")
    session.tui.boundary.feature.create_edge_zones("(*) fixed-angle 70 yes")

    # set up size field for wrapping
    session.tui.size_functions.set_global_controls(min_size, max_size, growth_rate_wrap)
    session.tui.scoped_sizing.compute("yes")

    # wrap objects
    session.tui.objects.wrap.wrap(
        "'(*)",
        "collectively",
        "wrapped-myocardium",
        "shrink-wrap",
        "external",
        "wrap",
        "hybrid",
        0.8,
    )

    if add_blood_pool:
        # if adding blood pool:
        # ; ------------------------------------------------------------------
        # ; script for auto-generating caps for all endocardial parts
        # ; and extracting blood pool volume
        # ; This first copies all the endocardial zones and septum and
        # ; and closes the the cavities based on the free faces
        # ; Note that the auto-patch utility may not work in all cases.
        # ; ------------------------------------------------------------------
        session.tui.objects.delete_all_geom()
        session.scheme_eval.scheme_eval(
            "(define zone-ids-endo (get-face-zones-of-filter '*endocardium*) )"
        )
        session.scheme_eval.scheme_eval(
            "(define zone-id-septum (get-face-zones-of-filter '*septum*) )"
        )
        a = session.scheme_eval.scheme_eval(
            "(define zone-ids-to-copy (append zone-ids-endo zone-id-septum) )"
        )
        session.scheme_eval.scheme_eval(
            '(ti-menu-load-string  ( format #f "/boundary/manage/copy ~a" zone-ids-to-copy) )'
        )
        session.scheme_eval.scheme_eval(
            "(define zone-ids-to-patch (get-unreferenced-face-zones-of-filter '*) )"
        )
        session.scheme_eval.scheme_eval(
            '(ti-menu-load-string  ( format #f "/boundary/modify/auto-patch-holes ~a"'
            "zone-ids-to-patch) )"
        )
        session.scheme_eval.scheme_eval(
            "(define zone-ids-patch (get-unreferenced-face-zones-of-filter '*patch*) )"
        )
        session.tui.objects.create("valves fluid 3 '(*-patch-*) '() mesh yes")

        # to here:
        session.tui.objects.delete_unreferenced_faces_and_edges()
        session.tui.objects.labels.create_label_per_zone("valves '(*)")
        session.tui.objects.labels.remove_zones("valves valves '(*)")
        session.tui.objects.merge("'(wrapped-myocardium valves) wrapped-myocardium")
        session.tui.diagnostics.face_connectivity.fix_duplicate_faces(
            "objects '(wrapped-myocardium)"
        )
        session.tui.boundary.remesh.controls.intersect.remesh_post_intersection("no")
        session.tui.boundary.remesh.intersect_all_face_zones("yes 40 0.05 no")
        session.tui.boundary.manage.remove_suffix("'(*)")
        session.tui.diagnostics.face_connectivity.fix_duplicate_faces(
            "objects '(wrapped-myocardium)"
        )

    # compute volumetric regions
    session.tui.objects.volumetric_regions.compute("wrapped-myocardium", "no")
    session.tui.objects.volumetric_regions.change_type("Wrapped-myocardium", "'(*)", "fluid")
    session.tui.objects.volumetric_regions.change_type("wrapped-myocardium", "(heart)", "solid")

    # start auto meshing
    session.tui.mesh.tet.controls.cell_sizing("size-field")
    session.tui.mesh.auto_mesh("wrapped-myocardium", "yes", "pyramids", "tet", "no")
    session.tui.mesh.modify.auto_node_move("(*)", "(*)", 0.3, 50, 120, "yes", 5)
    session.tui.objects.delete_all_geom()
    session.tui.mesh.zone_names_clean_up()
    session.tui.mesh.check_mesh()
    session.tui.mesh.check_quality()
    session.tui.boundary.manage.remove_suffix("(*)")

    # prepare for solve (removes unused nodes, faces, etc.)
    session.tui.mesh.prepare_for_solve("yes")

    # write to file
    session.tui.file.write_mesh(path_to_output)
    # session.meshing.tui.file.read_journal(script)
    session.exit()

    shutil.copy(path_to_output, path_to_output_old)

    # shutil.rmtree(work_dir_meshing)

    # change back to old directory
    os.chdir(old_directory)

    return


def _shrink_by_spaceclaim(input, output):
    """Use SpaceClaim shrinkwrapping to wrap surface and create high quality surface mesh."""
    try:
        from ansys.heart.preprocessor import SC_EXE
    except ImportError:
        LOGGER.error("Failed to import space claim path")
        return

    var_for_template = {
        "input": input,
        "output": output,
    }
    posixpath_template = str(os.path.join(_template_directory, "spaceclaim_shrink"))
    template = load_template("spaceclaim_shrink")

    script = "spaceclaim_shrink.py"

    with open(script, "w") as f:
        f.write(template.render(var_for_template))

    options = ["/Headless=False", "/Splash=False", "/ExitAfterScript=True"]
    subprocess.call([SC_EXE, "/RunScript=" + script, *options])

    return


def _run_gmsh(infile: str, outfile: str, mesh_size):
    """Run GMESH with specified in/output file and target mesh size.

    Arguments
    ---------
        infile (str): Path to .stl input file
        outfile (str): path to .vtk output
    """
    try:
        import gmsh
    except:
        ImportError("GMESH not installed. Install through pip install gmesh")

    gmsh.initialize()

    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size)
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)

    # load STL file
    gmsh.merge(infile)

    # split input surface mesh based on an angle threshold of 40 degrees between
    # triangles, and generate new discrete surfaces suitable for reparametrization
    gmsh.model.mesh.classifySurfaces(90 * np.pi / 180.0, True, True)

    # create a geometry (through reparametrization) for all discrete curves and
    # discrete surfaces
    gmsh.model.mesh.createGeometry()

    # add a volume
    s = gmsh.model.getEntities(2)
    l = gmsh.model.geo.addSurfaceLoop([s[i][1] for i in range(len(s))])
    gmsh.model.geo.addVolume([l])

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    gmsh.write(outfile)
    gmsh.finalize()
    # if '-nopopup' not in sys.argv:
    #     gmsh.fltk.run()
    #
    # gmsh.finalize()
    return


if __name__ == "__main__":
    LOGGER.info("Protected")
