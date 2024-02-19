"""Module contains methods for interaction with Fluent meshing."""

import glob
import logging
import os
from pathlib import Path
import shutil
from typing import List, Union

LOGGER = logging.getLogger("pyheart_global.preprocessor")
# from importlib.resources import files

from ansys.heart.preprocessor._load_template import load_template

try:
    from ansys.heart.preprocessor.input import _InputBoundary, _InputModel
except ImportError:
    LOGGER.debug("Failed to import _InputBoundary, _InputModel")

import numpy as np

# from pkg_resources import resource_filename
import pkg_resources
import pyvista as pv

import ansys.heart.preprocessor.mesh.fluenthdf5 as hdf5  # noqa: F401
from ansys.heart.preprocessor.mesh.fluenthdf5 import FluentCellZone, FluentMesh

_template_directory = pkg_resources.resource_filename("ansys.heart.preprocessor", "templates")

_fluent_version = "22.2.0"

try:
    _show_fluent_gui = bool(int(os.environ["SHOW_FLUENT_GUI"]))
except KeyError:
    _show_fluent_gui = False
LOGGER.info(f"Showing Fluent gui: {_show_fluent_gui}")

# check whether containerized version of Fluent is used
if os.getenv("PYFLUENT_LAUNCH_CONTAINER"):
    _uses_container = True
else:
    _uses_container = False


try:
    import ansys.fluent.core as pyfluent
except ImportError:
    LOGGER.info(
        "Failed to import PyFluent. Considering installing "
        "pyfluent with `pip install ansys-fluent-core`."
    )


def _get_face_zones_with_filter(pyfluent_session, prefixes: list) -> list:
    """Get list of available boundaries in Fluent session that use any of the prefixes."""
    face_zones = []
    # get unique prefixes
    prefixes = list(set(prefixes))
    for prefix in prefixes:
        face_zones_with_prefix = pyfluent_session.scheme_eval.scheme_eval(
            f'(tgapi-util-convert-zone-ids-to-name-strings (get-face-zones-of-filter "{prefix}"))'
        )
        if face_zones_with_prefix:
            face_zones += face_zones_with_prefix
    # get only unique and avoid duplicates:
    face_zones = list(set(face_zones))
    return face_zones


def _organize_connected_regions(grid: pv.UnstructuredGrid, scalar: str = "part-id"):
    """Ensure that cells that belong to same part are connected."""
    LOGGER.debug("Re-organize connected regions.")
    part_ids = np.unique(grid.cell_data[scalar])
    grid.cell_data.set_scalars(np.arange(0, grid.n_cells, dtype=int), name="orig-cell-ids")

    # grid1 = grid.copy(deep=False)
    grid1 = grid.copy(deep=True)

    tets = grid.cells_dict[10]

    # get a list of orphan cells
    orphan_cell_ids = []
    for part_id in part_ids:
        mask = grid1.cell_data[scalar] == part_id
        grid2 = grid1.extract_cells(mask)

        conn = grid2.connectivity()

        if np.unique(conn.cell_data["RegionId"]).shape == (1,):
            continue
        # for each region, find to what "main" region it is connected.
        # use point
        for region in np.unique(conn.cell_data["RegionId"])[1:]:
            orphan_cell_ids = conn.cell_data["orig-cell-ids"][conn.cell_data["RegionId"] == region]
            point_ids = np.array([grid1.get_cell(id).point_ids for id in orphan_cell_ids]).flatten()

            mask = np.isin(tets, point_ids)
            connected_cell_id = np.argwhere(
                np.all(np.vstack([np.sum(mask, axis=1) > 1, np.sum(mask, axis=1) < 4]), axis=0)
            )[0]

            grid.cell_data["part-id"][orphan_cell_ids] = grid.cell_data["part-id"][
                connected_cell_id
            ]

        # orphan_cell_ids += list(conn.cell_data["orig-cell-ids"][conn.cell_data["RegionId"] > 0])

    # for each orphan cell

    return grid


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

    # NOTE: when using containerized version - we need to copy all the files
    # to and from the mounted volume given by pyfluent.EXAMPLES_PATH (default)
    if _uses_container:
        mounted_volume = pyfluent.EXAMPLES_PATH
        work_dir_meshing = os.path.join(mounted_volume, "tmp_meshing")
        num_cpus = 1
        show_gui = False
    else:
        work_dir_meshing = os.path.abspath(os.path.join(workdir, "meshing"))
        show_gui = _show_fluent_gui

    if os.path.isdir(work_dir_meshing):
        shutil.rmtree(work_dir_meshing)
    os.mkdir(work_dir_meshing)

    LOGGER.debug(f"Path to meshing directory: {work_dir_meshing}")

    path_to_output_old = path_to_output
    path_to_output = os.path.join(work_dir_meshing, "volume-mesh.msh.h5")

    min_size = mesh_size
    max_size = mesh_size
    growth_rate = 1.2

    # clean up any stls in the directory
    stls = glob.glob(os.path.join(work_dir_meshing, "*.stl"))
    for stl in stls:
        os.remove(stl)

    # write all boundaries
    model.write_part_boundaries(work_dir_meshing)

    session = pyfluent.launch_fluent(
        mode="meshing",
        precision="double",
        processor_count=2,
        start_transcript=False,
        show_gui=_show_fluent_gui,
        product_version=_fluent_version,
    )

    session.transcript.start(
        os.path.join(work_dir_meshing, "fluent_meshing.log"), write_to_stdout=False
    )

    # import files
    session.tui.file.import_.cad('no "' + work_dir_meshing + '" "*.stl" yes 40 yes mm')
    session.tui.objects.merge("'(*) heart")
    session.tui.objects.labels.create_label_per_zone("heart '(*)")
    session.tui.diagnostics.face_connectivity.fix_free_faces("objects '(*) merge-nodes yes 1e-3")
    session.tui.diagnostics.face_connectivity.fix_self_intersections(
        "objects '(heart) fix-self-intersection"
    )
    # smooth all zones
    face_zone_names = _get_face_zones_with_filter(session, "*")

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

    session.tui.file.write_mesh('"' + path_to_output + '"')
    # session.meshing.tui.file.read_journal(script)
    session.exit()

    if path_to_output != path_to_output_old:
        shutil.copy(path_to_output, path_to_output_old)

    path_to_output = path_to_output_old

    mesh = FluentMesh()
    mesh.load_mesh(path_to_output)

    # use part definitions to find which cell zone belongs to which part.
    for input_part in model.parts:
        surface = input_part.combined_boundaries

        if surface.is_manifold:
            check_surface = True
        else:
            check_surface = False
            LOGGER.warning(
                "Part {0} not manifold - disabled surface check.".format(input_part.name)
            )

        for cz in mesh.cell_zones:
            # use centroid of first cell to find which input part it belongs to.
            centroid = pv.PolyData(np.mean(mesh.nodes[cz.cells[0, :], :], axis=0))
            if np.all(
                centroid.select_enclosed_points(surface, check_surface=False).point_data[
                    "SelectedPoints"
                ]
            ):
                cz.id = input_part.id

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

    import ansys.fluent.core as pyfluent

    # NOTE: when using containerized version - we need to copy all the files
    # to and from the mounted volume given by pyfluent.EXAMPLES_PATH (default)
    if _uses_container:
        mounted_volume = pyfluent.EXAMPLES_PATH
        work_dir_meshing = os.path.join(mounted_volume, "tmp_meshing")
        num_cpus = 1
        show_gui = False
    else:
        work_dir_meshing = os.path.abspath(os.path.join(workdir, "meshing"))
        show_gui = _show_fluent_gui

    if os.path.isdir(work_dir_meshing):
        shutil.rmtree(work_dir_meshing)
    os.mkdir(work_dir_meshing)

    path_to_output_old = path_to_output
    path_to_output = os.path.join(work_dir_meshing, "volume-mesh.msh.h5")

    min_size = mesh_size
    max_size = mesh_size
    growth_rate = 1.2

    # clean up any stls in the directory
    stls = glob.glob(os.path.join(work_dir_meshing, "*.stl"))
    for stl in stls:
        os.remove(stl)

    # change boundary names to ensure max length is not exceeded
    old_to_new_boundary_names = {}
    for ii, b in enumerate(model.boundaries):
        if b.name in old_to_new_boundary_names.keys():
            b.name = old_to_new_boundary_names[b.name]
        else:
            if "interface" in b.name:
                tmp_name = "input_interface{:03d}".format(ii)
            else:
                tmp_name = "input_boundary{:03d}".format(ii)
            old_to_new_boundary_names[b.name] = tmp_name
            b.name = tmp_name

    new_to_old_boundary_names = {v: k for k, v in old_to_new_boundary_names.items()}

    # change part names to avoid issues such as characters that are now allowed.
    old_to_new_partnames = {}
    for ii, p in enumerate(model.parts):
        old_part_name = p.name
        new_part_name = "part_{:03d}".format(ii)
        old_to_new_partnames[old_part_name] = new_part_name
        p.name = new_part_name
    new_to_old_partnames = {v: k for k, v in old_to_new_partnames.items()}

    # # find interface names
    # interface_boundary_names = [
    #     new_name
    #     for old_name, new_name in old_to_new_boundary_names.items()
    #     if "interface-" in old_name
    # ]

    # write all boundaries
    model.write_part_boundaries(work_dir_meshing)

    # launch pyfluent
    session = pyfluent.launch_fluent(
        mode="meshing",
        precision="double",
        processor_count=4,
        start_transcript=False,
        show_gui=_show_fluent_gui,
        product_version=_fluent_version,
    )

    session.transcript.start(
        os.path.join(work_dir_meshing, "fluent_meshing.log"), write_to_stdout=False
    )

    # import stls
    session.tui.file.import_.cad('no "' + work_dir_meshing + '" "*.stl" yes 40 yes mm')

    # each stl is imported as a separate object. Wrap the different collections of stls to create
    # new surface meshes for each of the parts.
    session.tui.size_functions.set_global_controls(min_size, max_size, growth_rate)
    session.tui.scoped_sizing.compute('"yes"')

    session.tui.objects.extract_edges("'(*) feature 40")
    visited_parts = []
    used_boundary_names = []

    for part in model.parts:
        LOGGER.info("Wrapping " + part.name + "...")
        # wrap object.
        session.tui.objects.wrap.wrap(
            "'({0}) collectively {1} shrink-wrap external wrapped hybrid".format(
                " ".join(part.boundary_names), part.name
            )
        )
        # manage boundary names of wrapped surfaces
        face_zone_names = _get_face_zones_with_filter(session, ["input_*:*"])
        if not face_zone_names:
            continue

        for face_zone_name in face_zone_names:
            old_name = face_zone_name

            # if old name contains part name before delimiter : -> skip this boundary
            prefix = old_name.split(":")[0]
            if any([prefix == p.lower() for p in visited_parts]):
                LOGGER.debug(f"Skip {face_zone_name}")
                continue

            new_name = part.name + ":" + prefix

            if new_name in used_boundary_names:
                LOGGER.debug(f"Name of boundary {new_name} already used.")
                continue

            session.tui.boundary.manage.name(old_name + " " + new_name)

            used_boundary_names += [new_name]

        visited_parts += [part.name]

    # print(_get_face_zones_with_filter(session, ["*"]))

    LOGGER.info("Wrapping model...")

    # wrap entire model in one pass so that we can create a single volume mesh. Use list of all
    # input boundaries are given as input. External material point for meshing.
    session.tui.objects.wrap.wrap(
        "'({0}) collectively {1} shrink-wrap external wrapped hybrid".format(
            " ".join(model.boundary_names), "model"
        )
    )

    # print(_get_face_zones_with_filter(session, ["*"]))

    # rename boundaries accordingly.
    prefixes = [bn + ":*" for bn in model.boundary_names]
    face_zone_names = _get_face_zones_with_filter(session, prefixes)

    if not face_zone_names:
        LOGGER.error("Expecting face zones to rename.")

    used_names = []
    for face_zone_name in face_zone_names:
        # Exclude renaming of boundaries that include name of visited parts
        old_name = face_zone_name

        new_name = "model" + ":" + old_name.split(":")[0]

        # find unique name
        rename_success = False
        ii = 0
        while not rename_success:
            if new_name not in used_names:
                break
            else:
                new_name = new_name = "model" + ":" + old_name.split(":")[0] + "_{:03d}".format(ii)
            ii += 1

        session.tui.boundary.manage.name(old_name + " " + new_name)

        used_names += [new_name]

    # mesh the entire model in one go.
    session.tui.objects.volumetric_regions.compute("model")
    session.tui.mesh.auto_mesh("model yes pyramids tet no")

    # clean up geometry objects
    session.tui.objects.delete_all_geom()

    # write mesh
    if os.path.isfile(path_to_output):
        os.remove(path_to_output)

    session.tui.file.write_mesh('"' + path_to_output + '"')
    session.exit()

    shutil.copy(path_to_output, path_to_output_old)

    path_to_output = path_to_output_old

    # Update the cell zones such that for each part we have a separate cell zone.
    mesh = FluentMesh()
    mesh.load_mesh(path_to_output)

    num_cells = mesh.cell_zones[0].cells.shape[0]

    # convert to unstructured grid.
    # cells = np.hstack([np.ones((num_cells, 1), dtype=int) * 4, mesh.cell_zones[0].cells])
    # celltypes = [pv.CellType.TETRA] * num_cells
    # grid = pv.UnstructuredGrid(cells.flatten(), celltypes, mesh.nodes)
    grid = mesh._to_vtk()

    # represent cell centroids as point cloud assign part-ids to cells.
    cell_centroids = grid.cell_centers()
    cell_centroids.point_data.set_scalars(name="part-id", scalars=0)

    # NOTE: should use wrapped surfaces to select part.
    # assign wrapped boundaries to input parts.
    for ii, part in enumerate(model.parts):
        face_zones_wrapped = [fz for fz in mesh.face_zones if part.name in fz.name]
        if len(face_zones_wrapped) == 0:
            LOGGER.error(f"Did not find any wrapped face zones for {part.name}")

        # replace with remeshed counterpart.
        for jj, boundary in enumerate(part.boundaries):
            for fz in face_zones_wrapped:
                if boundary.name in fz.name:
                    break
            remeshed_boundary = _InputBoundary(
                mesh.nodes,
                faces=np.hstack([np.ones(fz.faces.shape[0], dtype=int)[:, None] * 3, fz.faces]),
                id=boundary.id,
            )
            model.parts[ii].boundaries[jj] = remeshed_boundary

    # use individual wrapped parts to identify the parts of the wrapped model.
    for part in model.parts:
        if not part.is_manifold:
            LOGGER.warning("Part is not manifold.")

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
    try:
        cell_centroids_2.point_data.remove("orig_indices")
    except:
        KeyError
    try:
        cell_centroids_2.point_data.remove("cell-zone-ids")
    except:
        KeyError
    try:
        cell_centroids_2.cell_data.remove("cell-zone-ids")
    except:
        KeyError

    cell_centroids_1.point_data.remove("part-id")
    cell_centroids_1.cell_data.remove("cell-zone-ids")

    cell_centroids_1 = cell_centroids_1.interpolate(
        cell_centroids_2, n_points=1, pass_cell_data=False
    )
    cell_centroids.point_data["part-id"][orig_indices_1] = cell_centroids_1.point_data["part-id"]

    # assign part-ids to grid
    grid.cell_data["part-id"] = np.array(cell_centroids.point_data["part-id"], dtype=int)

    # Ensure that parts are continuous and well connected.
    grid = _organize_connected_regions(grid, scalar="part-id")

    if np.any(grid.cell_data["part-id"] == 0):
        raise ValueError("Invalid mesh, not all elements assigned to a part.")

    # change FluentMesh object accordingly.
    idx_sorted = np.argsort(np.array(grid.cell_data["part-id"], dtype=int))
    partids_sorted = np.sort(np.array(grid.cell_data["part-id"], dtype=int))

    new_mesh = mesh
    new_mesh.cells = new_mesh.cells[idx_sorted]
    new_mesh.cell_zones: List[FluentCellZone] = []

    for part in model.parts:
        part.name = new_to_old_partnames[part.name]
        cell_zone = FluentCellZone(
            min_id=np.argwhere(partids_sorted == part.id)[0][0],
            max_id=np.argwhere(partids_sorted == part.id)[-1][0],
            name=part.name,
            cid=part.id,
        )
        cell_zone.get_cells(new_mesh.cells)
        new_mesh.cell_zones.append(cell_zone)

    # remove any unused face zones.
    new_mesh.face_zones = [fz for fz in new_mesh.face_zones if "part" not in fz.name.lower()]

    # rename face zones - rename to original input names.
    for fz in new_mesh.face_zones:
        if "interior" in fz.name:
            continue
        fz.name = fz.name.replace("model:", "")
        if ":" in fz.name:
            fz.name = fz.name.split(":")[0]
        try:
            fz.name = new_to_old_boundary_names[fz.name]
        except KeyError:
            LOGGER.debug(f"Failed to rename {fz.name}")

    new_grid = new_mesh._to_vtk()
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
        path_to_blood_pool_script = os.path.join(
            _template_directory, "fluent_meshing_add_blood_mesh_template.jou"
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
        show_gui = _show_fluent_gui

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
        start_transcript=False,
        show_gui=show_gui,
        product_version=_fluent_version,
    )
    if session.health_check_service.status() != "SERVING":
        LOGGER.error("Fluent session failed. Exiting Fluent")
        session.exit()
        exit()

    min_size = mesh_size
    max_size = mesh_size
    growth_rate_wrap = 1.2

    session.transcript.start(
        os.path.join(work_dir_meshing, "fluent_meshing.log"), write_to_stdout=False
    )

    # import files
    session.tui.file.import_.cad("no " + work_dir_meshing + " part_*.stl yes 40 yes mm")
    # session.transcript.start
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


if __name__ == "__main__":
    LOGGER.info("Protected")
