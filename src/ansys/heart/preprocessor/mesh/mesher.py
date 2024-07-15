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

"""Module contains methods for interaction with Fluent meshing."""

import glob
import os
from pathlib import Path
import shutil
from typing import List, Union

import numpy as np
import pyvista as pv

from ansys.heart.core import LOG as LOGGER
from ansys.heart.preprocessor.input import _InputBoundary, _InputModel
import ansys.heart.preprocessor.mesh.fluenthdf5 as hdf5  # noqa: F401
from ansys.heart.preprocessor.mesh.fluenthdf5 import FluentCellZone, FluentMesh
from ansys.heart.preprocessor.mesh.objects import SurfaceMesh
from ansys.heart.preprocessor.mesh.vtkmethods import add_solid_name_to_stl

# os.environ["SHOW_FLUENT_GUI"] = "1"

_fluent_version = "24.1.0"
_show_fluent_gui: bool = False
_uses_container: bool = True

try:
    _show_fluent_gui = bool(int(os.environ["SHOW_FLUENT_GUI"]))
except KeyError:
    _show_fluent_gui = False
LOGGER.info(f"Showing Fluent gui: {_show_fluent_gui}")

# check whether containerized version of Fluent is used
if os.getenv("PYFLUENT_LAUNCH_CONTAINER"):
    _uses_container = True
    _show_fluent_gui = False
else:
    _uses_container = False

if _show_fluent_gui:
    _fluent_ui_mode = "gui"
else:
    _fluent_ui_mode = "hidden_gui"

try:
    import ansys.fluent.core as pyfluent
    from ansys.fluent.core.session_meshing import Meshing as MeshingSession
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
            connected_cell_ids = np.argwhere(
                np.all(np.vstack([np.sum(mask, axis=1) > 1, np.sum(mask, axis=1) < 4]), axis=0)
            ).flatten()
            unique_ids, counts = np.unique(
                grid.cell_data["part-id"][connected_cell_ids], return_counts=True
            )
            if unique_ids.shape[0] > 1:
                LOGGER.debug("More than 1 candidate.")

            grid.cell_data["part-id"][orphan_cell_ids] = unique_ids[np.argmax(counts)]

        # orphan_cell_ids += list(conn.cell_data["orig-cell-ids"][conn.cell_data["RegionId"] > 0])

    # for each orphan cell

    return grid


def _get_fluent_meshing_session() -> MeshingSession:
    """Get a Fluent Meshing session."""
    # NOTE: when using containerized version - we need to copy all the files
    # to and from the mounted volume given by pyfluent.EXAMPLES_PATH (default)
    if _uses_container:
        num_cpus = 1
    else:
        num_cpus = 2

    session = pyfluent.launch_fluent(
        mode="meshing",
        precision="double",
        processor_count=num_cpus,
        start_transcript=False,
        ui_mode=_fluent_ui_mode,
        product_version=_fluent_version,
        start_container=_uses_container,
    )

    return session


def _wrap_part(session: MeshingSession, boundary_names: list, wrapped_part_name: str) -> list:
    """Invoke the wrapper to wrap a part based on a list of boundary names."""
    pre_wrap_facezones = _get_face_zones_with_filter(session, ["*"])
    session.tui.objects.wrap.wrap(
        "'({0}) collectively {1} shrink-wrap external wrapped hybrid".format(
            " ".join(boundary_names), wrapped_part_name
        )
    )
    post_wrap_facezones = _get_face_zones_with_filter(session, ["*"])
    wrapped_face_zones = list(set(post_wrap_facezones) - set(pre_wrap_facezones))

    # rename the "new" face zones accordingly:
    wrapped_face_zone_names = []
    for face_zone_name in wrapped_face_zones:
        # Exclude renaming of boundaries that include name of visited parts
        old_name = face_zone_name
        new_name = wrapped_part_name + ":" + old_name.split(":")[0]
        # find unique name
        rename_success = False
        ii = 0
        while not rename_success:
            if new_name not in wrapped_face_zone_names:
                break
            else:
                new_name = new_name = (
                    wrapped_part_name + ":" + old_name.split(":")[0] + "_{:03d}".format(ii)
                )
            ii += 1
        session.tui.boundary.manage.name(old_name + " " + new_name)
        wrapped_face_zone_names += [new_name]

    return wrapped_face_zone_names


def mesh_fluid_cavities(
    fluid_boundaries: List[SurfaceMesh],
    caps: List[SurfaceMesh],
    workdir: str,
    remesh_caps: bool = True,
) -> FluentMesh:
    """Mesh the fluid cavities.

    Parameters
    ----------
    fluid_boundaries : List[SurfaceMesh]
        List of fluid boundaries used for meshing.
    caps : List[SurfaceMesh]
        List of caps that close each of the cavities.
    workdir : str
        Working directory
    remesh_caps : bool, optional
        Flag indicating whether to remesh the caps, by default True

    Returns
    -------
    Path
        Path to the .msh.h5 volume mesh.
    """
    if _uses_container:
        mounted_volume = pyfluent.EXAMPLES_PATH
        work_dir_meshing = os.path.join(mounted_volume, "tmp_meshing-fluid")
    else:
        work_dir_meshing = os.path.join(workdir, "meshing-fluid")

    if not os.path.isdir(work_dir_meshing):
        os.makedirs(work_dir_meshing)
    else:
        files = glob.glob(os.path.join(work_dir_meshing, "*.stl"))
        for f in files:
            os.remove(f)

    # write all boundaries
    for b in fluid_boundaries:
        filename = os.path.join(work_dir_meshing, b.name.lower() + ".stl")
        b.save(filename)
        add_solid_name_to_stl(filename, b.name.lower(), file_type="binary")

    for c in caps:
        filename = os.path.join(work_dir_meshing, c.name.lower() + ".stl")
        c.save(filename)
        add_solid_name_to_stl(filename, c.name.lower(), file_type="binary")

    session = _get_fluent_meshing_session()

    # import all stls
    if _uses_container:
        # NOTE: when using a Fluent container visible files
        # will be in /mnt/pyfluent. So need to use relative paths
        # or replace dirname by /mnt/pyfluent as prefix
        work_dir_meshing = "."
    session.tui.file.import_.cad(f"no {work_dir_meshing} *.stl")

    # merge objects
    session.tui.objects.merge("'(*)", "model-fluid")

    # fix duplicate nodes
    session.tui.diagnostics.face_connectivity.fix_free_faces("objects '(*)")

    # set size field
    session.tui.size_functions.set_global_controls(1, 1, 1.2)
    session.tui.scoped_sizing.compute("yes")

    # remesh all caps
    if remesh_caps:
        session.tui.boundary.remesh.remesh_constant_size("(cap_*)", "()", 40, 20, 1, "yes")

    # convert to mesh object
    session.tui.objects.change_object_type("(*)", "mesh", "yes")

    # compute volumetric regions
    session.tui.objects.volumetric_regions.compute("model-fluid")

    # mesh volume
    session.tui.mesh.auto_mesh("model-fluid")

    # clean up
    session.tui.objects.delete_all_geom()
    session.tui.objects.delete_unreferenced_faces_and_edges()

    # write
    file_path_mesh = os.path.join(workdir, "fluid-mesh.msh.h5")
    session.tui.file.write_mesh(file_path_mesh)

    mesh = FluentMesh(file_path_mesh)
    mesh.load_mesh()

    return mesh


def mesh_from_manifold_input_model(
    model: _InputModel,
    workdir: Union[str, Path],
    path_to_output: Union[str, Path],
    mesh_size: float = 2.0,
    overwrite_existing_mesh: bool = True,
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
    smooth_boundaries = False
    fix_intersections = False
    auto_improve_nodes = False

    if not isinstance(model, _InputModel):
        raise ValueError(f"Expecting input to be of type {str(_InputModel)}")

    if not os.path.isdir(workdir):
        os.makedirs(workdir)

    # NOTE: when using containerized version - we need to copy all the files
    # to and from the mounted volume given by pyfluent.EXAMPLES_PATH (default)
    if _uses_container:
        mounted_volume = pyfluent.EXAMPLES_PATH
        work_dir_meshing = os.path.join(mounted_volume)
    else:
        work_dir_meshing = os.path.abspath(os.path.join(workdir, "meshing"))

    if os.path.isdir(work_dir_meshing) and not _uses_container:
        shutil.rmtree(work_dir_meshing)

    try:
        os.makedirs(work_dir_meshing)
    except:
        LOGGER.debug("Failed to create working directory")

    LOGGER.debug(f"Path to meshing directory: {work_dir_meshing}")

    if not os.path.isfile(path_to_output) or overwrite_existing_mesh:

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
        LOGGER.debug(f"Writing input files in: {work_dir_meshing}")
        model.write_part_boundaries(work_dir_meshing)

        session = _get_fluent_meshing_session()

        session.transcript.start(
            os.path.join(work_dir_meshing, "fluent_meshing.log"), write_to_stdout=False
        )

        # import files
        if _uses_container:
            # NOTE: when using a Fluent container visible files
            # will be in /mnt/pyfluent. So need to use relative paths
            # or replace dirname by /mnt/pyfluent as prefix
            work_dir_meshing = "."

        session.tui.file.import_.cad('no "' + work_dir_meshing + '" "*.stl" yes 40 yes mm')
        session.tui.objects.merge("'(*) heart")
        session.tui.objects.labels.create_label_per_zone("heart '(*)")
        session.tui.diagnostics.face_connectivity.fix_free_faces(
            "objects '(*) merge-nodes yes 1e-3"
        )

        if fix_intersections:
            session.tui.diagnostics.face_connectivity.fix_self_intersections(
                "objects '(heart) fix-self-intersection"
            )

        # smooth all zones
        face_zone_names = _get_face_zones_with_filter(session, "*")

        if smooth_boundaries:
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
        if fix_intersections:
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

        if auto_improve_nodes:
            session.tui.mesh.modify.auto_node_move("(*)", "(*)", 0.3, 50, 120, "yes", 5)

        session.tui.objects.delete_all_geom()
        session.tui.mesh.zone_names_clean_up()
        # session.tui.mesh.check_mesh()
        # session.tui.mesh.check_quality()
        session.tui.boundary.manage.remove_suffix("(*)")

        session.tui.mesh.prepare_for_solve("yes")

        # write to file

        if _uses_container:
            session.tui.file.write_mesh(os.path.basename(path_to_output))
        else:
            session.tui.file.write_mesh('"' + path_to_output + '"')
        # session.meshing.tui.file.read_journal(script)
        session.exit()

        if path_to_output != path_to_output_old:
            shutil.copy(path_to_output, path_to_output_old)

        path_to_output = path_to_output_old
    else:
        LOGGER.debug(f"Reusing: {path_to_output}")

    mesh = FluentMesh()
    mesh.load_mesh(path_to_output)
    mesh._fix_negative_cells()

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
    overwrite_existing_mesh: bool = True,
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
        work_dir_meshing = os.path.join(mounted_volume)
    else:
        work_dir_meshing = os.path.abspath(os.path.join(workdir, "meshing"))

    if os.path.isdir(work_dir_meshing) and not _uses_container:
        shutil.rmtree(work_dir_meshing)

    try:
        os.makedirs(work_dir_meshing)
    except:
        LOGGER.debug("Failed to create working directory")

    if not os.path.isfile(path_to_output) or overwrite_existing_mesh:
        path_to_output_old = path_to_output
        path_to_output = os.path.join(work_dir_meshing, "volume-mesh.msh.h5")

        min_size = mesh_size
        max_size = mesh_size
        growth_rate = 1.2

        # clean up any stls in the directory
        stls = glob.glob(os.path.join(work_dir_meshing, "*.stl"))
        for stl in stls:
            os.remove(stl)

        for part in model.parts:
            part.name = part.name.lower().replace(" ", "_")

        # write all boundaries
        LOGGER.debug(f"Writing input files in: {work_dir_meshing}")
        model.write_part_boundaries(work_dir_meshing, add_name_to_header=False)

        # launch pyfluent
        session = _get_fluent_meshing_session()

        session.transcript.start(
            os.path.join(work_dir_meshing, "fluent_meshing.log"), write_to_stdout=False
        )

        # # import stls
        if _uses_container:
            # NOTE: when using a Fluent container visible files
            # will be in /mnt/pyfluent. So need to use relative paths
            # or replace dirname by /mnt/pyfluent as prefix
            work_dir_meshing = "."

        session.tui.file.import_.cad("no", work_dir_meshing, "*.stl", "yes", 40, "yes", "mm")

        # each stl is imported as a separate object. Wrap the different collections of stls to
        # create new surface meshes for each of the parts.
        session.tui.size_functions.set_global_controls(min_size, max_size, growth_rate)
        session.tui.scoped_sizing.compute('"yes"')

        session.tui.objects.extract_edges("'(*) feature 40")

        for part in model.parts:
            LOGGER.info("Wrapping " + part.name + "...")
            # wrap object.
            _wrap_part(session, part.boundary_names, part.name)

        # wrap entire model in one pass so that we can create a single volume mesh. Use list of all
        # input boundaries are given as input. External material point for meshing.
        # NOTE: this assumes that all the individually wrapped parts form a single
        # connected structure.
        LOGGER.info("Wrapping model...")
        _wrap_part(session, model.boundary_names, "model")

        # mesh the entire model in one go.
        session.tui.objects.volumetric_regions.compute("model")
        session.tui.mesh.auto_mesh("model yes pyramids tet no")

        # clean up geometry objects
        session.tui.objects.delete_all_geom()

        # write mesh
        if os.path.isfile(path_to_output):
            os.remove(path_to_output)

        if _uses_container:
            session.tui.file.write_mesh(os.path.basename(path_to_output))
        else:
            session.tui.file.write_mesh('"' + path_to_output + '"')
        session.exit()

        shutil.copy(path_to_output, path_to_output_old)

        path_to_output = path_to_output_old
    else:
        LOGGER.debug(f"Reusing {path_to_output}")
        for part in model.parts:
            part.name = part.name.replace(" ", "_").lower()

    # Update the cell zones such that for each part we have a separate cell zone.
    mesh = FluentMesh()
    mesh.load_mesh(path_to_output)
    mesh._fix_negative_cells()

    # convert to unstructured grid.
    grid = mesh._to_vtk()

    # represent cell centroids as point cloud assign part-ids to cells.
    cell_centroids = grid.cell_centers()
    cell_centroids.point_data.set_scalars(name="part-id", scalars=0)

    # NOTE: should use wrapped surfaces to select part.
    # assign wrapped boundaries to input parts.
    for ii, part in enumerate(model.parts):
        face_zones_wrapped = [fz for fz in mesh.face_zones if part.name + ":" in fz.name]
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

    # use individual wrapped parts to separate the parts of the wrapped model.
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
        part.name = part.name.replace("_", " ").capitalize()
        cell_zone = FluentCellZone(
            min_id=np.argwhere(partids_sorted == part.id)[0][0],
            max_id=np.argwhere(partids_sorted == part.id)[-1][0],
            name=part.name,
            cid=part.id,
        )
        cell_zone.get_cells(new_mesh.cells)
        new_mesh.cell_zones.append(cell_zone)

    # keep just the face zones of the entire wrapped model and the corresponding
    # interior face zone
    new_mesh.face_zones = [
        fz
        for fz in new_mesh.face_zones
        if "model:" in fz.name.lower() or "interior-" in fz.name.lower()
    ]

    # rename face zones - rename to original input names.
    for fz in new_mesh.face_zones:
        if "interior" in fz.name:
            continue
        fz.name = fz.name.replace("model:", "")
        if ":" in fz.name:
            fz.name = fz.name.split(":")[0]

    return new_mesh


if __name__ == "__main__":
    LOGGER.info("Protected")
