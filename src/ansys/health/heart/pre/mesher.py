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

"""Module containing methods for interaction with Fluent meshing."""

import glob
import os
from pathlib import Path
import shutil

import numpy as np
import pyvista as pv

import ansys.fluent.core as pyfluent
from ansys.fluent.core.launcher.launch_options import LaunchMode
from ansys.fluent.core.session_meshing import Meshing as MeshingSession
import ansys.fluent.core.utils.file_transfer_service as file_transfer_service
from ansys.health.heart import LOG as LOGGER
from ansys.health.heart.exceptions import SupportedFluentVersionNotFoundError
from ansys.health.heart.objects import Mesh, SurfaceMesh
from ansys.health.heart.pre.input import _InputBoundary, _InputModel
from ansys.health.heart.utils.fluent_reader import _FluentCellZone, _FluentMesh
from ansys.health.heart.utils.vtk_utils import (
    add_solid_name_to_stl,
    cell_ids_inside_enclosed_surface,
)
import ansys.platform.instancemanagement as pypim

_supported_fluent_versions = ["25.2", "25.1", "24.2", "24.1"]
"""List of supported Fluent versions."""
_supported_fluent_versions_container = ["25.2", "24.2", "24.1"]

_num_cpus: bool = 2
"""Number of CPUs to use for meshing."""
_extra_launch_kwargs = {}
"""Extra keyword arguments passed to ``pyfluent.launch_fluent()``."""
_fluent_version = None
"""Global variable to explicitly override the Fluent version used."""

_launch_mode: LaunchMode = None
"""Fluent Launch mode."""

_uses_container: bool = False
"""Global variable to switch to Fluent container mode."""

_fluent_ui_mode = pyfluent.UIMode(os.getenv("PYFLUENT_UI_MODE", pyfluent.UIMode.HIDDEN_GUI))


LOGGER.debug(f"Fluent user interface mode: {_fluent_ui_mode.value}")


def _get_supported_fluent_version() -> str:
    """Use PyFluent to get a supported Fluent version."""
    if os.getenv("PYANSYS_HEART_FLUENT_VERSION", None):
        version = os.getenv("PYANSYS_HEART_FLUENT_VERSION")
        if version not in _supported_fluent_versions:
            raise SupportedFluentVersionNotFoundError(
                f"Fluent version {version} is not supported. Supported versions are: {_supported_fluent_versions}"  # noqa: E501
            )
        return version

    for version in _supported_fluent_versions:
        try:
            pyfluent.launch_fluent(product_version=version, dry_run=True)
            LOGGER.info(
                f"Found Fluent {version} as latest compatible "
                + f"version from supported versions: {_supported_fluent_versions}."
            )
            return version
        except Exception:
            pass
    raise SupportedFluentVersionNotFoundError(
        f"""Did not find a supported Fluent version.
        Install one of these versions: {_supported_fluent_versions}"""
    )


def _get_face_zones_with_filter(pyfluent_session, prefixes: list) -> list[str]:
    """Get a list of available boundaries in a Fluent session that uses any of the prefixes."""
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


def _organize_connected_regions(
    grid: pv.UnstructuredGrid, scalar: str = "part-id"
) -> pv.UnstructuredGrid:
    """Ensure that cells that belong to same part are connected."""
    LOGGER.debug("Re-organize connected regions.")
    part_ids = np.unique(grid.cell_data[scalar])
    grid.cell_data["orig-cell-ids"] = np.arange(0, grid.n_cells, dtype=int)

    # grid1 = grid.copy(deep=False)
    grid1 = grid.copy(deep=True)
    tets = grid1.cells_dict[10]

    # get a list of orphan cells
    orphan_cell_ids = []
    for part_id in part_ids:
        mask = grid1.cell_data[scalar] == part_id
        grid2 = grid1.extract_cells(mask)

        #! Note that this may not suffice. Need to check if at least => 2 points are connected
        conn = grid2.connectivity()
        num_regions = np.unique(conn.cell_data["RegionId"]).shape[0]

        if num_regions == 1:
            continue

        LOGGER.debug(f"Found {num_regions - 1} unnconnected regions. Find connected candidate.")
        # for each region, find to what "main" region it is connected.
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
                LOGGER.debug("More than one candidate.")

            grid.cell_data["part-id"][orphan_cell_ids] = unique_ids[np.argmax(counts)]

    return grid


def _assign_part_id_to_orphan_cells(
    grid: pv.UnstructuredGrid, scalar="part-id"
) -> pv.UnstructuredGrid:
    """Use closest point interpolation to assign part ID to orphan cells."""
    grid.cell_data["_original-cell-ids"] = np.arange(0, grid.n_cells)
    orphans = grid.extract_cells(grid.cell_data[scalar] == 0)

    if orphans.n_cells == 0:
        LOGGER.debug("No orphan cells are detected.")
        return grid

    LOGGER.debug(f"Assigning part IDs to {orphans.n_cells} orphan cells...")

    grid_centers = grid.cell_centers()
    grid_centers = grid_centers.extract_points(grid_centers.cell_data["part-id"] != 0)
    part_ids = grid_centers.point_data["part-id"]

    grid_centers.point_data.clear()
    grid_centers.cell_data.clear()
    grid_centers.point_data["part-id"] = part_ids

    orphan_centers = orphans.cell_centers()
    orphan_centers.cell_data.clear()
    orphan_centers.point_data.clear()

    orphan_centers = orphan_centers.interpolate(grid_centers, n_points=1)

    interpolated_part_ids = np.array(orphan_centers.point_data["part-id"], dtype=np.int32)
    # assign interpolated part ids again to the original grid.
    grid2 = grid.copy()
    grid2.cell_data["part-id"][orphans.cell_data["_original-cell-ids"]] = interpolated_part_ids
    return grid2


def _get_cells_inside_wrapped_parts(model: _InputModel, mesh: _FluentMesh) -> pv.UnstructuredGrid:
    """Get cells inside each of the wrapped parts."""
    grid = mesh._to_vtk()

    # represent cell centroids as point cloud
    cell_centroids = grid.cell_centers()
    cell_centroids.point_data.set_scalars(name="part-id", scalars=0)
    cell_centroids.point_data["_original-cell-ids"] = np.arange(0, cell_centroids.n_points)
    cell_centroids1 = cell_centroids.copy()

    used_cell_ids = np.empty(shape=(0,), dtype=int)
    # use individual wrapped parts to separate the parts of the wrapped model.
    for part in model.parts:
        if not part.is_manifold:
            LOGGER.warning(f"Part {part.name} is not manifold.")
        LOGGER.debug(f"Redistributing cells based on surface of {part.name}...")

        # get surface
        surface = part.combined_boundaries.clean()

        point_ids_inside = cell_ids_inside_enclosed_surface(
            cell_centroids1, surface, tolerance=1e-9
        )

        # map back to original vtk object.
        point_ids_inside1 = cell_centroids1.point_data["_original-cell-ids"][point_ids_inside]
        cell_centroids.point_data["part-id"][point_ids_inside1] = part.id

        # reduce data set for efficiency. Extract points/cells is very slow for some reason.
        used_cell_ids = np.append(used_cell_ids, point_ids_inside1)
        point_ids_not_inside = np.setdiff1d(
            cell_centroids.point_data["_original-cell-ids"], used_cell_ids
        )
        # redefine cell centroids polydata. for processing.
        cell_centroids1 = pv.PolyData(cell_centroids.points[point_ids_not_inside, :])
        cell_centroids1.point_data["_original-cell-ids"] = cell_centroids.point_data[
            "_original-cell-ids"
        ][point_ids_not_inside]

    grid1 = grid.copy()
    grid1.cell_data["part-id"] = np.array(cell_centroids.point_data["part-id"], dtype=int)

    # Use closed point interpolation to assign part ids to any un-assigned (orphan) cells.
    grid2 = _assign_part_id_to_orphan_cells(grid1)

    if np.any(grid2.cell_data["part-id"] == 0):
        LOGGER.warning("Not all cells have a part ID assigned.")

    return grid2


def _get_fluent_meshing_session(working_directory: str | Path) -> MeshingSession:
    """Get a Fluent Meshing session."""
    # NOTE: There are three launch modes Fluent can be launched in:
    # 1. LaunchMode.PIM: Fluent is launched using the Product Instance Management (PIM) service.
    # 2. LaunchMode.CONTAINER: Fluent is launched in a container. (containerized mode)
    # 3. LaunchMode.STANDALONE: Fluent is launched as a standalone application. (fallback mode)
    # File transfer strategies are different for each mode.

    # check whether containerized version of Fluent is used
    global _uses_container
    global _launch_mode
    global _supported_fluent_versions

    _uses_container = bool(int(os.getenv("PYFLUENT_LAUNCH_CONTAINER", False)))
    if _uses_container:
        _supported_fluent_versions = _supported_fluent_versions_container

    num_cpus = int(os.getenv("PYANSYS_HEART_NUM_CPU", _num_cpus))

    if _fluent_version is None:
        product_version = _get_supported_fluent_version()
    else:
        product_version = _fluent_version

    # determine launch mode
    if pypim.is_configured():
        _launch_mode = LaunchMode.PIM
        transfer_strategy = None

    elif _uses_container:
        _launch_mode = LaunchMode.CONTAINER
        transfer_strategy = None
    else:
        _launch_mode = LaunchMode.STANDALONE
        transfer_strategy = file_transfer_service.StandaloneFileTransferStrategy()

    LOGGER.info(f"Launching meshing session with {product_version}...")

    launch_config = {
        "precision": pyfluent.Precision.DOUBLE,
        "processor_count": num_cpus,
        "start_transcript": False,
        "product_version": product_version,
        "ui_mode": _fluent_ui_mode,
        "file_transfer_service": transfer_strategy,
    }

    match _launch_mode:
        case LaunchMode.PIM:
            launch_config["ui_mode"] = None
            LOGGER.info(f"Launching Fluent in PIM-mode with config: {launch_config}")
            session = pyfluent.PureMeshing.from_pim(**launch_config, **_extra_launch_kwargs)

        case LaunchMode.CONTAINER:
            LOGGER.info(f"Launching Fluent in Container mode with config: {launch_config}")
            launch_config["container_dict"] = {
                "mount_source": f"{working_directory}",
                "mount_target": "/mnt/pyfluent/meshing",
            }
            launch_config["ui_mode"] = pyfluent.UIMode.NO_GUI_OR_GRAPHICS
            session = pyfluent.PureMeshing.from_container(**launch_config, **_extra_launch_kwargs)

        case LaunchMode.STANDALONE:
            LOGGER.info(f"Launching Fluent in Standalone mode with config: {launch_config}")
            session = pyfluent.PureMeshing.from_install(**launch_config, **_extra_launch_kwargs)

    return session


def _wrap_part(session: MeshingSession, boundary_names: list, wrapped_part_name: str) -> list[str]:
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


def _to_fluent_convention(string_to_convert: str) -> str:
    """Convert string to the Fluent-supported convention."""
    return string_to_convert.lower().replace(" ", "_")


def _update_size_per_part(
    part_names: list[str],
    global_size: float,
    size_per_part: dict = None,
) -> dict:
    """Update the dictionary containing the (wrap) size per part.

    Parameters
    ----------
    global_size : float
        Global size to use for parts that are not referenced.
    part_names : list[str]
        Part names involved in the model.
    size_per_part : dict, default: None
        Size per part used to override global size.
    """
    # convert both to Fluent-naming convention. Note: remove cases and spaces
    part_names = [_to_fluent_convention(part) for part in part_names]
    if size_per_part is not None:
        size_per_part = {_to_fluent_convention(part): size for part, size in size_per_part.items()}

    mesh_size_per_part = {part_name: global_size for part_name in part_names}

    if size_per_part is not None:
        for part, size in size_per_part.items():
            if part in part_names:
                mesh_size_per_part[part] = size

    return mesh_size_per_part


def _update_input_model_with_wrapped_surfaces(
    model: _InputModel, mesh: _FluentMesh, face_zone_ids_per_part: dict
) -> _InputModel:
    """Update the input model with the wrapped surfaces.

    Parameters
    ----------
    model : _InputModel
        Input model to bupdate.
    mesh : FluentMesh
        Fluent mesh containing all wrapped face zones.
    face_zone_ids_per_part : dict
        Face zone IDs for each part.

    Returns
    -------
    _InputModel
        Input model with wrapped boundaries assigned.
    """
    for ii, part in enumerate(model.parts):
        face_zone_ids_per_part[part.name]
        face_zones_wrapped = [
            fz for fz in mesh.face_zones if fz.id in face_zone_ids_per_part[part.name]
        ]
        if len(face_zones_wrapped) == 0:
            LOGGER.error(f"Did not find any wrapped face zones for {part.name}.")

        # replace with remeshed face zones, note that we may have more face zones now.
        remeshed_boundaries = []
        for fz in face_zones_wrapped:
            face_zone_name = fz.name.replace(part.name + ":", "")

            # try to maintain the input id.
            try:
                boundary_id = model.boundary_ids[model.boundary_names.index(face_zone_name)]
            except IndexError:
                boundary_id = fz.id

            remeshed_boundary = _InputBoundary(
                mesh.nodes,
                faces=np.hstack([np.ones(fz.faces.shape[0], dtype=int)[:, None] * 3, fz.faces]),
                id=boundary_id,
                name=face_zone_name,
            )
            remeshed_boundaries.append(remeshed_boundary)

        model.parts[ii].boundaries = remeshed_boundaries

    return model


def _post_meshing_cleanup(fluent_mesh: _FluentMesh) -> Mesh:
    """Clean up and retrieve VTK mesh after meshing."""
    # remove any unused nodes
    fluent_mesh.clean()

    # convert to vtk grid.
    vtk_grid = fluent_mesh._to_vtk()
    mesh = Mesh(vtk_grid)

    # get mapping from fluent mesh.
    mesh.cell_data["_volume-id"] = mesh.cell_data["cell-zone-ids"]
    mesh._volume_id_to_name = fluent_mesh.cell_zone_id_to_name

    # merge face zones based on fluent naming convention.
    fluent_mesh._merge_face_zones_based_on_connectivity(face_zone_separator=":")

    # add face zones to the mesh.
    for fz in fluent_mesh.face_zones:
        if "interior" not in fz.name:
            surface = SurfaceMesh(
                name=fz.name, triangles=fz.faces, nodes=fluent_mesh.nodes, id=fz.id
            )
            mesh.add_surface(surface, int(fz.id), name=fz.name)

    mesh = mesh.clean()

    return mesh


def _set_size_field_on_mesh_part(
    session: MeshingSession, mesh_size: float, part_name: str, growth_rate: float = 1.2
) -> MeshingSession:
    """Set the size field per part."""
    session.tui.scoped_sizing.create(
        f"boi-{part_name}",
        "boi",
        "object-faces-and-edges",
        "no",
        "yes",
        f"{part_name.lower()}",
        str(mesh_size),
        str(growth_rate),
    )

    return session


def _set_size_field_on_face_zones(
    session: MeshingSession,
    mesh_size: float,
    face_zone_names: list[str],
    boi_name: str,
    growth_rate: float = 1.2,
) -> MeshingSession:
    """Set the size field per list of boundaries."""
    session.tui.scoped_sizing.create(
        f"{boi_name}",
        "boi",
        "face-zone",
        "yes",
        "no",
        '"' + " ".join(face_zone_names).replace("'", '"') + '"',
        mesh_size,
        growth_rate,
    )
    return session


# TODO: fix method.
def _mesh_fluid_cavities(
    fluid_boundaries: list[SurfaceMesh],
    caps: list[SurfaceMesh],
    workdir: str,
    remesh_caps: bool = True,
) -> _FluentMesh:
    """Mesh the fluid cavities.

    Parameters
    ----------
    fluid_boundaries : List[SurfaceMesh]
        List of fluid boundaries used for meshing.
    caps : List[SurfaceMesh]
        List of caps that close each of the cavities.
    workdir : str
        Working directory.
    remesh_caps : bool, default: True
        whether to remesh the caps.

    Returns
    -------
    Path
        Path to the ``.msh.h5`` volume mesh.
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

    session = _get_fluent_meshing_session(work_dir_meshing)

    if _launch_mode == LaunchMode.PIM:
        # Upload files to session if in PIM or Container modes.
        LOGGER.info(f"Uploading files to session with working directory {work_dir_meshing}...")
        files = glob.glob(os.path.join(work_dir_meshing, "*.stl"))
        for file in files:
            session.upload(file)
        # In PIM mode files are uploaded to the Fluents working directory.
        work_dir_meshing = "."

    elif _launch_mode == LaunchMode.CONTAINER:
        # NOTE: when using a Fluent container visible files
        # will be in /mnt/pyfluent. (equal to mount target)
        work_dir_meshing = "/mnt/pyfluent/meshing"

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

    mesh = _FluentMesh(file_path_mesh)
    mesh.load_mesh()

    return mesh


def mesh_from_manifold_input_model(
    model: _InputModel,
    workdir: str | Path,
    path_to_output: str | Path,
    mesh_size: float = 2.0,
    overwrite_existing_mesh: bool = True,
) -> Mesh:
    """Create mesh from a good-quality manifold input model.

    Parameters
    ----------
    model : _InputModel
        Input model.
    workdir : Union[str, Path]
        Working directory.
    path_to_output : Union[str, Path]
        Path to the resulting Fluent mesh file.
    mesh_size : float, default: 2.0
        Uniform mesh size to use for both wrapping and filling the volume.

    Returns
    -------
    Mesh
        VTK mesh with both cell and face zones.
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
    except Exception as e:
        LOGGER.error(f"Failed to create working directory. {e}")

    LOGGER.info(f"Path to meshing directory: {work_dir_meshing}")

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
        LOGGER.info(f"Writing input files in: {work_dir_meshing}")
        model.write_part_boundaries(work_dir_meshing)

        session = _get_fluent_meshing_session(work_dir_meshing)

        session.transcript.start(
            os.path.join(work_dir_meshing, "fluent_meshing.log"), write_to_stdout=False
        )

        if _launch_mode == LaunchMode.PIM:
            # Upload files to session if in PIM or Container modes.
            LOGGER.info(f"Uploading files to session with working directory {work_dir_meshing}...")
            files = glob.glob(os.path.join(work_dir_meshing, "*.stl"))
            for file in files:
                session.upload(file)
            # In PIM mode files are uploaded to the Fluents working directory.
            work_dir_meshing = "."

        elif _launch_mode == LaunchMode.CONTAINER:
            # NOTE: when using a Fluent container visible files
            # will be in /mnt/pyfluent. (equal to mount target)
            work_dir_meshing = "/mnt/pyfluent/meshing"

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

        LOGGER.info(f"Writing mesh to {path_to_output}...")

        if _launch_mode in [LaunchMode.CONTAINER, LaunchMode.PIM]:
            session.tui.file.write_mesh(os.path.basename(path_to_output))
        else:
            session.tui.file.write_mesh('"' + path_to_output + '"')

        LOGGER.info(f"Copying {path_to_output} to {path_to_output_old}...")

        if _launch_mode == LaunchMode.PIM:
            session.download(os.path.basename(path_to_output), path_to_output_old)
        else:
            shutil.copy(path_to_output, path_to_output_old)

        if not os.path.isfile(path_to_output_old):
            raise FileNotFoundError(
                f"Failed to copy {os.path.basename(path_to_output)} to {path_to_output_old}. "
                "Please check the Fluent meshing log for errors."
            )

        session.exit()

        path_to_output = path_to_output_old
    else:
        LOGGER.info(f"Reusing: {path_to_output}")

    mesh = _FluentMesh()
    mesh.load_mesh(path_to_output)
    mesh._fix_negative_cells()

    # use part definitions to find which cell zone belongs to which part.
    for input_part in model.parts:
        surface = input_part.combined_boundaries

        if not surface.is_manifold:
            LOGGER.warning(
                "Part {0} is not manifold. Disabled surface check.".format(input_part.name)
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

    # Use only cell zones that are inside the parts defined in the input.
    mesh.cell_zones = [cz for cz in mesh.cell_zones if cz.id in model.part_ids]

    vtk_mesh = _post_meshing_cleanup(mesh)

    return vtk_mesh


def mesh_from_non_manifold_input_model(
    model: _InputModel,
    workdir: str | Path,
    path_to_output: str | Path,
    global_mesh_size: float = 2.0,
    _global_wrap_size: float = 1.5,
    overwrite_existing_mesh: bool = True,
    mesh_size_per_part: dict = None,
    _wrap_size_per_part: dict = None,
) -> Mesh:
    """Generate mesh from a non-manifold poor quality input model.

    Parameters
    ----------
    model : _InputModel
        Input model.
    workdir : Union[str, Path]
        Working directory.
    path_to_output : Union[str, Path]
        Path to the resulting Fluent mesh file.
    global_mesh_size : float, default: 2.0
        Uniform mesh size to use for all volumes and surfaces.
    _global_wrap_size : float, default: 1.5
        Global size used by the wrapper to reconstruct the geometry.
    overwrite_existing_mesh : bool, default: True
        Whether to overwrite an existing mesh.
    mesh_size_per_part : dict, default: None
        Dictionary specifying the mesh size that should be used for each part.
    _wrap_size_per_part : dict, default: None
        Dictionary specifying the wrap size that should be used to wrap each part.

    Notes
    -----
    This method uses Fluent wrapping technology to wrap the individual parts. First it
    creates manifold parts. Then, it consequently wraps the entire model and uses the manifold
    parts to split the wrapped model into the different cell zones.

    When specifying a mesh size per part, you can do that by either specifying the size for all
    parts or for specific parts. The default mesh size is used for any part not listed
    in the dictionary. This also applies to the wrapping step. You can control the wrap size
    per part or on a global level. By default, a size of 1.5 mm is used, but this value is not
    guaranteed to give good results.

    Note that a post-wrap remesh is triggered if the wrap size is not equal to the target mesh size.
    Remeshing might fail if the target mesh size deviates too much from the wrap size.

    Returns
    -------
    Mesh
        VTK mesh with both cell and face zones.
    """
    if not isinstance(model, _InputModel):
        raise ValueError(f"Expecting input to be of type {str(_InputModel)}.")

    mesh_size_per_part = _update_size_per_part(
        model.part_names, global_mesh_size, mesh_size_per_part
    )
    _wrap_size_per_part = _update_size_per_part(
        model.part_names, _global_wrap_size, _wrap_size_per_part
    )

    # Flag to determine whether to do a post-wrap remesh.
    if _wrap_size_per_part == mesh_size_per_part:
        post_wrap_remesh = False
    else:
        post_wrap_remesh = True

    min_mesh_size = np.min(list(mesh_size_per_part.values()))
    max_mesh_size = np.max(list(mesh_size_per_part.values()))

    min_mesh_size_wrap = np.min(list(_wrap_size_per_part.values()))
    max_mesh_size_wrap = np.max(list(_wrap_size_per_part.values()))

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
    except Exception as e:
        LOGGER.error(f"Failed to create working directory. {e}")

    if not os.path.isfile(path_to_output) or overwrite_existing_mesh:
        path_to_output_old = path_to_output
        path_to_output = os.path.join(work_dir_meshing, "volume-mesh.msh.h5")

        growth_rate = 1.2

        # clean up any stls in the directory
        stls = glob.glob(os.path.join(work_dir_meshing, "*.stl"))
        for stl in stls:
            os.remove(stl)

        # convert model names to Fluent-supported convention.
        for part in model.parts:
            part.name = _to_fluent_convention(part.name)

        # write all boundaries
        LOGGER.info(f"Writing input files in: {work_dir_meshing}")
        model.write_part_boundaries(work_dir_meshing, add_name_to_header=False)

        # launch pyfluent
        session = _get_fluent_meshing_session(work_dir_meshing)

        LOGGER.info(f"Starting Fluent Meshing in mode: {_launch_mode}")

        session.transcript.start(
            os.path.join(work_dir_meshing, "fluent_meshing.log"), write_to_stdout=False
        )

        if _launch_mode == LaunchMode.PIM:
            # Upload files to session if in PIM or Container modes.
            LOGGER.info(f"Uploading files to session with working directory {work_dir_meshing}...")
            files = glob.glob(os.path.join(work_dir_meshing, "*.stl"))
            for file in files:
                session.upload(file)
            # In PIM mode files are uploaded to the Fluents working directory.
            work_dir_meshing = "."

        elif _launch_mode == LaunchMode.CONTAINER:
            # NOTE: when using a Fluent container visible files
            # will be in /mnt/pyfluent. (equal to mount target)
            work_dir_meshing = "/mnt/pyfluent/meshing"

        session.tui.file.import_.cad("no", work_dir_meshing, "*.stl", "yes", 40, "yes", "mm")

        # each stl is imported as a separate object. Wrap the different collections of stls to
        # create new surface meshes for each of the parts.

        ## Set size field for wrapping process.
        #####################################################################
        session.tui.size_functions.set_global_controls(
            min_mesh_size_wrap, max_mesh_size_wrap, growth_rate
        )

        ## Set up boi's scoped to face zones of individual parts:
        part_names_input = [p.name for p in model.parts]
        for part, wrap_size in _wrap_size_per_part.items():
            idx = part_names_input.index(part)
            try:
                b_names = [b.name.replace("'", '"') for b in model.parts[idx].boundaries]
                _set_size_field_on_face_zones(
                    session, wrap_size, b_names, f"boi-{part}", growth_rate
                )
            except Exception as e:
                LOGGER.warning(f"Failed to set mesh size for {part}: {e}")

        session.tui.scoped_sizing.compute('"yes"')
        #####################################################################

        session.tui.objects.extract_edges("'(*) feature 40")

        part_face_zone_ids_post_wrap = {}
        for part in model.parts:
            LOGGER.info("Wrapping " + part.name + "...")
            # wrap object.
            _wrap_part(session, part.boundary_names, part.name)
            session.tui.objects.volumetric_regions.compute(part.name)

            part_face_zone_ids_post_wrap[part.name] = session.scheme_eval.scheme_eval(
                f"(get-face-zones-of-objects '({part.name}) )"
            )

        # NOTE: wrap entire model in one pass so that a single volume mesh can be created.
        # Use list of all input boundaries as input. Uses external material point for meshing.
        # This assumes that all the individually wrapped parts form a single
        # connected structure.
        LOGGER.info("Wrapping model...")
        _wrap_part(session, model.boundary_names, "model")

        ## Recompute size field for final target mesh size. This size-field is used for
        ## remeshing the wrapped model.
        if post_wrap_remesh:
            LOGGER.info("Post-wrap remeshing of model surfaces...")
            session.tui.objects.delete_all_geom()
            session.tui.scoped_sizing.delete_size_field()
            session.tui.scoped_sizing.delete_all()

            ## set size field for final mesh.
            #####################################################################
            session.tui.size_functions.set_global_controls(
                str(min_mesh_size), str(max_mesh_size), str(growth_rate)
            )

            for part, mesh_size in mesh_size_per_part.items():
                idx = part_names_input.index(part)
                try:
                    b_names = [b.name.replace("'", '"') for b in model.parts[idx].boundaries]
                    _set_size_field_on_mesh_part(session, mesh_size, part, growth_rate)
                except Exception as e:
                    LOGGER.warning(f"Failed to set mesh size for {part}: {e}")

            session.tui.scoped_sizing.compute('"yes"')
            #####################################################################

            session.tui.boundary.remesh.remesh_face_zones_conformally(
                "'(model*)", "()", 40, 20, "yes"
            )

        # mesh the entire model in one go.
        session.tui.objects.volumetric_regions.compute("model")
        session.tui.mesh.tet.controls.cell_sizing("size-field")

        LOGGER.info("Generating volume mesh...")
        session.tui.mesh.auto_mesh("model yes pyramids tet no")

        # clean up geometry objects
        session.tui.objects.delete_all_geom()
        session.tui.objects.delete_unreferenced_faces_and_edges()

        # write mesh
        if os.path.isfile(path_to_output):
            os.remove(path_to_output)

        LOGGER.info(f"Writing mesh to {path_to_output}...")

        if _launch_mode in [LaunchMode.CONTAINER, LaunchMode.PIM]:
            session.tui.file.write_mesh(os.path.basename(path_to_output))
        else:
            session.tui.file.write_mesh('"' + path_to_output + '"')

        LOGGER.info(f"Copying {path_to_output} to {path_to_output_old}...")

        if _launch_mode == LaunchMode.PIM:
            session.download(os.path.basename(path_to_output), path_to_output_old)
        else:
            shutil.copy(path_to_output, path_to_output_old)

        if not os.path.isfile(path_to_output_old):
            raise FileNotFoundError(
                f"Failed to copy {os.path.basename(path_to_output)} to {path_to_output_old}. "
                "Please check the Fluent meshing log for errors."
            )

        session.exit()

        path_to_output = path_to_output_old
    else:
        LOGGER.debug(f"Reusing {path_to_output}")
        for part in model.parts:
            part.name = _to_fluent_convention(part.name)

    LOGGER.info("Post Fluent-Meshing cleanup...")
    # update cell zones such that each part has a separate cell zone
    mesh = _FluentMesh()
    mesh.load_mesh(path_to_output)
    mesh._fix_negative_cells()

    # update input model with wrapped surfaces
    model = _update_input_model_with_wrapped_surfaces(model, mesh, part_face_zone_ids_post_wrap)

    # get cells inside each of the wrapped parts.
    grid = _get_cells_inside_wrapped_parts(model, mesh)

    # ensure parts are continuous and well connected.
    grid = _organize_connected_regions(grid, scalar="part-id")

    if np.any(grid.cell_data["part-id"] == 0):
        raise ValueError("Invalid mesh, not all elements assigned to a part.")

    # TODO: refactor and cleanup the following
    # TODO: stick with the VTK object instead of _FluentMesh.
    # change FluentMesh object accordingly.
    idx_sorted = np.argsort(np.array(grid.cell_data["part-id"], dtype=int))
    partids_sorted = np.sort(np.array(grid.cell_data["part-id"], dtype=int))

    new_mesh = mesh
    new_mesh.cells = new_mesh.cells[idx_sorted]
    new_mesh.cell_zones: list[_FluentCellZone] = []

    for part in model.parts:
        # convert back to original convention.
        # TODO: refactor so that we revert back to the original name.
        part.name = part.name.replace("_", " ").capitalize()
        cell_zone = _FluentCellZone(
            min_id=np.argwhere(partids_sorted == part.id)[0][0],
            max_id=np.argwhere(partids_sorted == part.id)[-1][0],
            name=part.name,
            cid=part.id,
        )
        cell_zone.get_cells(new_mesh.cells)
        new_mesh.cell_zones.append(cell_zone)

    # keep only the face zones of the entire wrapped model and the corresponding
    # interior face zone
    new_mesh.face_zones = [
        fz
        for fz in new_mesh.face_zones
        if "model:" in fz.name.lower() or "interior-" in fz.name.lower()
    ]

    # rename face zones to original input names
    for fz in new_mesh.face_zones:
        if "interior" in fz.name:
            continue
        fz.name = fz.name.replace("model:", "")
        if ":" in fz.name:
            fz.name = fz.name.split(":")[0]

    # Use only cell zones that are inside the parts defined in the input
    new_mesh.cell_zones = [cz for cz in new_mesh.cell_zones if cz.id in model.part_ids]

    vtk_mesh = _post_meshing_cleanup(new_mesh)

    return vtk_mesh
