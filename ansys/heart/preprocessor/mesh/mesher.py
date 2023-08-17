"""Module contains methods for interaction with Fluent meshing."""
import glob
import os
import shutil
import subprocess

from ansys.heart.custom_logging import LOGGER
from ansys.heart.preprocessor._load_template import load_template
import ansys.heart.preprocessor.mesh.fluenthdf5 as hdf5  # noqa: F401
import numpy as np

# from pkg_resources import resource_filename
import pkg_resources

_template_directory = pkg_resources.resource_filename("ansys.heart.preprocessor", "templates")

_fluent_version = "22.2.0"


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
