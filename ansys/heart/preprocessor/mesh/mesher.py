"""Module contains methods for interaction with Fluent meshing."""
import os
import subprocess
import glob
import shutil

from ansys.heart.custom_logging import LOGGER
from ansys.heart.preprocessor._load_template import load_template
import ansys.heart.preprocessor.mesh.fluenthdf5 as hdf5  # noqa: F401
import numpy as np

# from pkg_resources import resource_filename
import pkg_resources

_template_directory = pkg_resources.resource_filename("ansys.heart.preprocessor", "templates")


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
        meshing_mode=True,
        precision="double",
        processor_count=num_cpus,
        start_transcript=True,
        show_gui=show_gui,
    )
    assert session.check_health() == "SERVING"
    session.start_transcript()
    session._print_transcript
    # LOGGER.debug("Reading fluent journal file: {0}".format(script))
    min_size = mesh_size
    max_size = mesh_size
    growth_rate_wrap = 1.2
    add_blood_pool = False

    # import files
    session.meshing.tui.file.import_.cad("no " + work_dir_meshing + " part_*.stl yes 40 yes mm")
    session.meshing.tui.file.start_transcript(work_dir_meshing, "fluent_meshing.log")
    session.meshing.tui.objects.merge("'(*) heart")
    session.meshing.tui.objects.labels.create_label_per_zone("heart '(*)")
    session.meshing.tui.diagnostics.face_connectivity.fix_free_faces(
        "objects '(*) merge-nodes yes 1e-3"
    )
    session.meshing.tui.objects.create_intersection_loops("collectively '(*)")
    session.meshing.tui.boundary.feature.create_edge_zones("(*) fixed-angle 70 yes")

    # set up size field for wrapping
    session.meshing.tui.size_functions.set_global_controls(min_size, max_size, growth_rate_wrap)
    session.meshing.tui.scoped_sizing.compute("yes")

    # wrap objects
    session.meshing.tui.objects.wrap.wrap(
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
        # ; template script for auto-generating caps for all endocardial parts
        # ; and extracting blood pool volume
        # ; This first copies all the endocardial zones and septum and
        # ; and closes the the cavities based on the free faces
        # ; Note that the auto-patch utility may not work in all cases.
        # ; ------------------------------------------------------------------
        #
        session.meshing.tui.boundary.modify.auto_patch_holes()
        raise NotImplementedError("Adding blood pool is not implemented through PyFluent yet")
        tui = session.meshing.tui
        tui.objects.delete_all_geom()
        tui

    # compute volumetric regions
    session.meshing.tui.objects.volumetric_regions.compute("wrapped-myocardium", "no")
    session.meshing.tui.objects.volumetric_regions.change_type(
        "Wrapped-myocardium", "'(*)", "fluid"
    )
    session.meshing.tui.objects.volumetric_regions.change_type(
        "wrapped-myocardium", "(heart)", "solid"
    )

    # start auto meshing
    session.meshing.tui.mesh.tet.controls.cell_sizing("size-field")
    session.meshing.tui.mesh.auto_mesh("wrapped-myocardium", "yes", "pyramids", "tet", "no")
    session.meshing.tui.mesh.modify.auto_node_move("(*)", "(*)", 0.3, 50, 120, "yes", 5)
    session.meshing.tui.objects.delete_all_geom()
    session.meshing.tui.mesh.zone_names_clean_up()
    session.meshing.tui.mesh.check_mesh()
    session.meshing.tui.mesh.check_quality()
    session.meshing.tui.boundary.manage.remove_suffix("(*)")

    # write to file
    session.meshing.tui.file.write_mesh(path_to_output)
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
