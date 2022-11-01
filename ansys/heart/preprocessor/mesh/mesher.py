"""Module contains methods for interaction with Fluent meshing."""
import os
import pathlib
import subprocess
from typing import List

from ansys.heart.custom_logging import LOGGER
from ansys.heart.preprocessor import SC_EXE
from ansys.heart.preprocessor._load_template import load_template
import ansys.heart.preprocessor.mesh.fluenthdf5 as hdf5  # noqa: F401
from ansys.heart.preprocessor.mesh.objects import Cap
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
    Use Fluent meshing to wrap the surface and create tetrahedral mesh.

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

    # start Fluent session using PyFluent:
    # TODO: Catch errors in session
    session = pyfluent.launch_fluent(
        meshing_mode=True,
        precision="double",
        processor_count=num_cpus,
        start_transcript=False,
        show_gui=show_gui,
    )
    session.meshing.tui.file.read_journal(script)
    session.exit()

    # change back to old directory
    os.chdir(old_directory)

    return


def _shrink_by_spaceclaim(input, output):
    """Use SpaceClaim shrinkwrapping to wrap surface and create high quality surface mesh."""
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
    import ansys.fluent.core as pyfluent

    # from ansys.fluent.core.services import SurfaceDataType

    session = pyfluent.launch_fluent(start_instance=True, show_gui=True, meshing_mode=True)

    session.meshing.tui.file.read_mesh(
        "D:\\development\pyheart-lib\\pyheart-lib\downloads\\Strocchi2020_Demo1.2"
        "\\p05\\fluent_volume_mesh.msh.h5"
    )
    session.meshing.tui.switch_to_solution_mode("yes")
    field_data = session.field_data
    field_data.add_get_surfaces_request(
        surface_ids=[10296], provide_vertices=True, provide_faces=True
    )
    payload_data = field_data.get_fields()

    # zone_id = info["right-ventricle-epicardium"]["zone_id"]
    # data = session.field_data.get_surface_data("right-ventricle-epicardium",
    # session.field_data.add_get_surfaces_request([zone_id])
    data = session.field_data.get_fields()
    session.check_health()
    session.exit()

    LOGGER.info("Protected")
