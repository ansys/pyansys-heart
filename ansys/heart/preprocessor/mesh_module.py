import os
import subprocess
import numpy as np
import gmsh
import pathlib


from ansys.heart.preprocessor._load_template import load_template
from ansys.heart.preprocessor import SC_EXE, FLUENT_EXE
from ansys.heart.custom_logging import logger

# for fluent:
import ansys.fluent.core as pyfluent

_template_directory = os.path.join(os.path.dirname(__file__), "template")

"""Module contains methods for mesh operations"""

## for remeshing purposes


def mesh_by_fluentmeshing(
    path_to_input_stl: str,
    path_to_output: str,
    mesh_size: float = 2.0,
    journal_type: str = "original",
):
    """Uses Fluent meshing to wrap the surface and create
    tetrahedral mesh"""

    if journal_type not in ["original", "simplified_geometry"]:
        raise ValueError("Journal type %s not found " % journal_type)

    # change directory to directory of stl file
    old_directory = os.getcwd()
    working_directory = pathlib.Path(path_to_input_stl).parent
    os.chdir(working_directory)

    if journal_type == "original":
        var_for_template = {
            "input_path": path_to_input_stl,
            "output_path": path_to_output,
            "mesh_size": mesh_size,
        }
        template = load_template("fluent_meshing_template.jou")

    elif journal_type == "simplified_geometry":
        var_for_template = {
            "work_directory": working_directory,
            "output_path": path_to_output,
            "mesh_size": mesh_size,
        }
        template = load_template("fluent_meshing_template_simplified.jou")

    script = os.path.join(os.path.dirname(path_to_input_stl), "fluent_meshing.jou")

    with open(script, "w") as f:
        f.write(template.render(var_for_template))

    num_cpus = 2

    # start Fluent session using PyFluent:
    # TODO: Catch errors in session
    session = pyfluent.launch_fluent(
        meshing_mode=True, precision="double", processor_count=num_cpus, start_transcript=False
    )
    session.meshing.tui.file.read_journal(script)
    session.exit()

    # change back to old directory
    os.chdir(old_directory)

    return


def shrink_by_spaceclaim(input, output):
    """Uses SpaceClaim shrinkwrapping to wrap surface and
    create high quality surface mesh"""

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


def run_gmsh(infile: str, outfile: str, mesh_size):
    """Runs GMESH with specified in/output file
    and target mesh size

    Args:
        infile (str): Path to .stl input file
        outfile (str): path to .vtk output
    """
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


def add_solid_name_to_stl(filename, solid_name, file_type: str = "ascii"):
    """Adds name of solid to stl file. Supports only single block"""
    if file_type == "ascii":
        start_str = "solid"
        end_str = "endsolid"
        f = open(filename, "r")
        list_of_lines = f.readlines()
        f.close()
        list_of_lines[0] = "{0} {1}\n".format(start_str, solid_name)
        list_of_lines[-1] = "{0} {1}\n".format(end_str, solid_name)

        f = open(filename, "w")
        f.writelines(list_of_lines)
        f.close()
    # replace part name in binary file
    elif file_type == "binary":
        fid = open(filename, "r+b")
        fid.seek(0)
        data = fid.read(40)
        fid.seek(0)
        string_replace = "{:<40}".format(solid_name).encode()
        fid.write(string_replace)
        fid.close()

    return


if __name__ == "__main__":
    logger.info("Protected")
