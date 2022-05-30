from pathlib import Path
import os, subprocess
import numpy as np
import gmsh


from ansys.heart.preprocessor._load_template import load_template
from preprocessing import SC_EXE, FLUENT_EXE
from preprocessing.custom_logging import logger

_template_directory = os.path.join(os.path.dirname(__file__), "template")

"""Module contains methods for mesh operations"""

## for remeshing purposes


def mesh_by_fluentmeshing(
    path_to_input_stl: str, path_to_output: str, mesh_size: float = 2.0
):
    """Uses Fluent meshing to wrap the surface and create
    tetrahedral mesh"""

    var_for_template = {
        "input_path": path_to_input_stl,
        "output_path": path_to_output,
        "mesh_size": mesh_size,
    }

    template = load_template("fluent_meshing_template.jou")

    script = os.path.join(
        os.path.dirname(path_to_input_stl), "fluent_meshing.jou"
    )
    # script = "fluent_meshing.jou"

    with open(script, "w") as f:
        f.write(template.render(var_for_template))

    # subprocess.call( [SC_EXE, "/RunScript=" + script, *options] )
    num_cpus = 2

    # args = ['"' + FLUENT_EXE + '"', "-v3ddp", "-tm{:.0f}".format(num_cpus), "-meshing", "-i", script  ]
    args = [
        '"' + FLUENT_EXE + '"',
        "-v3ddp",
        "-hidden",
        "-tm{:.0f}".format(num_cpus),
        "-meshing",
        "-i",
        script,
    ]

    # subprocess.call( " ".join(args) )
    # p = subprocess.Popen(" ".join(args), shell=True)

    # TODO: need to add check on output
    logger.info("Launching Fluent Meshing...")
    p = subprocess.Popen(" ".join(args))
    p.communicate()
    p.wait(300)
    logger.info("done")

    # p = subprocess.call(" ".join(args) )

    return


def shrink_by_spaceclaim(input, output):
    """Uses SpaceClaim shrinkwrapping to wrap surface and 
    create high quality surface mesh"""

    var_for_template = {
        "input": input,
        "output": output,
    }
    posixpath_template = str(
        os.path.join(_template_directory, "spaceclaim_shrink")
    )
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


def add_solid_name_to_stl(filename, solid_name):
    """Adds name of solid to stl file. Supports only single block!"""
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

    return


if __name__ == "__main__":
    logger.info("Protected")
