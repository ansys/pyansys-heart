from ansys.heart.postprocessor.auto_process import mech_post, zerop_post
from ansys.heart.preprocessor.models import HeartModel

if __name__ == "__main__":
    """
    Example show how to post-process mechanical simulation results
    """
    # Get heart model
    model = HeartModel.load_model(".")
    # Define directory where zeropressure simulation is carried out
    # a folder "post" will be created under it with key simulation results (json, png, vtk...)
    # use Paraview state file post_zerop2.pvsm, you can easily visualize generated results
    dir = r""
    zerop_post(dir, model)

    # Define directory where main mechanical simulation is carried out
    # a folder "post" will be created under it with key simulation results (json, png, vtk...)
    # use Paraview state file post_main2.pvsm, you can easily visualize generated results
    dir2 = r""
    mech_post(dir2, model)
