import pytest
import pyvista as pv
from ansys.heart.preprocessor.mesh.vtkmethods import are_connected


def test_check_if_connected():
    """Test whether two pyvista objects are connected."""
    box1 = pv.Box()
    box2 = pv.Box()

    assert are_connected(box1, box2.translate((2, 0, 0))) == True
    assert are_connected(box1, box2.translate((3, 0, 0))) == False

    return
