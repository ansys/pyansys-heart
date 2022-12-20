import numpy as np
import pyvista as pv
from pyvista import CellType
import ansys.heart.preprocessor.models as models

# mesh points
vertices = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 0.5, -1]])
model: models.FullHeart = models.HeartModel.load_model(
    r"D:\development\pyheart-lib\pyheart-lib\downloads\Strocchi2020\01\FullHeart\heart_model.pickle"
)

CellType.TETRA

# visualize with pyvista
num_faces = model.left_atrium.endocardium.triangles.shape[0]
a = np.ones((num_faces, 1), dtype=int) * 3
faces = np.hstack([a, model.left_atrium.endocardium.triangles])
nodes = model.left_atrium.endocardium.nodes
faces = np.reshape(faces, (faces.size))

surf = pv.PolyData(nodes, faces)

mesh = surf
edges = mesh.extract_all_edges()

# pl = pv.Plotter()
# pl.add_mesh(mesh)
# act = pl.add_mesh(edges, line_width=2, color='k')
# pl.show()

# plot each face with a different color
surf.plot(
    # scalars=np.arange(3),
    cpos=[-1, 1, 0.5],
    show_scalar_bar=False,
    show_edges=True,
    line_width=2,
)
