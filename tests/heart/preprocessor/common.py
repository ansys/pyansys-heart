import pyvista as pv


def _create_simple_unstructured_grid() -> pv.UnstructuredGrid:
    """Create simple unstructured grid."""

    # Prep UnstructuredGrid input.
    cells = [4, 0, 1, 2, 3, 4, 0, 1, 2, 4]
    celltypes = [pv.CellType.TETRA, pv.CellType.TETRA]
    points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]

    ugrid = pv.UnstructuredGrid(cells, celltypes, points)
    ugrid.cell_data.set_scalars([1, 2], "tags")

    return ugrid
