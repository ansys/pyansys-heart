import os
import tempfile

import numpy as np

from ansys.heart.core.helpers.element_orth import (
    _read_orth_element_kfile,
)

_test_data = """*ELEMENT_SOLID_ORTHO
       1       1
       1       2       3       4       4       4       4       4
 -0.79687124E+00  0.60410495E+00  0.73101645E-02
 -0.21393461E+00 -0.27219762E+00 -0.93873706E+00
       2       1
       5       6       7       8       8       8       8       8
  0.43854175E+00  0.29742697E+00 -0.84806741E+00
 -0.59030128E+00 -0.61108034E+00 -0.52510346E+00
       3       2
       7       6       9       8       8       8       8       8
  0.46615814E+00  0.31190270E+00 -0.82789691E+00
 -0.61661678E+00 -0.57149497E+00 -0.54729073E+00
       4       3
      10      11      12      13      13      13      13      13
 -0.89163671E+00  0.43980761E+00 -0.10748604E+00
 -0.18581314E+00 -0.57499262E+00 -0.79806465E+00
*END
"""


def test_read_orth_element_kfile():
    with tempfile.TemporaryDirectory(prefix=".pyansys-heart") as tempdir:
        file = os.path.join(tempdir, "element_solid_ortho.k")
        with open(file, "w") as f:
            f.write(_test_data)

        el_ids, part_ids, conn, _, _ = _read_orth_element_kfile(file)

        assert np.allclose(el_ids, [1, 2, 3, 4])
        assert np.allclose(part_ids, [1, 1, 2, 3])
        assert np.allclose(
            conn[:, 0:4], np.array([[1, 2, 3, 4], [5, 6, 7, 8], [7, 6, 9, 8], [10, 11, 12, 13]])
        )
