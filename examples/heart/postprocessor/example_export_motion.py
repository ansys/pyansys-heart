# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""For information only, not yet stable
    require previously launched simulation
    change paths accordingly"""

import matplotlib.pyplot as plt

from ansys.heart.postprocessor.exporter import LVContourExporter
from ansys.heart.preprocessor.models import HeartModel

if __name__ == "__main__":
    """
    Show exporting LV motion in vtk
    require previously launched simulation
    change paths accordingly
    """
    model: HeartModel
    model = HeartModel.load_model(r"path_to_model.pickle")
    model.compute_left_ventricle_anatomy_axis(first_cut_short_axis=0.2)

    d3plot_file = r"my_path_to_d3plot"

    exporter = LVContourExporter(d3plot_file, model)

    exporter.export_contour_to_vtk("l4cv", model.l4cv_axis)

    exporter.export_contour_to_vtk("l2cv", model.l2cv_axis)

    normal = model.short_axis["normal"]
    p_start = model.short_axis["center"]
    for ap in model.left_ventricle.apex_points:
        if ap.name == "apex epicardium":
            p_end = ap.xyz

    for icut in range(2):
        p_cut = p_start + (p_end - p_start) * icut / 2
        cutter = {"center": p_cut, "normal": normal}
        exporter.export_contour_to_vtk(f"shor_{icut}", cutter)

    exporter.export_lvls_to_vtk("lvls")
    exporter.plot_lvls()
    plt.show()
