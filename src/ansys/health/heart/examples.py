# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
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

"""Utility functions to get various examples."""

from pathlib import Path
from typing import Literal

import pyvista as pv

dir_path = Path(__file__).parent
data_path = dir_path / "data_examples"


def get_preprocessed_fullheart(
    resolution: Literal["1.5mm", "2.0mm"] = "2.0mm",
) -> tuple[Path, Path, Path]:
    """Get a preprocessed full heart model.

    Parameters
    ----------
    resolution : Literal[&quot;1.5mm&quot;, &quot;2.0mm&quot;], default: "2.0mm"
        The resolution of the preprocessed full heart.

    Returns
    -------
    tuple[Path, Path, Path]
        Paths to the ``.vtu``, ``.partinfo.json``, and ``.namemap.json`` files.
    """
    return (
        str(data_path / f"rodero_01_fullheart_{resolution}.vtu"),
        str(data_path / f"rodero_01_fullheart_{resolution}.partinfo.json"),
        str(data_path / f"rodero_01_fullheart_{resolution}.namemap.json"),
    )


def get_input_leftventricle() -> tuple[Path, Path]:
    """Get the input of a left ventricle based on Rodero et al 01."""
    return (
        str(data_path / "rodero_01_leftventricle_surface.vtp"),
        str(data_path / "rodero_01_leftventricle_part_definition.json"),
    )


def get_fractal_tree_purkinje() -> tuple[pv.PolyData, pv.PolyData]:
    """Get the fractal tree Purkinje network based on Rodero et al 01.

    Returns
    -------
    tuple[pv.PolyData, pv.PolyData]
        The left and right Purkinje networks as PyVista PolyData objects.
    """
    return (
        pv.read(data_path / "left_purkinje.vtp"),
        pv.read(data_path / "right_purkinje.vtp"),
    )
