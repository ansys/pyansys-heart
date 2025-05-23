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

"""

Postprocess a Reaction-Eikonal model.
-------------------------------------
This example shows how to postprocess a full heart reaction eikonal model.
"""

###############################################################################
# .. warning::
#    When using a standalone version of the DPF Server, you must accept the `license terms
#    <https://dpf.docs.pyansys.com/version/stable/getting_started/licensing.html>`_. To
#    accept these terms, you can set this environment variable:
#
#    .. code-block:: python
#
#        import os
#        os.environ["ANSYS_DPF_ACCEPT_LA"] = "Y"

###############################################################################
# Perform the required imports
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import the required modules and set relevant paths.

from pathlib import Path

from ansys.health.heart.post.dpf_utils import EPpostprocessor

###############################################################################
# Create a postprocessor object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###############################################################################
# .. note::
#    This example assumes that you have you ran a full heart electrophysiology simulation
#    and that the d3plot files are located in ``data_path``.

# Import the required modules and set relevant paths.
workdir = Path.home() / "pyansys-heart" / "downloads" / "Rodero2021" / "01" / "FullHeart"

# Specify the path to the d3plot that contains the simulation results.
data_path = workdir / "simulation-EP" / "main-ep-ReactionEikonal" / "d3plot"

# Check if the file exists.
if not data_path.is_file():
    raise FileNotFoundError(f"File not found: {data_path}")

# Initialize the postprocessor.
post = EPpostprocessor(data_path)

###############################################################################
# Call methods to retrieve activation time
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get activation time of the full field at the last time step.
activation_times = post.get_activation_times()
print(activation_times.data)

activation_times.plot(show_edges=False, show_scalar_bar=True)

###############################################################################
# Create a clip view.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create a clip view of the activation time using ``pyvista``.
import pyvista as pv

# Retrieve the unstructured grid.
grid: pv.UnstructuredGrid = post.reader.model.metadata.meshed_region.grid
grid.point_data["activation_time"] = activation_times.data
grid.set_active_scalars("activation_time")

# Clip the model and plot.
grid.clip(
    normal=[0.7785200198880087, -0.027403237199259987, 0.6270212446357586],
    origin=[88.24004990770091, 54.41149629465821, 49.1801566480857],
).plot(show_scalar_bar=True)
