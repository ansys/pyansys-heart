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

"""PyAnsys Heart is a Python framework for heart modeling using Ansys tools."""

import os

try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:  # pragma: no cover
    import importlib_metadata

from ansys.health.heart.logger import LOG_LEVEL_FILE, LOG_LEVEL_STDOUT, Logger

LOG = Logger(LOG_LEVEL_STDOUT, to_file=False, to_stdout=True)
LOG.log_to_file(os.path.join(os.getcwd(), "PyAnsys-Heart.log"), LOG_LEVEL_FILE)
LOG.debug("Loaded logging module as LOG")

try:
    import pyvista as pv

    pv.OFF_SCREEN = bool(os.environ["PYVISTA_OFF_SCREEN"])
    LOG.debug(f"Pyvista OFF_SCREEN: {pv.OFF_SCREEN}")
except KeyError:
    pass

__version__ = importlib_metadata.version("ansys-health-heart")
"""Version of PyAnsys Heart."""
