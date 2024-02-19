"""PyAnsys Heart is a Python framework for heart modeling using ANSYS tools."""

try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:  # pragma: no cover
    import importlib_metadata

import logging

from ansys.heart.core.logger import Logger

__version__ = importlib_metadata.version("pyansys-heart")
"""The version of pyansys-heart."""

LOG = Logger(level=logging.DEBUG, to_file=False, to_stdout=True)
