"""Collection of modules for heart modeling."""

try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:  # pragma: no cover
    import importlib_metadata

# Read from the pyproject.toml
# major, minor, patch

import logging

from ansys.heart.core.logging import Logger

__version__ = importlib_metadata.version("pyansys-heart")

LOG = Logger(level=logging.DEBUG, to_file=False, to_stdout=True)
