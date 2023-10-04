"""Collection of modules for heart modeling."""
import logging

from ansys.heart.core.logging import Logger

LOG = Logger(level=logging.DEBUG, to_file=False, to_stdout=True)
LOG.debug("Loaded logging module as LOG")
