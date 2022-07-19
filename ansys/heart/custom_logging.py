"""Logging module.

This module supplies a general framework for logging in pyheartlib.
This is adopted and modified from PyFluent
"""

import logging
import os
from pathlib import Path
import tempfile
from typing import Any


class Logger:
    """Logger class.

    Methods
    -------
    set_level(level)
        Set logging level
    enable_logging_to_stdout()
        Enable logging to stdout
    disable_logging_to_stdout()
        Disable logging to stdout
    enable_logging_to_file(filepath)
        Enable logging to file
    disable_logging_to_file()
        Disable logging to file
    """

    def __init__(self, level: Any = logging.DEBUG):
        self.logger = logging.getLogger("pyheart")
        self.logger.propagate = False
        self.stream_handler = None
        self.file_handler = None
        self.log_filepath = None
        dateformat = "%Y/%m/%d %H:%M:%S"
        self.formatter = logging.Formatter(
            "%(asctime)s - %(levelname)s - %(message)s", datefmt=dateformat
        )
        self.enable_logging_to_stdout()
        self.logger.setLevel(level)

        # Writing logging methods.
        self.debug = self.logger.debug
        self.info = self.logger.info
        self.warning = self.logger.warning
        self.error = self.logger.error
        self.critical = self.logger.critical
        self.log = self.logger.log

    def set_level(self, level: Any) -> None:
        """Set logging level.

        Parameters
        ----------
        level : Any
            Any of the logging level (CRITICAL, ERROR, WARNING, INFO, DEBUG)
            in string or enum format
        """
        self.logger.setLevel(level)

    def _get_default_log_filepath(self):
        fd, filepath = tempfile.mkstemp(
            suffix=f"-{os.getpid()}.txt",
            prefix="pyfluent-",
            dir=str(Path.cwd()),
        )
        os.close(fd)
        return Path(filepath)

    def enable_logging_to_stdout(self) -> None:
        """Enable logging to stdout."""
        if self.stream_handler is None:
            self.stream_handler = logging.StreamHandler()
            self.stream_handler.setFormatter(self.formatter)
        self.logger.addHandler(self.stream_handler)

    def disable_logging_to_stdout(self) -> None:
        """Disable logging to stdout."""
        self.logger.removeHandler(self.stream_handler)

    def enable_logging_to_file(self, filepath: str = None) -> None:
        """Enable logging to file.

        Parameters
        ----------
        filepath : str, optional
            filapath, a default filepath will be chosen if filepath is not
            passed
        """
        self.logger.removeHandler(self.file_handler)
        if not filepath and not self.log_filepath:
            self.log_filepath = self._get_default_log_filepath()
        elif filepath:
            self.log_filepath = Path(filepath)
        self.file_handler = logging.FileHandler(self.log_filepath)
        self.file_handler.setFormatter(self.formatter)
        self.logger.addHandler(self.file_handler)

    def disable_logging_to_file(self) -> None:
        """Disable logging to file."""
        self.logger.removeHandler(self.file_handler)


LOGGER = Logger()
LOGGER.enable_logging_to_stdout()

if __name__ == "__main__":
    print()
    # log = CustomLogger( to_stdout = True )
