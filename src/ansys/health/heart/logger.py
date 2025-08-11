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

"""Logging module.

The logging module supplies a general framework for logging in PyAnsys Heart.
This module is built on the Python `logging <https://docs.python.org/3/library/logging.html>`_
library. It does not intend to replace it but rather provide a way to interact between
the Python ``logging`` library and PyAnsys Heart.

The loggers used in the module include the name of the instance, which
is intended to be unique. This name is printed in all the active
outputs and is used to track the different PyAnsys Heart modules.


Usage
-----

Global logger
~~~~~~~~~~~~~
There is a global logger named ``PyAnsys Heart_global`` that is created when
``ansys.health.heart.__init__`` is called.  If you want to use this global logger,
you must call it at the top of your module:

.. code:: python

   from ansys.health.heart import LOG

You can rename this logger to avoid conflicts with other loggers (if any):

.. code:: python

   from ansys.health.heart import LOG as logger


The default logging level of ``LOG`` is ``ERROR``.
You can change this level and output lower-level messages:

.. code:: python

   LOG.logger.setLevel("DEBUG")
   LOG.file_handler.setLevel("DEBUG")  # If present.
   LOG.std_out_handler.setLevel("DEBUG")  # If present.


Alternatively, you can ensure all the handlers are set to the input log level
with this code:

.. code:: python

   LOG.setLevel("DEBUG")

This logger does not log to a file by default. If you want, you can
add a file handler with this code:

.. code:: python

   import os

   file_path = os.path.join(os.getcwd(), "pymapdl.log")
   LOG.log_to_file(file_path)

This also sets the logger to be redirected to this file. If you want
to change the characteristics of this global logger from the beginning
of the execution, you must edit the file ``__init__`` file in the
``ansys.health.heart`` directory.

To log using this logger, call the desired method as a normal logger:

.. code:: pycon

    >>> import logging
    >>> from ansys.health.heart.logging import Logger
    >>> LOG = Logger(level=logging.DEBUG, to_file=False, to_stdout=True)
    >>> LOG.debug("This is LOG debug message.")

    DEBUG -  -  <ipython-input-24-80df150fe31f> - <module> - This is LOG debug message.

Other loggers
~~~~~~~~~~~~~
You can create your own loggers using the Python ``logging`` library as
you would do in any other script.  There would be no conflicts between
these loggers.

"""

from copy import copy
from datetime import datetime
import logging
import sys
from types import TracebackType
from typing import Any, Dict, Literal, Mapping, MutableMapping, Optional, Type, Union, cast

## Default configuration
LOG_LEVEL_STDOUT = logging.INFO
LOG_LEVEL_FILE = logging.DEBUG
FILE_NAME = "PyAnsys Heart.log"

# For convenience
DEBUG = logging.DEBUG
INFO = logging.INFO
WARN = logging.WARN
ERROR = logging.ERROR
CRITICAL = logging.CRITICAL

## Formatting

STDOUT_MSG_FORMAT = (
    "%(asctime)s - %(levelname)s - %(instance_name)s - %(module)s - %(funcName)s - %(message)s"
)

DATEFORMAT = "%Y/%m/%d %H:%M:%S"

FILE_MSG_FORMAT = STDOUT_MSG_FORMAT

DEFAULT_STDOUT_HEADER = """
LEVEL - INSTANCE NAME - MODULE - FUNCTION - MESSAGE
"""
DEFAULT_FILE_HEADER = DEFAULT_STDOUT_HEADER

NEW_SESSION_HEADER = f"""
===============================================================================
       NEW SESSION - {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}
==============================================================================="""

LOG_LEVEL_STRING_TYPE = Literal["DEBUG", "INFO", "WARN", "WARNING", "ERROR", "CRITICAL"]
LOG_LEVEL_TYPE = Union[LOG_LEVEL_STRING_TYPE, int]

string_to_loglevel: Dict[LOG_LEVEL_STRING_TYPE, int] = {
    "DEBUG": DEBUG,
    "INFO": INFO,
    "WARN": WARN,
    "WARNING": WARN,
    "ERROR": ERROR,
    "CRITICAL": CRITICAL,
}


class PyAnsysHeartCustomAdapter(logging.LoggerAdapter):
    """Keeps the reference to the PyAnsys Heart service instance dynamic.

    If you use the standard approach, which is supplying *extra* input
    to the logger, you must input PyAnsys Heart service instances
    each time that you log a message.

    Using adapters, you only need to specify the PyAnsys Heart service instance
    that you are referring to once.
    """

    level = (
        None  # This is maintained for compatibility with ``supress_logging``, but it does nothing.
    )
    file_handler: Optional[logging.FileHandler] = None
    std_out_handler: Optional[logging.StreamHandler] = None

    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self.file_handler = logger.file_handler
        self.std_out_handler = logger.std_out_handler

    def process(self, msg: str, kwargs: MutableMapping[str, Dict[str, str]]):
        """Process extra arguments."""
        kwargs["extra"] = {}
        # This are the extra parameters sent to log
        kwargs["extra"]["instance_name"] = (
            self.extra.name
        )  # here self.extra is the argument pass to the log records.
        return msg, kwargs

    def log_to_file(
        self, filename: str = FILE_NAME, level: LOG_LEVEL_TYPE = LOG_LEVEL_FILE
    ) -> None:
        """Add a file handler to the logger.

        Parameters
        ----------
        filename : str, default: FILE_NAME
            Name of the file to record logs to.
        level : str or int, default: LOG_LEVEL
            Level of logging, such as ``DEBUG``.
        """
        addfile_handler(self.logger, filename=filename, level=level, write_headers=True)
        self.file_handler = self.logger.file_handler

    def log_to_stdout(self, level: LOG_LEVEL_TYPE = LOG_LEVEL_STDOUT) -> None:
        """Add a stdout handler to the logger.

        Parameters
        ----------
        level : str or int, default: LOG_LEVEL
            Level of the logging record.
        """
        if self.std_out_handler:
            raise Exception("Stdout logger is already defined.")

        add_stdout_handler(self.logger, level=level)
        self.std_out_handler = self.logger.std_out_handler

    def setLevel(self, level: Union[int, str] = "DEBUG"):  # noqa: N802
        """Change the log level of the object and the attached handlers."""
        if isinstance(level, str):
            level = string_to_loglevel[cast(LOG_LEVEL_STRING_TYPE, level.upper())]
        self.logger.setLevel(level)
        for each_handler in self.logger.handlers:
            each_handler.setLevel(level)
        self.level = level


class _PyAnsysHeartPercentStyle(logging.PercentStyle):
    def __init__(self, fmt, *, defaults=None):
        self._fmt = fmt or self.default_format
        self._defaults = defaults

    def _format(self, record) -> str:
        defaults = self._defaults
        if defaults:
            values = defaults | record.__dict__
        else:
            values = record.__dict__

        # We can make any changes that we want in the record here. For example, adding a key.

        # We could create an ``if`` here if we want conditional formatting, and even
        # change the record.__dict__.
        # Because we don't want to create conditional fields now, it is fine to keep
        # the same MSG_FORMAT for all of them.

        # For the case of logging exceptions to the logger.
        values.setdefault("instance_name", "")

        return STDOUT_MSG_FORMAT % values


class _PyAnsysHeartFormatter(logging.Formatter):
    """Provides a ``Formatter`` class for overwriting default format styles."""

    def __init__(
        self,
        fmt: str = STDOUT_MSG_FORMAT,
        datefmt: Optional[str] = DATEFORMAT,
        style: Literal["%", "{", "$"] = "%",
        validate: bool = True,
        defaults: Optional[Mapping[str, Any]] = None,
    ):
        if sys.version_info[1] < 8:
            super().__init__(fmt, datefmt, style)
        else:
            # 3.8: The validate parameter was added
            super().__init__(fmt, datefmt, style, validate)
        self._style = _PyAnsysHeartPercentStyle(fmt, defaults=defaults)  # overwriting


class InstanceFilter(logging.Filter):
    """Ensures that the ``instance_name`` record always exists."""

    def filter(self, record: logging.LogRecord):
        """Filter record."""
        if not hasattr(record, "instance_name") and hasattr(record, "name"):
            record.instance_name = record.name
        elif not hasattr(record, "instance_name"):  # pragma: no cover
            record.instance_name = ""
        return True


class Logger:
    """Provides the logger used for each PyAnsys Heart session.

    This class lets you add handlers to the logger to output messages to a file or
    to the standard output (stdout).

    Parameters
    ----------
    level : int, default: logging.DEBUG
        Logging level to filter the message severity allowed in the logger.
    to_file : bool, default: False
        Whether to write log messages to a file.
    to_stdout : bool, default: True
        Whether to write the log messages to stdout.
    filename : str, default: FILE_NAME
        Name of the file to write log messages to.

    Examples
    --------
    Demonstrate logger usage from a PyAnsys Heart instance, which is automatically
    created when a PyAnsys Heart instance is created.

    Import the global PyAnsys Heart logger and add a file output handler.

    >>> import os
    >>> from ansys.health.heart import LOG
    >>> file_path = os.path.join(os.getcwd(), "PyAnsys Heart.log")
    >>> LOG.log_to_file(file_path)
    """

    file_handler: Optional[logging.FileHandler] = None
    std_out_handler: Optional[logging.StreamHandler] = None
    _level = logging.DEBUG
    _instances: Dict[str, Any] = {}

    def __init__(
        self,
        level: LOG_LEVEL_TYPE = logging.DEBUG,
        to_file: bool = False,
        to_stdout: bool = True,
        filename: str = FILE_NAME,
    ):
        """Initialize the main logger class for PyAnsys Heart.

        Parameters
        ----------
        level : str or int, default: logging.DEBUG
            Level of logging as defined in the ``logging`` package.
        to_file : bool, default: False
            Whether to write log messages to a file.
        to_stdout : bool, default: True
            Whether to write log messages to the standard output (stdout), which is the
            command line.
        filename : str, default: FILE_NAME
            Name of the output file, which is ``'PyAnsys Heart.log'`` by default.
        """
        # create default main logger
        self.logger: logging.Logger = logging.getLogger("PyAnsys Heart_global")
        self.logger.addFilter(InstanceFilter())
        if isinstance(level, str):
            level = cast(LOG_LEVEL_STRING_TYPE, level.upper())

        self.logger.setLevel(level)
        self.logger.propagate = True
        self.level = self.logger.level  # TODO: TO REMOVE

        # Writing logging methods.
        self.debug = self.logger.debug
        self.info = self.logger.info
        self.warning = self.logger.warning
        self.error = self.logger.error
        self.critical = self.logger.critical
        self.log = self.logger.log

        if to_file:
            # We record to file
            self.log_to_file(filename=filename, level=level)

        if to_stdout:
            self.log_to_stdout(level=level)

        # Using logger to record unhandled exceptions
        self.add_handling_uncaught_expections(self.logger)

    def log_to_file(
        self,
        filename: str = FILE_NAME,
        level: LOG_LEVEL_TYPE = LOG_LEVEL_FILE,
        remove_other_file_handlers: bool = False,
    ) -> None:
        """Add a file handler to the logger.

        Parameters
        ----------
        filename : str, default:
            Name of the file to record logs to, which is ``'PyAnsys Heart.log'`` by default.
        level : str or int, default: LOG_LEVEL_FILE
            Level of logging, which is ``'DEBUG'`` by default.
        remove_other_file_handlers : bool, default: False
            Whether to remove all other file handlers.

        Examples
        --------
        Write to the ``PyAnsys Heart.log`` file in the current working directory.

        >>> from ansys.health.heart import LOG
        >>> import os
        >>> file_path = os.path.join(os.getcwd(), "PyAnsys Heart.log")
        >>> LOG.log_to_file(file_path)
        """
        if remove_other_file_handlers:
            _clear_all_file_handlers(self)

        addfile_handler(self, filename=filename, level=level, write_headers=True)

    def log_to_stdout(self, level: LOG_LEVEL_TYPE = LOG_LEVEL_STDOUT):
        """Add a stdout handler to the logger.

        Parameters
        ----------
        level : str or int, default: LOG_LEVEL_STDOUT
            Level of logging record, which is ``'DEBUG'`` by default.
        """
        add_stdout_handler(self, level=level)

    def setLevel(self, level: LOG_LEVEL_TYPE = "DEBUG"):  # noqa: N802
        """Set the log level for the logger and its handlers.

        Parameters
        ----------
        level : str or int, default: "DEBUG"
            Logging level to set.
        """
        if isinstance(level, str):
            level = string_to_loglevel[cast(LOG_LEVEL_STRING_TYPE, level.upper())]
        self.logger.setLevel(level)
        for each_handler in self.logger.handlers:
            each_handler.setLevel(level)
        self._level = level

    def _make_child_logger(self, suffix: str, level: Optional[LOG_LEVEL_TYPE]) -> logging.Logger:
        """Create a child logger.

        This method uses the ``getChild()``method or copies attributes between the
        ``pymapdl_global`` logger and the new one.
        """
        logger = logging.getLogger(suffix)
        logger.std_out_handler = None
        logger.file_handler = None

        if self.logger.hasHandlers():
            for each_handler in self.logger.handlers:
                new_handler = copy(each_handler)
                if each_handler == self.file_handler:
                    logger.file_handler = new_handler
                elif each_handler == self.std_out_handler:
                    logger.std_out_handler = new_handler

                if level:
                    # The logger handlers are copied and changed. The loglevel is
                    # the specified log level is lower than the one of the
                    # global.
                    if isinstance(level, str):
                        new_loglevel = string_to_loglevel[
                            cast(LOG_LEVEL_STRING_TYPE, level.upper())
                        ]
                    elif isinstance(level, int):  # pragma: no cover
                        new_loglevel = level

                    if each_handler.level > new_loglevel:
                        new_handler.setLevel(level)

                logger.addHandler(new_handler)

        if level:
            if isinstance(level, str):
                level = string_to_loglevel[cast(LOG_LEVEL_STRING_TYPE, level.upper())]
            logger.setLevel(level)

        else:
            logger.setLevel(self.logger.level)

        logger.propagate = True
        return logger

    def add_child_logger(self, suffix: str, level: Optional[LOG_LEVEL_TYPE] = None):
        """Add a child logger to the main logger.

        This logger is more general than an instance logger, which is designed to
        track the state of PyAnsys Heart instances.

        If the logging level is in the arguments, a new logger with a reference
        to the ``_global`` logger handlers is created instead of a child logger.

        Parameters
        ----------
        suffix : str
            Name of the logger.
        level : str or int, default: None
            Level of logging.

        Returns
        -------
        logging.logger
            Logger class.
        """
        name = self.logger.name + "." + suffix
        self._instances[name] = self._make_child_logger(name, level)
        return self._instances[name]

    def __getitem__(self, key):
        """Overload the access method by item for the ``Logger`` class."""
        if key in self._instances.keys():
            return self._instances[key]
        else:
            raise KeyError(f"There is no instances with name {key}.")

    def add_handling_uncaught_expections(self, logger: logging.Logger):
        """Redirect the output of an exception to a logger.

        Parameters
        ----------
        logger : str
            Name of the logger.
        """

        def handle_exception(
            exc_type: Type[BaseException],
            exc_value: BaseException,
            exc_traceback: Optional[TracebackType],
        ):
            if issubclass(exc_type, KeyboardInterrupt):
                sys.__excepthook__(exc_type, exc_value, exc_traceback)
                return
            logger.critical(
                "Uncaught exception",
                exc_info=(exc_type, exc_value, exc_traceback),
            )

        sys.excepthook = handle_exception


def addfile_handler(logger, filename=FILE_NAME, level=LOG_LEVEL_STDOUT, write_headers=False):
    """
    Add a file handler to the input.

    Parameters
    ----------
    logger : logging.Logger
        Logger to add the file handler to.
    filename : str, default: FILE_NAME
        Name of the output file, which is ``'pyconv-de.log'`` by default.
    level : int, default: 10
        Level of logging. The default is ``10``, in which case the
        ``logging.DEBUG`` level is used.
    write_headers : bool, default: False
        Whether to write headers to the file.

    Returns
    -------
    Logger
        :class:`Logger` or :class:`logging.Logger` object.
    """
    file_handler = logging.FileHandler(filename)
    file_handler.setLevel(level)
    file_handler.setFormatter(logging.Formatter(FILE_MSG_FORMAT))

    if isinstance(logger, Logger):
        logger.file_handler = file_handler
        logger.logger.addHandler(file_handler)

    elif isinstance(logger, logging.Logger):
        logger.file_handler = file_handler
        logger.addHandler(file_handler)

    if write_headers:
        file_handler.stream.write(NEW_SESSION_HEADER)
        file_handler.stream.write(DEFAULT_FILE_HEADER)

    return logger


def _clear_all_file_handlers(logger: Logger) -> Logger:
    """Clear all file handlers from the logger.

    Parameters
    ----------
    logger : Logger
        Logger to clear file handlers from.

    Returns
    -------
    Logger
        Logger without file handlers.
    """
    file_handlers = [
        handler for handler in logger.logger.handlers if isinstance(handler, logging.FileHandler)
    ]
    for handler in file_handlers:
        logger.logger.removeHandler(handler)
        handler.close()
    return logger


def add_stdout_handler(logger, level=LOG_LEVEL_STDOUT, write_headers=False):
    """
    Add a stdout handler to the logger.

    Parameters
    ----------
    logger : logging.Logger
        Logger to add the stdout handler to.
    level : int, default: ``10``
        Level of logging. The default is ``10``, in which case the
        ``logging.DEBUG`` level is used.
    write_headers : bool, default: False
        Whether to write headers to the file.

    Returns
    -------
    Logger
        :class:`Logger` or :class:`logging.Logger` object.
    """
    std_out_handler = logging.StreamHandler()
    std_out_handler.setLevel(level)
    std_out_handler.setFormatter(_PyAnsysHeartFormatter(STDOUT_MSG_FORMAT))

    if isinstance(logger, Logger):
        logger.std_out_handler = std_out_handler
        logger.logger.addHandler(std_out_handler)

    elif isinstance(logger, logging.Logger):
        logger.addHandler(std_out_handler)

    if write_headers:
        std_out_handler.stream.write(DEFAULT_STDOUT_HEADER)

    return logger


# ===============================================================
# Finally define logger
# ===============================================================

# LOG = Logger(level=logging.INFO, to_file=False, to_stdout=True)
# LOG.debug("Loaded logging module as LOG")
