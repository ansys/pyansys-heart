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

"""Custom exceptions for PyAnsys Heart."""


class LSDYNATerminationError(BaseException):
    """Exception raised when ``Normal Termination`` is not found in the LS-DYNA logs."""

    def __init__(self, message: str | list = ""):
        if isinstance(message, list):
            message = "".join(message)
        super().__init__(f"The LS-DYNA process did not terminate as expected: {message}")


class DatabaseNotSupportedError(NotImplementedError):
    """Exception raised when the database is not supported."""

    def __init__(self, db_type, message):
        super().__init__(f"Database type '{db_type}' is not supported. {message}.")


class SupportedDPFServerNotFoundError(Exception):
    """Exception raised when no supported DPF server is found."""


class SupportedFluentVersionNotFoundError(Exception):
    """Exception raised when no supported Fluent version is found."""


class InvalidInputModelTypeError(TypeError):
    """Exception raised when the input heart model type is invalid."""

    def __init__(self, message: str):
        super().__init__(message)


class InvalidHeartModelError(Exception):
    """Exception raised when the heart model is invalid."""


class LSDYNANotFoundError(FileNotFoundError):
    """Exception raised when the LS-DYNA executable file is not found."""


class D3PlotNotSupportedError(IOError):
    """Exception raised when the D3plot file is not supported."""


class MPIProgamNotFoundError(FileNotFoundError):
    """Exception raised when MPI program is not found."""


class WSLNotFoundError(FileNotFoundError):
    """Exception raised when WSL executable is not found."""


class MissingEnvironmentVariableError(EnvironmentError):
    """Exception raised when a required environment variable is missing."""
