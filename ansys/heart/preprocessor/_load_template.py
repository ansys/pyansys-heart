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

#
# Methods to configure jinja2 and load templates.
#
# import preprocessing
import logging
import os
from pathlib import Path
import posixpath

LOGGER = logging.getLogger("pyheart_global.preprocessor")
import jinja2

# Shared Jinja environment
_environment = None


def get_this_folder():
    """Get folder of this file."""
    return Path(__file__).parent


def get_loader():
    """Get the loader."""
    template_folder = get_this_folder() / "templates"
    return jinja2.FileSystemLoader(str(template_folder.resolve()))


def _jinja_environment():
    """Return a shared Jinja environment to create templates from."""
    global _environment
    if _environment is None:
        _environment = jinja2.Environment(
            # Automatic loading of templates stored in the module
            # This also enables template inheritance
            loader=get_loader(),
            # Keep a single trailing newline, if present
            keep_trailing_newline=True,
            # Don't replace undefined template variables by an empty string
            # but raise a jinja2.UndefinedError instead.
            undefined=jinja2.StrictUndefined,
        )
    return _environment


def load_template(*name):
    """
    Load a template from the local template directory.

    Templates can be specified as a single filename, e.g.
    ``load_template('temp.txt')``, or loaded from subdirectories using e.g.
    ``load_template('subdir_1', 'subdir_2', 'file.txt')``.

    """
    # Due to a Jinja2 convention, posixpaths must be used, regardless of the
    # user's operating system!
    path = posixpath.join(*name)
    if os.path.sep != "/" and os.path.sep in path:  # pragma: no linux cover
        LOGGER.warning("Paths to templates must be specified as posix paths.")

    env = _jinja_environment()
    return env.get_template(path)


if __name__ == "__main__":
    LOGGER.info("Protected")
