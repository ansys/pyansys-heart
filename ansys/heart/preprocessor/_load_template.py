#
# Methods to configure jinja2 and load templates.
#
import os
from pathlib import Path
import posixpath

# import preprocessing
from ansys.heart.custom_logging import LOGGER
import jinja2

# Shared Jinja environment
_environment = None


def get_this_folder():
    return Path(__file__).parent


def get_loader():
    template_folder = get_this_folder() / "templates"
    return jinja2.FileSystemLoader(str(template_folder.resolve()))


def _jinja_environment():
    """
    Returns a shared Jinja environment to create templates from.
    """
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
    Loads a template from the local template directory.

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
