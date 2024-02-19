"""Version of pyansys-heart library.

On the ``main`` branch, use 'dev0' to denote a development version.
For example:

version_info = 0, 1, 'dev0'

Examples
--------
Print the version

>>> from ansys.heart import _version
>>> print(library.__version__)
0.1.dev0

"""

try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:  # pragma: no cover
    import importlib_metadata

# Read from the pyproject.toml
# major, minor, patch
__version__ = importlib_metadata.version("pyansys-heart")
