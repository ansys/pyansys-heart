"""Sphinx documentation configuration file."""
from datetime import datetime

# add sys path for sphinx to find the packages/modules
# NOTE: could be replaced?
import os
import sys

abspath = os.path.dirname(os.path.abspath(__file__))
base_path = os.path.abspath(os.path.join(abspath, "..", ".."))
sys.path.insert(0, base_path)

# dynalib_path = r"D:\development\dynalib\dynalib"
# sys.path.insert(0, dynalib_path)


from ansys.heart._version import __version__
from pyansys_sphinx_theme import pyansys_logo_black

# Project information
project = "ansys-heart-lib"
copyright = f"(c) {datetime.now().year} ANSYS, Inc. All rights reserved"
author = "ANSYS, Inc."
release = version = __version__

# use the default pyansys logo
html_logo = pyansys_logo_black
html_theme = "pyansys_sphinx_theme"

html_short_title = html_title = "ansys-heart-lib"

# specify the location of your github repo
html_theme_options = {
    "github_url": "https://github.com/pyansys/pyheart",
    "show_prev_next": False,
    "show_breadcrumbs": True,
    "additional_breadcrumbs": [
        ("PyAnsys", "https://docs.pyansys.com/"),
    ],
}

# Sphinx extensions
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "numpydoc",
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",
]

autodoc_member_order = "bysource"

# Intersphinx mapping
intersphinx_mapping = {
    "python": ("https://docs.python.org/dev", None),
    # kept here as an example
    # "scipy": ("https://docs.scipy.org/doc/scipy/reference", None),
    # "numpy": ("https://numpy.org/devdocs", None),
    # "matplotlib": ("https://matplotlib.org/stable", None),
    # "pandas": ("https://pandas.pydata.org/pandas-docs/stable", None),
    # "pyvista": ("https://docs.pyvista.org/", None),
    # "grpc": ("https://grpc.github.io/grpc/python/", None),
}

# numpydoc configuration
# numpydoc_show_inherited_class_members = True
numpydoc_show_class_members = False
numpydoc_xref_param_type = True

# Consider enabling numpydoc validation. See:
# https://numpydoc.readthedocs.io/en/latest/validation.html#
numpydoc_validate = True
numpydoc_validation_checks = {
    "GL06",  # Found unknown section
    "GL07",  # Sections are in the wrong order.
    "GL08",  # The object does not have a docstring
    "GL09",  # Deprecation warning should precede extended summary
    "GL10",  # reST directives {directives} must be followed by two colons
    "SS01",  # No summary found
    "SS02",  # Summary does not start with a capital letter
    # "SS03", # Summary does not end with a period
    "SS04",  # Summary contains heading whitespaces
    # "SS05", # Summary must start with infinitive verb, not third person
    "RT02",  # The first line of the Returns section should contain only the
    # type, unless multiple values are being returned"
}


# static path
html_static_path = ["_static"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# # Configuration for Sphinx autoapi
# autoapi_type = "python"
# autoapi_dirs = ["../../ansys/heart"]
# autoapi_options = [
#     "members",
#     "undoc-members",
#     "show-inheritance",
#     "show-module-summary",
#     "special-members",
# ]
# autoapi_template_dir = "_autoapi_templates"
# suppress_warnings = ["autoapi.python_import_resolution"]
# exclude_patterns = ["_autoapi_templates/index.rst"]
# autoapi_python_use_implicit_namespaces = True

# # Examples gallery customization
# nbsphinx_execute = "always"
# nbsphinx_custom_formats = {
#     ".mystnb": ["jupytext.reads", {"fmt": "mystnb"}],
# }
# nbsphinx_thumbnails = {
#     "examples/01_getting_started/01_math": "_static/thumbnails/101_getting_started.png",
#     "examples/01_getting_started/02_units": "_static/thumbnails/101_getting_started.png",
#     "examples/01_getting_started/03_sketching": "_static/thumbnails/101_getting_started.png",
#     "examples/01_getting_started/04_modeling": "_static/thumbnails/101_getting_started.png",
#     "examples/02_sketching/basic_usage": "_static/thumbnails/basic_usage.png",
#     "examples/02_sketching/dynamic_sketch_plane": "_static/thumbnails/dynamic_sketch_plane.png",
#     "examples/02_sketching/advanced_sketching_gears": "_static/thumbnails/advanced_sketching_gears.png",  # noqa: E501
#     "examples/03_modeling/add_design_material": "_static/thumbnails/add_design_material.png",
#     "examples/03_modeling/plate_with_hole": "_static/thumbnails/plate_with_hole.png",
#     "examples/03_modeling/tessellation_usage": "_static/thumbnails/tessellation_usage.png",
#     "examples/03_modeling/design_organization": "_static/thumbnails/design_organization.png",
# }
# nbsphinx_epilog = """
# ----
# .. admonition:: Download this example!
#     Download this example as a `Jupyter Notebook <{cname_pref}/{ipynb_file_loc}>`_
#     or as a `Python script <{cname_pref}/{py_file_loc}>`_ from the previous links.
# """.format(
#     cname_pref=f"https://{cname}/version/{switcher_version}",
#     ipynb_file_loc="{{ env.docname }}.ipynb",
#     py_file_loc="{{ env.docname }}.py",
# )

# nbsphinx_prolog = """
# .. admonition:: Download this example!
#     Download this example as a `Jupyter Notebook <{cname_pref}/{ipynb_file_loc}>`_
#     or as a `Python script <{cname_pref}/{py_file_loc}>`_ from the previous links.
# ----
# """.format(
#     cname_pref=f"https://{cname}/version/{switcher_version}",
#     ipynb_file_loc="{{ env.docname }}.ipynb",
#     py_file_loc="{{ env.docname }}.py",
# )

# typehints_defaults = "comma"
# simplify_optional_unions = False

# # additional logos for the latex coverpage
# latex_additional_files = [watermark, ansys_logo_white, ansys_logo_white_cropped]

# # change the preamble of latex with customized title page
# # variables are the title of pdf, watermark
# latex_elements = {"preamble": latex.generate_preamble(html_title)}
