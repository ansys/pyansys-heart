"""Sphinx documentation configuration file."""

from datetime import datetime
import os
from pathlib import Path
import subprocess

from ansys_sphinx_theme import ansys_favicon, get_version_match

from ansys.heart import __version__

# Project information
project = "pyansys-heart"
copyright = f"(c) {datetime.now().year} ANSYS, Inc. All rights reserved"
author = "ANSYS, Inc."
release = version = __version__
cname = os.getenv("DOCUMENTATION_CNAME", "heart.docs.pyansys.com")

# use the default pyansys logo
html_theme = "ansys_sphinx_theme"

html_short_title = html_title = "PyAnsys Heart"
html_favicon = ansys_favicon

# specify the location of your github repo
html_theme_options = {
    "github_url": "https://github.com/ansys/pyansys-heart",
    "show_prev_next": False,
    "show_breadcrumbs": True,
    "additional_breadcrumbs": [
        ("PyAnsys", "https://docs.pyansys.com/"),
    ],
    "switcher": {
        "json_url": f"https://{cname}/versions.json",
        "version_match": get_version_match(__version__),
    },
    "check_switcher": False,
    "logo": "pyansys",
    "ansys_sphinx_theme_autoapi": {
        "project": project,
        "ignore": [
            "*writer*",
            "*misc*",
            "*custom_dynalib_keywords*",
            "*system_models.py",
            "*define_function_strings.py",
            "*heart_decks.py",
            "*keyword_module.py",
            "*material_keywords.py",
        ],
    },
}

# Sphinx extensions
extensions = [
    "numpydoc",
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",
    "sphinx_autodoc_typehints",
    "sphinx_gallery.gen_gallery",
    "sphinxcontrib.video",
    "sphinx_design",
    "sphinx_jinja",
    "sphinx.ext.autodoc",
    "ansys_sphinx_theme.extension.autoapi",
]

# static path
html_static_path = ["_static"]

# custom css file
html_css_files = ["custom.css"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# Configuration for InterSphinx
# -----------------------------------------------------------------------------
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

# Configuration for Numpydoc
# -----------------------------------------------------------------------------
numpydoc_validate = True
numpydoc_xref_param_type = True
numpydoc_show_class_members = False
# https://numpydoc.readthedocs.io/en/latest/validation.html
numpydoc_validation_checks = {
    "GL06",  # Found unknown section
    "GL07",  # Sections are in the wrong order.
    # "GL08",  # The object does not have a docstring
    "GL09",  # Deprecation warning should precede extended summary
    "GL10",  # reST directives {directives} must be followed by two colons
    "SS01",  # No summary found
    "SS02",  # Summary does not start with a capital letter
    # "SS03",  # Summary does not end with a period
    "SS04",  # Summary contains heading whitespaces
    # "SS05",  # Summary must start with infinitive verb, not third person
    "RT02",  # The first line of the Returns section should contain only the
    # type, unless multiple values are being returned"
}

# Configuration for Sphinx gallery
# -----------------------------------------------------------------------------
sphinx_gallery_conf = {
    "examples_dirs": "../../examples",
    "gallery_dirs": "examples",
}

# Configuration for Sphinx autoapi
# -----------------------------------------------------------------------------

suppress_warnings = [
    "autoapi.python_import_resolution",
    "autosectionlabel.*",
    "config.cache",
    "design.fa-build",
]

# Configuration for Sphinx Autodoc Typehints
# -----------------------------------------------------------------------------
typehints_defaults = "comma"
simplify_optional_unions = False

# Common content for every RST file such as links
rst_epilog = ""
links_filepath = Path(__file__).parent.absolute() / "links.rst"
rst_epilog += links_filepath.read_text(encoding="utf-8")

exclude_patterns = ["links.rst"]

# Configuration for Jinja
# -----------------------------------------------------------------------------
jinja_globals = {
    "PYANSYS_HEART_VERSION": version,
}

# Get list of tox environments and add to jinja context
envs = subprocess.run(["tox", "list", "-q"], capture_output=True, text=True).stdout.splitlines()
envs.remove("default environments:")
envs.remove("additional environments:")
envs.remove("")

jinja_contexts = {
    "toxenvs": {
        "envs": envs,
    }
}
