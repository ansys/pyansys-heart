"""Project installation script."""

"""pyheart-lib setup file."""
import codecs
from io import open as io_open
import os

from setuptools import find_namespace_packages, setup

_THIS_FILE = os.path.abspath(os.path.dirname(__file__))
__version__ = None
version_file = os.path.join(_THIS_FILE, "ansys", "heart", "_version.py")
with io_open(version_file, mode="r") as fd:
    exec(fd.read())


# Get the long description from the README file
# This is needed for the description on PyPI
def read(rel_path):
    """Get long description from the README file."""
    with codecs.open(os.path.join(_THIS_FILE, rel_path), "r") as fp:
        return fp.read()


with open(os.path.join(_THIS_FILE, "README.rst"), encoding="utf-8") as f:
    long_description = f.read()

packages = []
for package in find_namespace_packages(include="ansys*"):
    if package.startswith("ansys.heart"):
        packages.append(package)

# list of required packages
install_requires = (
    [
        "h5py",
        "Jinja2>=3.1.2",
        "meshio==5.3.4",
        "numpy==1.21.6",
        # "pandas==1.3.5", pandas can be replaced by dynalib
        "scipy==1.7.3",
        "vtk==9.1.0",
        # "tqdm==4.64.0", # optional?
    ],
)

# add these files as package data
# can test if files are indeed added to distribution by: python setup.py sdist bdist_wheel
package_data = {
    "ansys.heart.preprocessor.templates": [
        "fluent_meshing_template_improved_2.jou",
        "fluent_meshing_add_blood_mesh_template.jou",
    ],
    "ansys.heart.preprocessor": ["*.json"],
}

setup(
    name="ansys-heart-lib",
    description="Python framework for heart modeling using ansys tools",
    packages=packages,
    # package_dir={"": "ansys"},
    include_package_data=True,
    package_data=package_data,
    version=__version__,
    long_description=long_description,
    install_requires=install_requires,
    tests_require=["pytest"],
    # long_description_content_type='text/x-rst',
    url="https://github.com/pyansys/pyheart-lib",
    license="MIT",
    author="ANSYS, Inc.",
    maintainer="PyAnsys developers",
    maintainer_email="pyansys.support@ansys.com",
    python_requires=">=3.7",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
"""Setup installation."""
