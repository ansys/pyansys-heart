"""Project installation script."""

"""pyheart-lib setup file."""
import codecs
from io import open as io_open
import os

from setuptools import find_namespace_packages, setup

HERE = os.path.abspath(os.path.dirname(__file__))
__version__ = None
version_file = os.path.join(HERE, "ansys", "heart", "_version.py")
with io_open(version_file, mode="r") as fd:
    exec(fd.read())


# Get the long description from the README file
# This is needed for the description on PyPI
def read(rel_path):
    """Get long description from the README file."""
    with codecs.open(os.path.join(HERE, rel_path), "r") as fp:
        return fp.read()


with open(os.path.join(HERE, "README.rst"), encoding="utf-8") as f:
    long_description = f.read()

packages = []
for package in find_namespace_packages(include="ansys*"):
    if package.startswith("ansys.heart"):
        packages.append(package)

setup(
    name="ansys-heart-lib",
    packages=packages,
    version=__version__,
    description="Python framework for heart modeling using ansys tools",
    long_description=long_description,
    # long_description_content_type='text/x-rst',
    url="https://github.com/pyansys/pyheart-lib",
    license="MIT",
    author="ANSYS, Inc.",
    maintainer="PyAnsys developers",
    maintainer_email="pyansys.support@ansys.com",
    # how to add dynalib?
    install_requires=[
        "gmsh==4.10.3",
        "h5py==3.6.0",
        "Jinja2==3.1.2",
        "matplotlib==3.5.2",
        "meshio==5.3.4",
        "numpy==1.21.6",
        "pandas==1.3.5",
        "scipy==1.7.3",
        "vtk==9.1.0",
        "tqdm==4.64.0",
    ],
    python_requires=">=3.7",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
"""Setup installation."""

# setup(
#     name="ansys-heart-lib",
#     version="0.1.dev0",
#     url="https://github.com/pyansys/pyheart-lib",
#     author="ANSYS, Inc.",
#     author_email="pyansys.support@ansys.com",
#     maintainer="PyAnsys developers",
#     maintainer_email="pyansys.maintainers@ansys.com",
#     classifiers=[
#         "Development Status :: 4 - Beta",
#         "Programming Language :: Python :: 3",
#         "License :: OSI Approved :: MIT License",
#         "Operating System :: OS Independent",
#     ],
#     license="MIT",
#     license_file="LICENSE",
#     description="Python framework for heart modeling using ansys tools",
#     long_description=open("README.rst").read(),
#     install_requires=["importlib-metadata >=4.0"],
#     python_requires=">=3.8",
#     packages=find_namespace_packages(where="src", include="ansys*"),
#     package_dir={"": "src"},
# )
