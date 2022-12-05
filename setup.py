"""Project installation script."""
import subprocess
import sys

"""pyheart-lib setup file."""
import codecs
from io import open as io_open
import os

from setuptools import find_namespace_packages, setup
from setuptools.command.develop import develop


class PostDevelopCommand(develop):
    """Post-installation for development mode."""

    def run(self):
        """Post run to install dynalib."""
        develop.run(self)
        print("Installing dynalib...")
        subprocess.call("git clone https://github.com/pyansys/dynalib.git")
        subprocess.call("python -m pip install -e dynalib")

        if sys.version_info.minor == 7 or sys.version_info.minor == 8:
            print("Installing qd...")
            subprocess.call("python -m pip install qd")


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

with open("requirements_build.txt") as f:
    install_requires = f.read().splitlines()

# add these files as package data
# can test if files are indeed added to distribution by: python setup.py sdist bdist_wheel
package_data = {
    "ansys.heart.preprocessor.templates": [
        "fluent_meshing_template_improved_2.jou",
        "fluent_meshing_add_blood_mesh_template.jou",
    ],
    "ansys.heart.preprocessor": ["*.json"],
    "ansys.heart.writer.templates": ["*.json"],
    "ansys.heart.misc.paraview_macro": ["*.pvpy"],
    "ansys.heart.writer": ["calcium_from_EP.txt"],
    "ansys.heart.calibration": ["material.k", "PassiveCalibration.lsopt", "run.template"],
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
    # install dynalib
    cmdclass={"develop": PostDevelopCommand},
    python_requires=">=3.7",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
"""Setup installation."""
