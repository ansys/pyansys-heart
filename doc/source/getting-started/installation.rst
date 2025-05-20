Installation
============

This page explains how to install PyAnsys Heart from PyPI, the GitHub source code, or a wheel file.

.. note::

    If you do not have access to PyAnsys Heart, follow the instructions in :ref:`install_from_wheel`.
    You might need to request the wheel file from your Ansys contact.

.. warning::

    Consider installing using a virtual environment to avoid conflicts with other packages. For more information,
    see `Creation of virtual environments <Python documentation_>`_ in the Python documentation.

Install from PyPI
-----------------

Before you install PyAnsys Heart, ensure that you have the latest version
of the `pip`_ package manager:

.. code:: bash

    python -m pip install --upgrade pip

Then, install PyAnsys Heart:

.. code:: bash

    python -m pip install ansys-health-heart

Install from GitHub source code
-------------------------------

To install the latest version of PyAnsys Heart from the GitHub source code,
clone the repository:

.. code:: bash

    git clone https://github.com/ansys/pyansys-heart
    cd pyansys-heart
    pip install -e .

Then, verify the installation:

.. code:: bash

    python -m pip install tox
    tox

.. _install_from_wheel:

Install from a wheel file
-------------------------

If you do not have an internet connection, you can install PyAnsys Heart from a wheel file.
Download the wheelhouse archive for your corresponding machine architecture
from the `repository’s Releases page <PyAnsys Heart release page_>`_.

Each release contains a wheel file for the corresponding Python version and
machine architecture. For example, to install the wheel file for
Python 3.10 on a Windows machine, run the following commands:

.. code:: bash

    unzip pyansys-heart-v0.6.1-wheelhouse-windows-latest-3.12.zip wheelhouse
    pip install pyansys-heart -f wheelhouse --no-index --upgrade --ignore-installed

If you are on Windows with Python 3.12, unzip the wheelhouse archive to a wheelhouse
directory and then install using the same ``pip install`` command as in the preceding example.
