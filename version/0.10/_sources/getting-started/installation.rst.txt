Installation
============

This guide helps you install PyAnsys Heart. It provides information on how to install
the package from PyPI, from a wheel file, or from the source code.

.. Note::

    If you do not have access to *PyAnsys Heart* you can follow the instructions under *Install from a wheel file*. You may need to request the wheel files from your Ansys contact.

.. Warning::

    Consider installing using a virtual environment to avoid conflicts with other packages.

PyPI
----

Before installing PyAnsys Heart ensure that you have the latest version
of the `pip <https://pip.pypa.io/en/stable/installation/>`_
package manager, run the following command:

.. code:: bash

    python -m pip install --upgrade pip

Then, to install PyAnsys Heart, run the following command:

.. code:: bash

    python -m pip install pyansys-heart

GitHub
------
To install the latest version of PyAnsys Heart from the source code,
clone the repository:

.. code:: bash

    git clone https://github.com/ansys/pyansys-heart
    cd pyansys-heart
    pip install -e .

to verify the installation, run the following command:

.. code:: bash

    python -m pip install tox
    tox

Install from a wheel file
-------------------------

if you lack the internet connection, you can install PyAnsys Heart from a wheel file.
You should install PyAnsys Heart by downloading the wheelhouse archive for your
corresponding machine architecture from the repository’s
`Release page <https://github.com/ansys/pyansys-heart/releases>`_.

Each release contains a wheel file for the corresponding Python version and
machine architecture. For example, to install the wheel file for
Python 3.10 on a Windows machine, run the following command:

.. code:: bash

    unzip pyansys-heart-v0.6.1-wheelhouse-windows-latest-3.12.zip wheelhouse
    pip install pyansys-heart -f wheelhouse --no-index --upgrade --ignore-installed

If you are on Windows with Python 3.12, unzip the wheelhouse archive to a wheelhouse
directory and then install using the same pip install command as in the preceding example.

Consider installing using a virtual environment to avoid conflicts with other packages. For more information,
refer to the `Python documentation <https://docs.python.org/3/library/venv.html>`_.

