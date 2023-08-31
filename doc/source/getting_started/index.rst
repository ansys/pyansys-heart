Getting Started
===============

.. include:: ../../../README.rst

.. Getting Started
.. ===============

.. Python framework for heart modeling using Ansys tools


.. How to install
.. --------------

.. At least two installation modes are provided: user and developer.

.. For users (Not yet available!)
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. User installation can be performed by running:

.. .. code:: bash

..     python -m pip install ansys-heart-lib

.. For developers
.. ^^^^^^^^^^^^^^

.. Installing Pyheart lib in developer mode allows
.. you to modify the source and enhance it.

.. Before contributing to the project, please refer to the `PyAnsys Developer's guide`_. You will
.. need to follow these steps:

.. 1. Start by cloning this repository:

..     .. code:: bash

..         git clone https://github.com/pyansys/pyheart-lib

.. 2. Create a fresh-clean Python environment and activate it. Refer to the
..    official `venv`_  or `conda`_ documentation if you require further information:

..     Through `venv`_

..     .. code:: bash

..         # Create a virtual environment
..         python -m venv .venv

..         # Activate environment:
..         # - in POSIX system:
..         source .venv/bin/activate
..         # - in Windows cmd shell:
..         .venv\Scripts\activate.bat
..         # or in Windows powershell
..         .venv\Scripts\Activate.ps1

..     Through the virtual environment manager `conda`_

..     .. code:: bash

..         # Create virtual environment with specific python version
..         conda create --name my-venv python=3.8

..         # Activate environment
..         conda activate my-venv

.. 3. Make sure you have the latest version of `pip`_

..     .. code:: bash

..         python -m pip install -U pip

.. 4. Install the project in editable mode:

..     .. code:: bash

..         python -m pip install --editable pyheart-lib

..     Install version of dynalib manually by

..     .. code:: bash

..         # latest version
..         pip install git+https://github.com/pyansys/dynalib.git@main

..         # or if encountering issues with dynalib you can install a specific working version
..         pip install git+https://github.com/pyansys/dynalib.git@afce06ba178888d992ff51838ca521abb824c8ab

..         # Otherwise you can install it in editable mode:
..         git clone https://github.com/pyansys/dynalib.git
..         cd dynalib
..         pip install -e .

..     Alternatively, use

..     .. code:: bash

..         python setup.py develop

..     which will install dynalib (and qd) automatically.

..     Note 1: qd will be installed only if Python is 3.7 or 3.8.
..     Note 2: this option may fail in some cases, please resort back to pip install --editable and manual installation of `dynalib`_

.. 5. Install additional requirements (if needed):

..      .. code:: bash

..         python -m pip install -r requirements_build.txt
..         python -m pip install -r requirements_docs.txt
..         python -m pip install -r requirements_tests.txt

.. 6. Finally, verify your development version after installation by running:

..     .. code:: bash

..         python -m pip install -r requirements_tests.txt
..         pytest tests -v


.. Style and Testing
.. -----------------

.. If required, you can always call the style commands (`black`_, `isort`_,
.. `flake8`_...) or unit testing ones (`pytest`_) from the command line. However,
.. this does not guarantee that your project is being tested in an isolated
.. environment, which is another reason to consider using `tox`_.


.. Documentation
.. -------------

.. Visit the `documentation`_ for a
.. detailed description of the library.
.. For building documentation, you can either run the usual rules provided in the
.. `Sphinx`_ Makefile, such us:

.. .. code:: bash

..     python -m pip install -r requirements_docs.txt
..     # Linux
..     make -C doc/ html
..     # Windows
..     cd doc/
..     make html

..     # subsequently open the documentation with (under Linux):
..     your_browser_name doc/html/index.html

.. Distributing
.. ------------

.. If you would like to create either source or wheel files, start by installing
.. the building requirements:

.. .. code:: bash

..     python -m pip install -r requirements_build.txt

.. Then, you can execute:

..     .. code:: bash

..         python -m build
..         python -m twine check dist/*


.. .. LINKS AND REFERENCES
.. .. _black: https://github.com/psf/black
.. .. _flake8: https://flake8.pycqa.org/en/latest/
.. .. _isort: https://github.com/PyCQA/isort
.. .. _PyAnsys Developer's guide: https://dev.docs.pyansys.com/
.. .. _pre-commit: https://pre-commit.com/
.. .. _pytest: https://docs.pytest.org/en/stable/
.. .. _Sphinx: https://www.sphinx-doc.org/en/master/
.. .. _pip: https://pypi.org/project/pip/
.. .. _tox: https://tox.wiki/
.. .. _venv: https://docs.python.org/3/library/venv.html
.. .. _dynalib: https://github.com/pyansys/dynalib
.. .. _conda: https://docs.conda.io/en/latest/
.. .. _documentation: https://heart.docs.pyansys.com/