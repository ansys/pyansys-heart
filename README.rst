Pyheart lib
===========

Python framework for heart modeling using Ansys tools


How to install
--------------

At least two installation modes are provided: user and developer.

For users
^^^^^^^^^

User installation can be performed by running:

.. code:: bash

    python -m pip install ansys-heart-lib

For developers
^^^^^^^^^^^^^^

Installing Pyheart lib in developer mode allows
you to modify the source and enhance it.

Before contributing to the project, please refer to the `PyAnsys Developer's guide`_. You will 
need to follow these steps:

1. Start by cloning this repository:

    .. code:: bash

        git clone https://github.com/pyansys/pyheart-lib

2. Create a fresh-clean Python environment and activate it. Refer to the
   official `venv`_ documentation if you require further information:

    .. code:: bash

        # Create a virtual environment
        python -m venv .venv

        # Activate it in a POSIX system
        source .venv/bin/activate

        # Activate it in Windows CMD environment
        .venv\Scripts\activate.bat

        # Activate it in Windows Powershell
        .venv\Scripts\Activate.ps1

3. Make sure you have the latest version of `pip`_:

    .. code:: bash

        python -m pip install -U pip

4. Install the project in editable mode:

    .. code:: bash
    
        python -m pip install --editable ansys-heart-lib

5. Install additional requirements (if needed):

     .. code:: bash

        python -m pip install -r requirements/requirements_build_pip.txt
        python -m pip install -r requirements/requirements_doc.txt
        python -m pip install -r requirements/requirements_tests.txt


   Download and install `dynalib`_ in development mode in your environment:
     .. code:: bash

        python -m pip install -e /path/to/dynalib-root

6. Finally, verify your development installation by running:

    .. code:: bash
        
        python -m pip install -r requirements/requirements_tests.txt
        pytest tests -v


Style and Testing
-----------------

If required, you can always call the style commands (`black`_, `isort`_,
`flake8`_...) or unit testing ones (`pytest`_) from the command line. However,
this does not guarantee that your project is being tested in an isolated
environment, which is another reason to consider using `tox`_.


Documentation
-------------

For building documentation, you can either run the usual rules provided in the
`Sphinx`_ Makefile, such us:

.. code:: bash

    python -m pip install -r requirements/requirements_doc.txt
    make -C doc/ html

    # subsequently open the documentation with (under Linux):
    your_browser_name doc/html/index.html

Distributing
------------

If you would like to create either source or wheel files, start by installing
the building requirements:

.. code:: bash

    python -m pip install -r requirements/requirements_build.txt

Then, you can execute:

    .. code:: bash

        python -m build
        python -m twine check dist/*


.. LINKS AND REFERENCES
.. _black: https://github.com/psf/black
.. _flake8: https://flake8.pycqa.org/en/latest/
.. _isort: https://github.com/PyCQA/isort
.. _PyAnsys Developer's guide: https://dev.docs.pyansys.com/
.. _pre-commit: https://pre-commit.com/
.. _pytest: https://docs.pytest.org/en/stable/
.. _Sphinx: https://www.sphinx-doc.org/en/master/
.. _pip: https://pypi.org/project/pip/
.. _tox: https://tox.wiki/
.. _venv: https://docs.python.org/3/library/venv.html
.. _dynalib: https://github.com/pyansys/dynalib
