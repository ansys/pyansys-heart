###########
Pyheart lib
###########
|pyansys| |python| |MIT|

.. |pyansys| image:: https://img.shields.io/badge/Py-Ansys-ffc107.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAIAAACQkWg2AAABDklEQVQ4jWNgoDfg5mD8vE7q/3bpVyskbW0sMRUwofHD7Dh5OBkZGBgW7/3W2tZpa2tLQEOyOzeEsfumlK2tbVpaGj4N6jIs1lpsDAwMJ278sveMY2BgCA0NFRISwqkhyQ1q/Nyd3zg4OBgYGNjZ2ePi4rB5loGBhZnhxTLJ/9ulv26Q4uVk1NXV/f///////69du4Zdg78lx//t0v+3S88rFISInD59GqIH2esIJ8G9O2/XVwhjzpw5EAam1xkkBJn/bJX+v1365hxxuCAfH9+3b9/+////48cPuNehNsS7cDEzMTAwMMzb+Q2u4dOnT2vWrMHu9ZtzxP9vl/69RVpCkBlZ3N7enoDXBwEAAA+YYitOilMVAAAAAElFTkSuQmCC
   :target: https://docs.pyansys.com/
   :alt: PyAnsys

.. |python| image:: https://img.shields.io/badge/Python-3.8-blue
   :target: https://pypi.org/project/pyheartlib/
   :alt: Python

.. |MIT| image:: https://img.shields.io/badge/license-MIT-yellow
   :target: https://opensource.org/license/mit/
   :alt: MIT

Pyheart lib is a python framework for heart modeling using Ansys tools.

Note: Please read LICENSE file before using this package.


Prerequisites
=============

Operation system
----------------

- Windows 10
- Linux Ubuntu

Python
------

Officially Python 3.8 version. Before starting the installation run ``python --version`` and check that it fits with the supported versions.

Software
--------

  .. list-table:: Required flagship products
    :widths: 300 200 300
    :header-rows: 1

    * - Product
      - Versions
      - Link to download

    * - Ansys Fluent
      - R22.2
      - `Ansys Customer Portal`_

    * - Ansys DPF Server
      - 2023 R2-pre1
      - `DPF`_

    * - Ansys LSDYNA
      - ???
      - ?

Notes: Fluent is required for meshing and Ansys DPF Server for post-process module and calibration module

How to install
==============

At least two installation modes are provided: user and developer.

For users (Not yet available!)
------------------------------

User installation can be performed by running:

.. code:: bash

    python -m pip install ansys-heart-lib

For developers
--------------

Installing Pyheart lib in developer mode allows
you to modify the source and enhance it.

Before contributing to the project, please refer to the `PyAnsys Developer's guide`_. You will 
need to follow these steps:

1. Start by cloning this repository:

    .. code:: bash

        git clone https://github.com/pyansys/pyheart-lib

2. Create a fresh-clean Python environment and activate it. Refer to the
   official `venv`_  or `conda`_ documentation if you require further information:

    Through `venv`_

    .. code:: bash

        # Create a virtual environment
        python -m venv .venv

        # Activate environment: 
        # - in POSIX system:
        source .venv/bin/activate
        # - in Windows cmd shell:
        .venv\Scripts\activate.bat
        # or in Windows powershell
        .venv\Scripts\Activate.ps1

    Through the virtual environment manager `conda`_

    .. code:: bash

        # Create virtual environment with specific python version
        conda create --name my-venv python=3.8

        # Activate environment
        conda activate my-venv        

3. Make sure you have the latest version of `pip`_

    .. code:: bash

        python -m pip install -U pip

4. Install the project in editable mode:

    .. code:: bash
    
        python -m pip install --editable pyheart-lib
    
    Install version of dynalib manually by

    .. code:: bash
        
        # latest version
        pip install git+https://github.com/pyansys/dynalib.git@main

        # or if encountering issues with dynalib you can install a specific working version
        pip install git+https://github.com/pyansys/dynalib.git@afce06ba178888d992ff51838ca521abb824c8ab

        # Otherwise you can install it in editable mode:
        git clone https://github.com/pyansys/dynalib.git
        cd dynalib
        pip install -e .


5. Install additional requirements (if needed):

     .. code:: bash

        python -m pip install -r requirements_build.txt
        python -m pip install -r requirements_docs.txt
        python -m pip install -r requirements_tests.txt

6. Finally, verify your development version after installation by running:

    .. code:: bash
        
        python -m pip install -r requirements_tests.txt
        pytest tests -v


Style and Testing
=================

If required, you can always call the style commands (`black`_, `isort`_,
`flake8`_...) or unit testing ones (`pytest`_) from the command line. However,
this does not guarantee that your project is being tested in an isolated
environment, which is another reason to consider using `tox`_.


Documentation
=============

Visit the `documentation`_ for a
detailed description of the library.
For building documentation, you can either run the usual rules provided in the
`Sphinx`_ Makefile, such us:

.. code:: bash

    python -m pip install -r requirements_docs.txt
    # Linux
    make -C doc/ html 
    # Windows
    cd doc/
    make html

    # subsequently open the documentation with (under Linux):
    your_browser_name doc/html/index.html

Distributing
============

If you would like to create either source or wheel files, start by installing
the building requirements:

.. code:: bash

    python -m pip install -r requirements_build.txt

Then, you can execute:

    .. code:: bash

        python -m build
        python -m twine check dist/*


.. LINKS AND REFERENCES
.. _Ansys Customer Portal: https://support.ansys.com/Home/HomePage
.. _dpf: https://download.ansys.com/Others/DPF%20Pre-Release
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
.. _conda: https://docs.conda.io/en/latest/
.. _documentation: https://heart.docs.pyansys.com/