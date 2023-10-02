|pyansys| |python| |MIT|

.. |pyansys| image:: https://img.shields.io/badge/Py-Ansys-ffc107.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAIAAACQkWg2AAABDklEQVQ4jWNgoDfg5mD8vE7q/3bpVyskbW0sMRUwofHD7Dh5OBkZGBgW7/3W2tZpa2tLQEOyOzeEsfumlK2tbVpaGj4N6jIs1lpsDAwMJ278sveMY2BgCA0NFRISwqkhyQ1q/Nyd3zg4OBgYGNjZ2ePi4rB5loGBhZnhxTLJ/9ulv26Q4uVk1NXV/f///////69du4Zdg78lx//t0v+3S88rFISInD59GqIH2esIJ8G9O2/XVwhjzpw5EAam1xkkBJn/bJX+v1365hxxuCAfH9+3b9/+////48cPuNehNsS7cDEzMTAwMMzb+Q2u4dOnT2vWrMHu9ZtzxP9vl/69RVpCkBlZ3N7enoDXBwEAAA+YYitOilMVAAAAAElFTkSuQmCC
   :target: https://docs.pyansys.com/
   :alt: PyAnsys

.. |python| image:: https://img.shields.io/badge/Python-3.8-blue
   :target: https://www.python.org/downloads/release/python-380/
   :alt: Python

.. |MIT| image:: https://img.shields.io/badge/license-MIT-yellow
   :target: https://opensource.org/license/mit/
   :alt: MIT

PyAnsys-Heart is a `Python`_ framework for heart modeling using Ansys tools.

.. Note::

    Please read LICENSE file before using this package.


Prerequisites
=============

Operating system
----------------

- Windows 10
- Linux Ubuntu


Ansys tools
-----------

This framework was developed and tested under Python 3.8 version. Before starting the installation run ``python --version`` and check that it fits with the supported versions.

Software
--------

  .. list-table:: Required Ansys products
    :widths: 200 300 200 400
    :header-rows: 1

    * - Product
      - Versions
      - Scope
      - Link to download

    * - Ansys Fluent
      - R22 R2
      - Pre-processor
      - `Ansys Customer Portal`_

    * - Ansys DPF Server
      - R24 R1-pre0
      - Post-processor
      - `DPF-Server`_

    * - Ansys LSDYNA
      - DEV-97584 or greater
      - Simulator
      - Contact us for latest working version

.. note::
    Fluent is required for meshing and the Ansys DPF Server for post-processing electrophysiology
    and mechanical results. Currently we advice to use the pre-release version of `DPF-Server`_ since support
    for `d3plot` result files is updated frequently.

How to install
==============

In user mode
------------

.. warning::

    Installing as user-only is not yet supported.

.. User installation can be performed by running:

.. .. code:: bash

..     python -m pip install ansys-heart-lib

In editable mode
----------------

Installing PyAnsys-Heart in developer mode allows
you to modify the source and enhance it.

Before contributing to the project, please refer to the `PyAnsys Developer's guide`_. You will
need to follow these steps:

1. Start by cloning this repository:

    .. code:: bash

        git clone https://github.com/ansys/pyansys-heart

   Since this is a private repository you may need to provide your github username.
   Alternatively you can download and unpack the zip file from `PyAnsys-Heart`_

2. Create a fresh-clean Python environment and activate it. Make sure you use one of the supported Python versions. Refer to the
   official `venv`_  or `conda`_ documentation if you require further information:

   Through `venv`_:

    .. code:: bash

        # Create a virtual environment
        python -m venv .venv

        # Activate environment:

        # POSIX systems:
        source .venv/bin/activate

        # Windows cmd shell:
        .venv\Scripts\activate.bat

        # or in Windows powershell
        .venv\Scripts\Activate.ps1

   Through the virtual environment manager `conda`_

    .. code:: bash

        # Create virtual environment with a given Python version
        conda create --name my-venv python=3.8

        # Activate the environment
        conda activate my-venv

3. Make sure you have the latest version of `pip`_ installed in your virtual environment.

    .. code:: bash

        python -m pip install -U pip

4. Install the project in editable mode by pointing to the right location:

    .. code:: bash

        python -m pip install --editable .

   Install a version of dynalib into your virtual environment.

    .. code:: bash

        # latest version
        pip install git+https://github.com/ansys/dynalib.git@main

   or if encountering issues with dynalib you can install a specific version

        pip install git+https://github.com/ansys/dynalib.git@afce06ba178888d992ff51838ca521abb824c8ab


5. Install additional requirements (if needed):

     .. code:: bash

        # dependencies for local doc building
        python -m pip install .[doc]
        # dependencies needed for (unit) testing
        python -m pip install .[tests]

6. You may verify your development version by running all or a set of unit-tests:

    .. code:: bash

        python -m pip install .[tests]

        # run quick tests
        python -m pytest -v -m "not requires_fluent and not local"

        # run tests requiring Fluent
        python -m pytest -v -m requires_fluent

        # run all tests
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
detailed description of the library or for specific examples.

For building documentation, you can either run the usual rules provided in the
`Sphinx`_ Makefile, such us:

.. code:: bash

    # install any dependencies for building the documentation.
    python -m pip install .[doc]

    # Linux
    make -C doc/ html

    # Windows
    cd doc/
    make.bat html

subsequently open the documentation by opening `doc/html/index.html`:


Distributing
============

If you would like to create either source or wheel files, you can execute:

.. code:: bash

    python -m pip install .


Licensing terms
===============

PyAnsys-Heart is licensed under the MIT license:

    MIT License

    Copyright (c) 2023 ANSYS, Inc. All rights reserved.

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

PyAnsys-Heart makes no commercial claim over any Ansys products whatsoever. This library extends the functionality of the listed Ansys products by adding a Python interface without changing the core behavior or licensing of the original products. This library requires legally licensed copies of the involved Ansys products.


.. LINKS AND REFERENCES
.. _Python: https://www.python.org/
.. _PyAnsys-Heart: https://github.com/ansys/pyansys-heart
.. _Ansys Customer Portal: https://support.ansys.com/Home/HomePage
.. _dpf-server: https://download.ansys.com/Others/DPF%20Pre-Release
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
.. _dynalib: https://github.com/ansys/dynalib
.. _conda: https://docs.conda.io/en/latest/
.. _documentation: https://heart.docs.pyansys.com/