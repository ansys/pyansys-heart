PyAnsys Heart
=============
|pyansys| |python310| |python311| |GH-CI| |MIT| |black| |pre-commit|

.. |pyansys| image:: https://img.shields.io/badge/Py-Ansys-ffc107.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAIAAACQkWg2AAABDklEQVQ4jWNgoDfg5mD8vE7q/3bpVyskbW0sMRUwofHD7Dh5OBkZGBgW7/3W2tZpa2tLQEOyOzeEsfumlK2tbVpaGj4N6jIs1lpsDAwMJ278sveMY2BgCA0NFRISwqkhyQ1q/Nyd3zg4OBgYGNjZ2ePi4rB5loGBhZnhxTLJ/9ulv26Q4uVk1NXV/f///////69du4Zdg78lx//t0v+3S88rFISInD59GqIH2esIJ8G9O2/XVwhjzpw5EAam1xkkBJn/bJX+v1365hxxuCAfH9+3b9/+////48cPuNehNsS7cDEzMTAwMMzb+Q2u4dOnT2vWrMHu9ZtzxP9vl/69RVpCkBlZ3N7enoDXBwEAAA+YYitOilMVAAAAAElFTkSuQmCC
   :target: https://docs.pyansys.com/
   :alt: PyAnsys

.. |python310| image:: https://img.shields.io/badge/Python-3.10-blue
   :target: https://www.python.org/downloads/release/python-3100/
   :alt: Python310

.. |python311| image:: https://img.shields.io/badge/Python-3.11-blue
   :target: https://www.python.org/downloads/release/python-3110/
   :alt: Python311

.. |python312| image:: https://img.shields.io/badge/Python-3.12-blue
   :target: https://www.python.org/downloads/release/python-3120/
   :alt: Python312

.. |GH-CI| image:: https://github.com/ansys/pyansys-heart/actions/workflows/ci_cd.yml/badge.svg
   :target: https://github.com/ansys/pyansys-heart/actions/workflows/ci_cd.yml
   :alt: GH-CI

.. |MIT| image:: https://img.shields.io/badge/license-MIT-yellow
   :target: https://opensource.org/blog/license/mit
   :alt: MIT

.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg?style=flat
   :target: https://github.com/psf/black
   :alt: Black

.. |pre-commit| image:: https://results.pre-commit.ci/badge/github/ansys/pyansys-heart/main.svg
   :target: https://results.pre-commit.ci/latest/github/ansys/pyansys-heart/main
   :alt: pre-commit.ci

PyAnsys Heart is a `Python`_ framework for heart modeling using Ansys tools.

.. Note::

    Please read LICENSE file before using this package.


Prerequisites
--------------

Operating system
^^^^^^^^^^^^^^^^

- Windows 10
- Linux Ubuntu


Ansys tools
^^^^^^^^^^^

This framework was developed and tested under |Python310|, |Python311|, and |Python312|. Before starting the
installation run ``python --version`` and check that it fits with the supported versions.

Software
^^^^^^^^

.. list-table:: Required Ansys products
  :widths: 200 300 200 400
  :header-rows: 1

  * - Product
    - Versions
    - Scope
    - Link to download

  * - Ansys Fluent
    - R24 R1
    - Pre-processor
    - `Ansys Customer Portal`_

  * - Ansys LS-DYNA
    - R16.0
    - Simulator
    - Contact us for appropriate version

.. Note::

    Fluent is required for meshing.

How to install
--------------

In user mode
^^^^^^^^^^^^

Request the value of the ``PYANSYS_PYPI_PRIVATE_READ_PAT`` token by sending an
email to `pyansys.core@ansys.com <mailto:pyansys.core@ansys.com>`_,
then fill in the its value in place of ``<TOKEN>`` to export the following environment variables:

1. Create a fresh Python environment and activate it. Make sure you use one of
    the supported Python versions. Refer to the official `venv`_  or `conda`_ documentation
    if you require further information:

Using `venv`_:

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

Using `conda`_:

.. code:: bash

      # Create virtual environment with a given Python version
      conda create --name my-venv python=3.10
      # Activate the environment
      conda activate my-venv

2. Install the project in your virtual environment, specifying the right token:

.. code:: bash

    pip install ansys-dpf-core==0.10.0
    pip install pyansys-heart==0.4.0 --index-url=https://<TOKEN>@pkgs.dev.azure.com/pyansys/_packaging/pyansys/pypi/simple/

.. note::

    You can also use environment variables to format the `-index-url` with the token value and URL.


In editable mode
^^^^^^^^^^^^^^^^

Installing PyAnsys-Heart in developer mode allows
you to modify the source and enhance it.

Before contributing to the project, please refer to the `PyAnsys Developer's guide`_. You will
need to follow these steps:

1. Start by cloning this repository:

.. code:: bash

    git clone https://github.com/ansys/pyansys-heart

Since this is a private repository you may need to provide your github username.
Alternatively you can download and unpack the zip file from `PyAnsys Heart`_

2. Create a fresh Python environment and activate it. Make sure you use one of
    the supported Python versions. Refer to the official `venv`_  or `conda`_ documentation
    if you require further information:

Using `venv`_:

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

Using `conda`_:

.. code:: bash

    # Create virtual environment with a given Python version
    conda create --name my-venv python=3.10
    # Activate the environment
    conda activate my-venv

3. Make sure you have the latest version of `pip`_ installed in your virtual environment.

.. code:: bash

    python -m pip install -U pip

4. Install dynalib 0.1.0 into your virtual environment with the following command. Request the appropriate private pypi token from pyansys.core@ansys.com.

.. code:: bash

    # latest version
    pip install dynalib==0.1.0 --index-url=https://<TOKEN>@pkgs.dev.azure.com/pyansys/_packaging/pyansys/pypi/simple/

Install the project in editable mode by pointing to the right location:

.. code:: bash

    python -m pip install --editable .

5. Install additional requirements (if needed):

.. code:: bash

    # dependencies for local doc building
    python -m pip install -e .[doc]
    # dependencies needed for (unit) testing
    python -m pip install -e .[tests]

6. You may verify your development version by running all or a set of tests:

.. code:: bash

    python -m pip install -e .[tests]

    # run quick tests
    python -m pytest -v -m "not requires_fluent or (not extract_models)"

    # run tests requiring Fluent
    python -m pytest -v -m requires_fluent

    # run all tests
    pytest tests -v


Style and Testing
-----------------

If required, you can always call the style commands (`black`_, `isort`_,
`flake8`_...) or unit testing ones (`pytest`_) from the command line. Alternatively, you can
use `pre-commit`_, which will ensure that all style requirements are met. However,
this does not guarantee that your project is being tested in an isolated
environment, which is another reason to consider using `tox`_.


Documentation and issues
------------------------
Documentation for the latest stable release of PyAnsys Heart is hosted at `documentation`_.

In the upper right corner of the documentation's title bar, there is an option for switching from
viewing the documentation for the latest stable release to viewing the documentation for the
development version or previously released versions.

On the `PyAnsys Heart Issues <https://github.com/ansys/pyansys-heart/issues>`_ page,
you can create issues to report bugs and request new features. On the `PyAnsys Heart Discussions
<https://github.com/ansys/pyansys-heart/discussions>`_ page or the `Discussions <https://discuss.ansys.com/>`_
page on the Ansys Developer portal, you can post questions, share ideas, and get community feedback.

To reach the project support team, email `pyansys.core@ansys.com <mailto:pyansys.core@ansys.com>`_.


License
-------

PyAnsys Heart is licensed under the MIT license. Please refer to the `LICENSE` file for more information.
PyAnsys Heart makes no commercial claim over any Ansys products whatsoever.
This library extends the functionality of the listed Ansys products by adding a Python interface
without changing the core behavior or licensing of the original products. This library requires
legally licensed copies of the involved Ansys products.


.. LINKS AND REFERENCES
.. _Python: https://www.python.org/
.. _PyAnsys Heart: https://github.com/ansys/pyansys-heart
.. _Ansys Customer Portal: https://support.ansys.com/Home/HomePage
.. _dpf-server: https://download.ansys.com/Others/DPF%20Pre-Release
.. _black: https://github.com/psf/black
.. _flake8: https://flake8.pycqa.org/en/latest/
.. _isort: https://github.com/PyCQA/isort
.. _pre-commit: https://pre-commit.com/
.. _PyAnsys Developer's guide: https://dev.docs.pyansys.com/
.. _pre-commit: https://pre-commit.com/
.. _pytest: https://docs.pytest.org/en/stable/
.. _Sphinx: https://www.sphinx-doc.org/en/master/
.. _pip: https://pypi.org/project/pip/
.. _tox: https://tox.wiki/
.. _venv: https://docs.python.org/3/library/venv.html
.. _conda: https://docs.conda.io/en/latest/
.. _documentation: https://heart.docs.pyansys.com/
