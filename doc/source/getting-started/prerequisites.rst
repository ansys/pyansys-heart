Prerequisites
=============

Operating system
----------------

- Windows 10
- Linux


Ansys tools
-----------

This framework was developed and tested under `Python 3.10 <Python310_>`_, `Python 3.11 <Python311_>`_, and `Python 3.12 <Python312_>`_.
Before starting the installation, run the ``python --version`` command and check that you are using one of the supported versions.

Software
--------

.. list-table:: **Required Ansys products**
  :widths: 200 300 200 400
  :header-rows: 1

  * - Product
    - Supported versions
    - Scope
    - Link to download

  * - Ansys Fluent
    - 2024 R1, 2024 R2, 2025 R1
    - Preprocessor
    - `Ansys Customer Portal`_

  * - Ansys DPF Server
    - 2024.1 (comes with the 2024 R1 installation), 2024.1rc1, 2024.2rc0
    - Postprocessor
    - `Ansys Customer Portal`_

  * - Ansys LS-DYNA
    - 16.0.0 IntelMPI double precision, 16.0.0 MSMPI double precision
    - Simulator
    - `Ansys LSDYNA Product Space_` or contact the `PyAnsys Core team <pyansys_core_>`_ to get more information.

.. note::

  Ansys Fluent is required for meshing. Also note that currently the postprocessor module is only compatible with Ansys DPF Servers 2024.1 (comes with the 2024 R1 installation), 2024.1rc1, and 2024.2rc0. Later versions are currently not supported. Hence, installing Ansys Fluent 2024 R1 is currently the most convenient.

