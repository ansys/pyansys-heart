Prerequisites
=============

Operating system
----------------

- Windows 10
- Linux Ubuntu


Ansys tools
-----------

This framework was developed and tested under `Python310`_, `Python311`_ and `Python312`_ versions.
Before starting the
installation run ``python --version`` and check that it fits with the supported versions.

Software
--------

.. list-table:: Required Ansys products
  :widths: 200 300 200 400
  :header-rows: 1

  * - Product
    - Supported versions
    - Scope
    - Link to download

  * - Ansys Fluent
    - R24R1, R24R2, R25R1
    - Pre-processor
    - `Ansys Customer Portal`_

  * - Ansys DPF Server
    - 2024.1 (R24R1 install), 2024.1rc1, 2024.2rc0
    - Post-processor
    - `Ansys Customer Portal`_

  * - Ansys LS-DYNA
    - R16.0
    - Simulator
    - Contact `PyAnsys Core team <pyansys_core_>`_ to get more information

.. Note::

    Fluent is required for meshing. Also note that currently the postprocessor module is only compatible with Ansys DPF Servers 2024.1 (comes with R24R1 installation), 2024.1rc1 and 2024.2rc0. Later versions are currently not supported. Hence installing Ansys Fluent R24R1 is currently the most convenient.

