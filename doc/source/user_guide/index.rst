==========
User guide
==========
This section can contain a basic user guide for `PyAnsys-Heart <https://github.com/ansys/PyAnsys-Heart>`_

..
   This toctreemust be a top level index to get it to show up in
   pydata_sphinx_theme

.. toctree::
   :maxdepth: 1
   :hidden:

   prepocessor
   writer
   simulator
   postprocessor
   calibration



Overview
========

PyAnsys Heart is a high-level interface to LS-DYNA for heart modeling.
We introduce the relevant heart modeling features that are available and introduce the modular python library to set up and drive these simulations.


`Paper1`_

.. code:: python

    import os
    from ansys.mapdl.core import launch_mapdl

    path = os.getcwd()
    mapdl = launch_mapdl(run_location=path)

...~

Brief theory
============

Default option is explained. Literature review.
LS-DYNA manual is necessary to understand model


Misc
~~~~
- Ventrical fiber
- Atrial fiber
- UHC

Electrophisology
~~~~~~~~~~~~~~~~
- Karim 1
- Karim 2

Mechanics
~~~~~~~~~
More explanation in Writer ??

- Material:
   Cardiac tissue is modelled with *MAT_295. Users are suggested to read LS-DYNA manual for
more details. Generally, it can be split into 2 parts: passive and active.
For passive part, we use by default Holzapfel model for both isotropic part and anisotropic part.
For active part, different models are also available. If Mechanical model is employed, we use xx model (ACTYPE=1) and if
ep-mech coupled model is asked, the default model is XX (ACTYPE=3)

- Boundary conditions:
   valve
   pericardium

- Circulation model
   a lumped-parameter approach is available in LSDYNA, where blood flow in/out of the cavity – represented by an in-compressible volume – can be added as an additional constraint to update the pressure of the cavity

- Stress free configuration


Coupling
~~~~~~~~
Coupling is automatic, see EPMECA writer....
() effect is not considered


References
==========
.. _Paper1:
https://www.python.org/
