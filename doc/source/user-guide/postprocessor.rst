
.. _ref_postprocessor:


Postprocessor
=============

The post-processor :attr:`ansys.health.heart.post` module contains methods to post-process
the results of a cardiac simulation. The post-processor can be used to extract and visualize results
from the simulation output files, such as LS-DYNA's ``d3plot`` files and other relevant information.

The two main classes include:

- :attr:`ansys.health.heart.post.dpf_utils.D3plotReader`: a class for reading LS-DYNA ``d3plot`` files.
- :attr:`ansys.health.heart.post.dpf_utils.EPpostprocessor`: a class for post-processing electrophysiology simulation results.

These two classes add additional functionality on top of `PyDPF <https://dpf.docs.pyansys.com/>`_ to extract relevant
simulation results in a convenient way. For example, the :attr:`ansys.health.heart.post.dpf_utils.EPpostprocessor` class can be used to extract the transmembrane potential
and other relevant quantities from the simulation results.

>>> from ansys.health.heart.post.dpf_utils import D3plotReader, EPpostprocessor
>>> reader = D3plotReader("path-to-d3plot-file")
>>> # Get ep fields at step 10
>>> ep_fields = reader.get_ep_fields(at_step=10)
>>> # Get displacement at 200 ms
>>> displacement = reader.get_displacement(time=200)