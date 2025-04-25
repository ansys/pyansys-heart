.. _ref_postprocessor:

Postprocessor
=============

The :attr:`postprocessor <ansys.health.heart.post>` module provides methods to postprocess the results of a cardiac simulation. Use this module to extract and visualize results from simulation output files, such as LS-DYNA's ``d3plot`` files and other relevant information.

This module includes two main classes:

- :attr:`ansys.health.heart.post.dpf_utils.D3plotReader`: This class provides methods to read LS-DYNA ``d3plot`` files.
- :attr:`ansys.health.heart.post.dpf_utils.EPpostprocessor`: This class provides methods to postprocesses electrophysiology simulation results.

These classes build on `PyDPF <https://dpf.docs.pyansys.com/>`_ to make it easier to extract relevant simulation results. For example, you can use the :attr:`ansys.health.heart.post.dpf_utils.EPpostprocessor` class to extract the transmembrane potential and other quantities from simulation results.

Here is an example of how to use these classes:

>>> from ansys.health.heart.post.dpf_utils import D3plotReader, EPpostprocessor
>>> # Load a d3plot file
>>> reader = D3plotReader("path-to-d3plot-file")
>>> # Extract electrophysiology fields at step 10
>>> ep_fields = reader.get_ep_fields(at_step=10)
>>> # Extract displacement at 200 ms
>>> displacement = reader.get_displacement(time=200)