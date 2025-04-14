
.. _ref_simulator:

*********
Simulator
*********

:attr:`Simulator <ansys.heart.simulator.simulator>` is used to link up different simulation steps for cardiac modeling. For example, for electrophysiology simulations, fiber orientation :attr:`BaseSimulator.compute_fibers` and Purkinje network :attr:`EPSimulator.compute_purkinje` are computed before launching the physical simulation. In mechanical analysis, it is necessary to compute the stress free configuration :attr:`MechanicsSimulator.compute_stress_free_configuration` before running the simulation.


Based on different applications, different simulators need to be created.

    - :attr:`BaseSimulator`, parent class for all other Simulators, it holds general methods, like fiber generation.
    - :attr:`EPSimulator`, used for running electrophysiology cardiac simulation
    - :attr:`MechanicsSimulator`, used for running mechanical cardiac simulation
    - :attr:`EPMechanicsSimulator`, used for running electrical-mechanical coupled cardiac simulation

A simple usage example is given in the following:

>>> # Get a heart model
>>> import ansys.health.heart.models as models
>>> model = models.HeartModel.load_model("path_to_model")

>>> # Set up a LS-DYNA executable
>>> from ansys.heart.simulator.simulator import DynaSettings, MechanicsSimulator
>>> dyna_settings = DynaSettings(
    lsdyna_path=lsdyna_path,
    dynatype="intelmpi",
    num_cpus=8)

>>> # instantiate simulator
>>> simulator = EPSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory="output-path")

Default modeling parameters are saved :attr:`here <ansys.heart.simulator.settings.defaults>`, you can load them to the simulator:

.. code:: pycon

   >>> simulator.settings.load_defaults()
   >>> # we can print settings
   >>> print(simulator.settings.mechanics.analysis.end_time)
   800 millisecond
   >>> # let's change it to 1600 ms
   >>> simulator.settings.mechanics.analysis.end_time = Quantity(1600, "ms")

Alternatively, settings can be load from a yaml file as follow

>>> simulator.settings.load("a-yaml-file")

