.. _ref_simulator:

Simulator
=========

:attr:`Simulator <ansys.heart.simulator.simulator>` links different simulation steps for cardiac modeling. For example, in electrophysiology simulations, fiber orientation :attr:`BaseSimulator.compute_fibers` and Purkinje network :attr:`EPSimulator.compute_purkinje` are computed before launching the physical simulation. In mechanical analysis, you must compute the stress-free configuration :attr:`MechanicsSimulator.compute_stress_free_configuration` before running the simulation.

Different simulators must be created based on the application:

    - :attr:`BaseSimulator`: The parent class for all other simulators. It holds general methods, such as fiber generation.
    - :attr:`EPSimulator`: Runs electrophysiology cardiac simulations.
    - :attr:`MechanicsSimulator`: Runs mechanical cardiac simulations.
    - :attr:`EPMechanicsSimulator`: Runs electrical-mechanical coupled cardiac simulations.

Here is a simple code example:

>>> # get a heart model
>>> import ansys.health.heart.models as models
>>> model = models.HeartModel.load_model("path_to_model")

>>> # set up an LS-DYNA executable
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

Default modeling parameters are saved to the :attr:`<ansys.heart.simulator.settings.defaults>`_ attribute.
You can load them into the simulator:

.. code:: pycon

   >>> simulator.settings.load_defaults()
   >>> # we can print settings
   >>> print(simulator.settings.mechanics.analysis.end_time)
   800 millisecond
   >>> # Change it to 1600 ms
   >>> simulator.settings.mechanics.analysis.end_time = Quantity(1600, "ms")

Alternatively, you can load settings from a YAML file:

>>> simulator.settings.load("a-yaml-file")

