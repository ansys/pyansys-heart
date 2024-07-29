
.. _ref_components:

*********
simulator
*********

:attr:`Simulator <ansys.heart.simulator.simulator>` is used to link up different simulation steps for cardiac modeling. For example, for electrophysiological simulations, we typically need to firstly generate fiber orientation :attr:`\ <BaseSimulator.compute_fibers>` and Purkinje network :attr:`\ <EPSimulator.compute_purkinje>`, then start physical simulation. For mechanical simulations, it's necessary to compute stress free configuration :attr:`\ <MechanicsSimulator.compute_stress_free_configuration>` before running heart beat simulation.


Based on different applications, different simulators need to be created.

    - :attr:`BaseSimulator`, parent class for all other Simulators, and it holds general methods, like to generate fibers.
    - :attr:`EPSimulator`, used for running electrophysiological coupled cardiac simulation
    - :attr:`MechanicsSimulator`, used for running mecahnical coupled cardiac simulation
    - :attr:`EPMechanicsSimulator`, used for running electrical-mecahnical coupled cardiac simulation

A simple use example is given as the following:

>>> # Get a heart model
>>> import ansys.heart.preprocessor.models as models
>>> model = models.HeartModel.load_model("a-model-path")

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

Default modeling parameters are save :attr:`here <ansys.heart.simulator.settings.defaults>`, you can load them to the simulator:

.. code:: pycon

   >>> simulator.settings.load_defaults()
   >>> # we can print settings
   >>> print(simulator.settings.mechanics.analysis.end_time)
   800 millisecond
   >>> # let's change it to 1600 ms
   >>> simulator.settings.mechanics.analysis.end_time = Quantity(1600, "ms")

Alternatively, settings can be load from a yaml file as follow

>>> simulator.settings.load("a-yaml-file")

