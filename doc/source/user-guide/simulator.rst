.. _ref_simulator:

Simulator
=========

The Simulator :attr:`<ansys.heart.simulator` module links different simulation steps for cardiac modeling. For example, in electrophysiology simulations, you compute fiber orientation and the Purkinje network using the :attr:`BaseSimulator.compute_fibers` and Purkinje network :attr:`EPSimulator.compute_purkinje` methods before you run the physical simulation. In mechanical simulations, you must compute the stress-free configuration using the :attr:`MechanicsSimulator.compute_stress_free_configuration` method before running the simulation.

Different simulators must be created based on the application:

- :attr:`BaseSimulator`: The parent class for all other simulators. It holds general methods, such as fiber generation.
- :attr:`EPSimulator`: Derives from :attr:`BaseSimulator` and has specific methods dedicated for cardiac electrophysiology simulations.
- :attr:`MechanicsSimulator`: Derives from :attr:`BaseSimulator` and has specific methods for mechanical cardiac simulations.
- :attr:`EPMechanicsSimulator`: Derives from :attr:`BaseSimulator` and has specific methods for electrical-mechanical coupled cardiac simulations.

Here is a simple code example:

>>> # load a heart model
>>> import ansys.health.heart.models as models
>>> model = models.HeartModel.load_model("path_to_model.vtu", "path_to_info.partinfo.json")

>>> # set up an LS-DYNA executable
>>> from ansys.heart.simulator.simulator import DynaSettings, MechanicsSimulator
>>> dyna_settings = DynaSettings(
    lsdyna_path="path-to-lsdyna-exe.exe",
    dynatype="intelmpi",
    platform="windows",
    num_cpus=8)

>>> # instantiate the simulator
>>> simulator = EPSimulator(
    model=model,
    dyna_settings=dyna_settings,
    simulation_directory="output-path")

Default modeling parameters are saved from the :attr:`ansys.heart.settings.defaults` module.
You can load them into the simulator:

.. code:: pycon

   >>> simulator.settings.load_defaults()
   >>> # Print settings
   >>> print(simulator.settings.mechanics.analysis.end_time)
   800 millisecond
   >>> # Change it to 1600 ms
   >>> simulator.settings.mechanics.analysis.end_time = Quantity(1600, "ms")
   >>> # Save to YAML file
   >>> simulator.settings.save("a-yaml-file.yml")

Alternatively, you can load settings from a YAML file:

>>> simulator.settings.load("a-yaml-file.yml")

Finally you can run relevant steps prior to running the final simulation of the physics of interest:

>>> simulator.compute_fibers()
>>> simulator.compute_purkinje()
>>> simulator.simulate()
