.. _ref_simulator:

Simulator
=========

The :attr:`simulator <ansys.health.heart.simulator>` module links different simulation steps for cardiac modeling.
For example, in electrophysiology simulations, you compute fiber orientation and the Purkinje network using the :attr:`BaseSimulator.compute_fibers` and :attr:`EPSimulator.compute_purkinje` methods before running the physical simulation.
In mechanical simulations, you compute the stress-free configuration using the :attr:`MechanicsSimulator.compute_stress_free_configuration` method before running the simulation.

You create different simulators based on the application and physics of interest:

- :attr:`BaseSimulator <ansys.health.heart.simulator.BaseSimulator>`: This parent class provides general methods, such as fiber generation.
- :attr:`EPSimulator <ansys.health.heart.simulator.EPSimulator>`: This class derives from :attr:`BaseSimulator <ansys.health.heart.simulator.BaseSimulator>` and includes specific methods for cardiac electrophysiology simulations.
- :attr:`MechanicsSimulator <ansys.health.heart.simulator.MechanicsSimulator>`: This class derives from :attr:`BaseSimulator <ansys.health.heart.simulator.BaseSimulator>` and includes specific methods for mechanical cardiac simulations.
- :attr:`EPMechanicsSimulator <ansys.health.heart.simulator.EPMechanicsSimulator>`: This class derives from :attr:`BaseSimulator <ansys.health.heart.simulator.BaseSimulator>` and includes specific methods for electrical-mechanical coupled cardiac simulations.

Here is a simple code example:

Load a heart model.

.. code-block:: python

    import ansys.health.heart.models as models
    model = models.HeartModel.load_model("path_to_model.vtu", "path_to_info.partinfo.json")

Set up the LS-DYNA settings.

.. code-block:: python

    from ansys.health.heart.simulator.simulator import DynaSettings, MechanicsSimulator
    dyna_settings = DynaSettings(
        lsdyna_path="path-to-lsdyna-exe.exe",
        dynatype="intelmpi",
        platform="windows",
        num_cpus=8
    )

Instantiate the simulator.

.. code-block:: python

    simulator = EPSimulator(
        model=model,
        dyna_settings=dyna_settings,
        simulation_directory="output-path"
    )

The :attr:`settings <ansys.health.heart.settings.defaults>` module saves default modeling parameters. You can load these parameters into the simulator:

.. code-block:: python

    simulator.settings.load_defaults()
    # Print settings
    print(simulator.settings.mechanics.analysis.end_time)
    # Output: 800 millisecond
    # Change it to 1600 ms
    simulator.settings.mechanics.analysis.end_time = Quantity(1600, "ms")
    # Save to a YAML file
    simulator.settings.save("a-yaml-file.yml")

Alternatively, you can load settings from a YAML file:

.. code-block:: python

    simulator.settings.load("a-yaml-file.yml")

Finally, run the relevant steps before running the final simulation of the physics of interest:

.. code-block:: python

    simulator.compute_fibers()
    simulator.compute_purkinje()
    simulator.simulate()