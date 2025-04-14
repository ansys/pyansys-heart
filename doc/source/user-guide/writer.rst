
.. _ref_writer:

Writer
======

The :attr:`DynaWriter <ansys.heart.writer.dynawriter>` base class generates LS-DYNA input files for different simulations.

Based on different applications, different writers must be created.

- :attr:`PurkinjeGenerationDynaWriter`: Generates an LS-DYNA input deck for creating a Purkinje network.
- :attr:`FiberGenerationDynaWriter`: Generates an LS-DYNA input deck for creating fiber orientation vectors.
- :attr:`MechanicsDynaWriter`: Generates an LS-DYNA input deck for mechanical simulations.
- :attr:`ZeroPressureMechanicsDynaWriter`: Generates an LS-DYNA input deck for stress-free configuration simulations.
- :attr:`ElectrophysiologyDynaWriter`: Generates an LS-DYNA input deck for electrophysiology simulations.
- :attr:`ElectroMechanicsDynaWriter`: Generates an LS-DYNA input deck for electrical-mechanical coupled simulations.

Here is a simple code example:

>>> # Get a heart model
>>> import ansys.health.heart.models as models
>>> model = models.HeartModel.load_model("path_to_model")

>>> import ansys.heart.writer.dynawriter as writers
>>> import copy
>>> # instantiate a Writer with the model and default settings
>>> settings = SimulationSettings().load_defaults()
>>> writer = writers.MechanicsDynaWriter(copy.deepcopy(model),settings=settings) # Writers may change the model, it's better to pass the copy of load_model
>>> writer.update()
>>> writer.export("path_to_folder")

