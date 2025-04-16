
.. _ref_writer:

******
Writer
******

:attr:`DynaWriter <ansys.heart.writer.dynawriter>` is used to generate LS-DYNA input files for different simulations.

Based on different applications, different Writers need to be created.

    - :attr:`PurkinjeGenerationDynaWriter`, to generate a LS-DYNA input deck for creating Purkinje network.
    - :attr:`FiberGenerationDynaWriter`, to generate a LS-DYNA input deck for creating fibers orientation vectors.
    - :attr:`MechanicsDynaWriter`, to generate a LS-DYNA input deck for mechanical simulations
    - :attr:`ZeroPressureMechanicsDynaWriter`, to generate a LS-DYNA input deck for stress free configuration simulations
    - :attr:`ElectrophysiologyDynaWriter`, to generate a LS-DYNA input deck for electrophysiology simulations
    - :attr:`ElectroMechanicsDynaWriter`, to generate a LS-DYNA input deck for electrical-mecahnical coupled simulations

A simple use example is given as the following:

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

