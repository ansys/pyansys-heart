
.. _ref_components:

******
writer
******

:attr:`DynaWriter <ansys.heart.writer.dynawriter>` is used to generate LS-DYNA input files for different simulations.

Based on different applications, different Writers need to be created.

    - :attr:`PurkinjeGenerationDynaWriter`, used to generate LS-DYNA input deck for creating Purkinje network.
    - :attr:`FiberGenerationDynaWriter`, used to generate LS-DYNA input deck for creating fibers orientatin vectors.
    - :attr:`MechanicsDynaWriter`, used to generate LS-DYNA input deck for mechanical simulations
    - :attr:`ZeroPressureMechanicsDynaWriter`, used to generate LS-DYNA input deck for stress_free_configuration simulations
    - :attr:`ElectrophysiologyDynaWriter`, used to generate LS-DYNA input deck for electrophysiological simulations
    - :attr:`ElectroMechanicsDynaWriter`, used to generate LS-DYNA input deck for electrical-mecahnical coupled simulations

A simple use example is given as the following:

>>> # Get a heart model
>>> import ansys.heart.preprocessor.models as models
>>> model = models.HeartModel.load_model("a-model-path")

>>> import ansys.heart.writer.dynawriter as writers
>>> import copy
>>> # instantiate a Writer with the model and default settings
>>> settings = SimulationSettings().load_defaults()
>>> writer = writers.MechanicsDynaWriter(copy.deepcopy(model),settings=settings) # Writers may change the model, it's better to pass the copy of load_model
>>> writer.update()
>>> writer.export("a-folder-path")

