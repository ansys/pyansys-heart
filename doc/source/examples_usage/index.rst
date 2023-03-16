Examples 
========

This page contains examples of pyheart-lib usage.

Downloading a case
^^^^^^^^^^^^^^^^^^
Download and extract a full heart mesh from the public database of 24 pathological hearts
by Strocchi et al (2020).

.. literalinclude:: ../../../examples/heart/preprocessor/example_downloader.py
.. figure:: ../images/mesh.jpg

Pre-process data for heart simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example deals with the preprocessing part of the heart simulation

Import all necessary modules from ansys.heart
"""""""""""""""""""""""""""""""""""""""""""""
.. code:: 

    import os
    import pathlib

    import ansys.heart.preprocessor.models as models
    from ansys.heart.simulator.support import run_preprocessor
    import ansys.heart.writer.dynawriter as writers

Initialization of working directories
"""""""""""""""""""""""""""""""""""""
.. code::

    path_to_case = os.path.join(
        pathlib.Path(__file__).parents[3], "downloads\\Strocchi2020\\01\\01.case"
    )
    workdir = os.path.join(pathlib.Path(path_to_case).parent, "BiVentricle")

    path_to_model = os.path.join(workdir, "heart_model.pickle")


Here we have a case from Strocchi2020 database. 
The preprocessor creates a path to a new folder called "BiVentricle" where the preprocessing results will be downloaded.

Preprocessing
"""""""""""""
.. code::

    model = run_preprocessor(
        model_type=models.BiVentricle,
        database="Strocchi2020",
        path_original_mesh=path_to_case,
        work_directory=workdir,
        path_to_model=path_to_model,
        mesh_size=2.0,
    )

run_preprocessor is filled with parameters like the database corresponding to the case studied, the mesh size and the working directories.

Loading the model
"""""""""""""""""
.. code::

    model = models.HeartModel.load_model(path_to_model)
    if not isinstance(model, models.HeartModel):
        exit()
    model.info.workdir = workdir

Can be performed when a model has already been created by the preprocessor. 

Loop for creating 4 folders containing LS-Dyna files
""""""""""""""""""""""""""""""""""""""""""""""""""""
This part uses the dynawriter module to create folders with LS-Dyna files associated.
The 4 folders created are about:

1. Electrophysiology

2. Mechanics

3. Zero pressure configuration computation

4. Fiber generation

5. Purkinje generation

.. code::

    for writer in (
        writers.ElectrophysiologyDynaWriter(model),
        writers.MechanicsDynaWriter(model),
        writers.ZeroPressureMechanicsDynaWriter(model),
        writers.FiberGenerationDynaWriter(model),
        writers.PurkinjeGenerationDynaWriter(model),
    ):
        exportdir = os.path.join(
            writer.model.info.workdir,
            writer.__class__.__name__.lower().replace("dynawriter", ""),
        )

        writer.model.mesh.write_to_vtk(
            os.path.join(writer.model.info.workdir, "volume_mesh.vtk")
        )
        writer.update()
        writer.export(exportdir)
    print("done")

Electrophysiology simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example of an electrophysiology simulation comprising: 

1. Fiber orientation computation and plot

2. Purkinje construction and plot

3. Electrophysiology simulation

.. literalinclude:: ../../../examples/heart/simulator/example_simulator_EP.py


Mechanics simulation
^^^^^^^^^^^^^^^^^^^^

Example of a mechanics simulation comprising: 

1. Fiber orientation computation and plot

2. Stress free configuration

3. Mechanics simulation

.. literalinclude:: ../../../examples/heart/simulator/example_simulator_Mechanics.py