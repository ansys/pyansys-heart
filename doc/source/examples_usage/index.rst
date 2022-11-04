Examples 
---------

This page contains examples of pyheart-lib usage

Extracting a simulation mesh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Use the preprocessor to extract the simulation files of a Bi-Ventricular heart model which includes 
mechanics and a system model:

.. literalinclude:: ../../../examples/heart/preprocessor/example_doc_extract_biventricular_model.py

Pre-process data for heart simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example deals with the preprocessing part of the heart simulation

Import all necessary modules from ansys.heart
"""""""""""""""""""""""""""""""""""""""""""""
.. code:: 

    import os
    import pathlib
    import ansys.heart.preprocessor.models as models
    from ansys.heart.workflow.support import run_preprocessor
    import ansys.heart.writer.dynawriter as writers

Initialization of working directories
"""""""""""""""""""""""""""""""""""""
.. code::

    path_to_case = "D:\\Cristobal_cases\\01.vtk"
    workdir = os.path.join(pathlib.Path(path_to_case).parent, "BiVentricle")
    path_to_model = os.path.join(workdir, "heart_model.pickle")

Here we have a case from Cristobal2021 database. It has been downloaded as .vtk file. It also could be a .case file.
In path_to_case, indicate the path to the downloaded file.

Workdir creates a path to a new folder called "BiVentricle" where the preprocessing results will be downloaded.

Loop for preprocessing data
"""""""""""""""""""""""""""
.. code::

    model = run_preprocessor(
        model_type=models.BiVentricle,
        database="Cristobal2021",
        path_original_mesh=path_to_case,
        work_directory=workdir,
        path_to_model=path_to_model,
        mesh_size=2.0,
    )

run_preprocessor is filled with right parameters like the database corresponding to the case studied, the mesh size and the working directories.

Loop for creating 4 folders containing LS-Dyna files
""""""""""""""""""""""""""""""""""""""""""""""""""""
This part uses the dynawriter module to create folders with LS-Dyna files associated.
The 4 folders created are about :

1. the fiber generation

2. the Purkinje network generation

3. the zero-pressure mechanics

4. the complete mechanics


.. code::

    for writer in (
        writers.MechanicsDynaWriter(model, "ConstantPreloadWindkesselAfterload"),
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