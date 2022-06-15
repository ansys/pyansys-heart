Examples
========
This page contains examples of pyheart-lib usage

Extracting a simulation mesh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Use the preprocessor to extract the simulation files of a Bi-Ventricular heart model which includes 
mechanics and a system model:

.. code:: python

    import os
    from ansys.heart.preprocessor.heart_model import HeartModel
    from ansys.heart.preprocessor.model_information import ModelInformation
    from ansys.heart.writer.dynawriter import MechanicsDynaWriter

    case_path = os.path.join(
        "pyheart-lib",
        "tests",
        "heart", 
        "assets", 
        "cases", 
        "01", 
        "01.case" )
    
    work_dir = './bi_ventricle_model'


    # Valid models include: 
    #   [ "LeftVentricle", 
    #     "BiVentricle", 
    #     "FourChamber" ]
    # Define the necessary model information:
    model_info = ModelInformation(
        model_type = "BiVentricle",
        database_name = "Strocchi2020",
        path_original_mesh = case_path,
        working_directory = work_dir,
    )    
    model_info.mesh_size = 2.0
    model_info_path = os.path.join( work_dir, "model_info.json" )    

    # initialize HeartModel object and extract simulation mesh
    model = HeartModel( model_info )
    model.extract_simulation_mesh()

    # store model information
    model.dump_model( 
        model_info_path, 
        clean_working_directory = True 
        )
    
    # Use writer class to write ls-dyna .k files
    export_dir = os.path.join( work_dir, "ls-dyna-files" )
    writer = MechanicsDynaWriter ( model )
    writer.update()
    writer.export( export_dir )


    
    
    

