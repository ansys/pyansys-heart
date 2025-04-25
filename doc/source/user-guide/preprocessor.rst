.. _ref_preprocessor:

Preprocessor
=============

This topic provides an overview of the :attr:`Preprocessor <ansys.health.heart.pre>` module. Use this module to preprocess a case from the following two public datasets:

- `A Publicly Available Virtual Cohort of Four-chamber Heart Meshes for Cardiac Electro-mechanics Simulations <https://zenodo.org/records/3890034>`_
- `Virtual cohort of adult healthy four-chamber heart meshes from CT images <https://zenodo.org/records/4590294>`_

The :attr:`Preprocessor <ansys.health.heart.pre>` module provides methods to extract the necessary information and VTP object from these two databases. You can then process this data further and use it with the :attr:`HeartModel <ansys.health.heart.models.HeartModel>` class.
The module also includes methods to conveniently download data from these two public sources.

>>> from ansys.health.heart.utils.download import download_case_from_zenodo, unpack_case
>>> tar_file = download_case_from_zenodo("Rodero2021", 1, "my-download-dir")
>>> file_path = unpack_case(tar_file)

The ``file_path`` variable contains the path to the downloaded and unpacked CASE or VTK file.

Alternatively, you can provide your own set of input files. Specify a path to the VTP/VTK file that describes the input geometry and a JSON file that describes the parts.

The part definitions JSON file has the following format:

>>> part_definitions = {
...    "Left ventricle": {
...        "id": 1,
...        "enclosed_by_boundaries": {
...            "left-ventricle-endocardium": 1,
...            "left-ventricle-epicardium": 2,
...            "interface_left-ventricle-myocardium_mitral-valve": 3,
...        },
...    }
...}

The ``id`` represents the volumetric part ID, and ``enclosed_by_boundaries`` contains the IDs of the boundaries that enclose the volumetric part. Add a cell data array called ``surface-id`` so that the preprocessor can identify the IDs of the enclosing boundaries. The :attr:`HeartModel <ansys.health.heart.models.HeartModel>` reads the input model and part definitions JSON file.

>>> # Initialize left-ventricular heart model.
>>> model = models.LeftVentricle(working_directory="my-working-dir")
>>> # Load input model and part definitions:
>>> model.load_input(heart, part_definitions, "surface-id")

Use the :attr:`HeartModel <ansys.health.heart.models.HeartModel.mesh_volume>` method to generate the volumetric meshes from the input model.

>>> # Remesh the model using wrapping
>>> model.mesh_volume(use_wrapper=True, global_mesh_size=1.0)

For working examples, see:
- https://heart.health.docs.pyansys.com/version/dev/examples/preprocessor/preprocess_truncated_LV_pr.html#sphx-glr-examples-preprocessor-preprocess-truncated-lv-pr-py
- https://heart.health.docs.pyansys.com/version/dev/examples/preprocessor/doc_preprocess_fullheart_rodero_01.html#sphx-glr-examples-preprocessor-doc-preprocess-fullheart-rodero-01-py