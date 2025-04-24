
.. _ref_preprocessor:

Pre-processor
=============

This topic provides an overview of the :attr:`Preprocessor <ansys.health.heart.pre>` module. This module can be used to pre-process a case from the following two public datasets:

- `A Publicly Available Virtual Cohort of Four-chamber Heart Meshes for Cardiac Electro-mechanics Simulations <https://zenodo.org/records/3890034>`_
- `Virtual cohort of adult healthy four-chamber heart meshes from CT images <https://zenodo.org/records/4590294>`_

The :attr:`Pre-processor <ansys.health.heart.pre>` module contains methods to extract the necessary information and VTP object from these two databases for further processing and consumption by the :attr:`HeartModel <ansys.health.heart.models.HeartModel>` class. This module also contains methods to conveniently download from these two sources:

>>> from ansys.health.heart.utils.download import download_case_from_zenodo, unpack_case
>>> tar_file = download_case_from_zenodo("Rodero2021", 1, "my-download-dir")
>>> file_path = unpack_case(tar_file)

Here ``file_path`` gives you the path to the downloaded and unpacked CASE or VTK file.

Alternatively, you can provide your own set of input files. In this case, you must specify a path to the VTP/VTK that describes the input geometry and a JSON file that describes the parts.

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

Here ``id`` is the volumetric part ID, and ``enclosed_by_boundaries`` is a dictionary that contains the IDs of boundaries
that enclose the volumetric part. The ID's of the enclosing boundaries should be identifiable from a cell data array,
for instance by adding a cell data array called ``surface-id``. Consequently this input model and part definitions JSON file
can be read into by the main class :attr:`HeartModel <ansys.health.heart.models.HeartModel>`.

>>> # Initialize left-ventricular heart model.
>>> model = models.LeftVentricle(working_directory="my-working-dir")
>>> # Load input model and part definitions:
>>> model.load_input(heart, part_definitions, "surface-id")

Finally the :attr:`HeartModel <ansys.health.heart.models.HeartModel.mesh_volume>` method can be used to generate the volumetric meshes from the input model.

>>> # remesh the model using wrapping
>>> model.mesh_volume(use_wrapper=True, global_mesh_size=1.0)

For working examples see https://heart.health.docs.pyansys.com/version/dev/examples/preprocessor/preprocess_truncated_LV_pr.html#sphx-glr-examples-preprocessor-preprocess-truncated-lv-pr-py and https://heart.health.docs.pyansys.com/version/dev/examples/preprocessor/doc_preprocess_fullheart_rodero_01.html#sphx-glr-examples-preprocessor-doc-preprocess-fullheart-rodero-01-py