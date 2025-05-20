.. _ref_release_notes:

Release notes
#############

.. vale off

.. towncrier release notes start

`0.13.1 <https://github.com/ansys/pyansys-heart/releases/tag/v0.13.1>`_ (May 20, 2025)
======================================================================================

.. tab-set::


  .. tab-item:: Fixed

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - bump ansys-sphinx-theme
          - `#1102 <https://github.com/ansys/pyansys-heart/pull/1102>`_


  .. tab-item:: Documentation

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - update contributing.rst
          - `#1105 <https://github.com/ansys/pyansys-heart/pull/1105>`_

        * - add acknowledgments
          - `#1110 <https://github.com/ansys/pyansys-heart/pull/1110>`_

        * - fix links in README
          - `#1116 <https://github.com/ansys/pyansys-heart/pull/1116>`_

        * - update install instructions
          - `#1117 <https://github.com/ansys/pyansys-heart/pull/1117>`_


  .. tab-item:: Dependencies

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - update scipy requirement from <=1.15.2 to <=1.15.3
          - `#1108 <https://github.com/ansys/pyansys-heart/pull/1108>`_


  .. tab-item:: Maintenance

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - bump the actions group across 1 directory with 3 updates
          - `#1088 <https://github.com/ansys/pyansys-heart/pull/1088>`_

        * - bump version to 0.13.dev0
          - `#1096 <https://github.com/ansys/pyansys-heart/pull/1096>`_

        * - update CHANGELOG for v0.12.1
          - `#1100 <https://github.com/ansys/pyansys-heart/pull/1100>`_

        * - update CHANGELOG for v0.12.2
          - `#1104 <https://github.com/ansys/pyansys-heart/pull/1104>`_

        * - release to public pypi
          - `#1112 <https://github.com/ansys/pyansys-heart/pull/1112>`_


  .. tab-item:: Miscellaneous

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - refactor dyna writer to improve maintainability
          - `#1101 <https://github.com/ansys/pyansys-heart/pull/1101>`_

        * - refactor and further cleanup for release
          - `#1109 <https://github.com/ansys/pyansys-heart/pull/1109>`_


`0.12.2 <https://github.com/ansys/pyansys-heart/releases/tag/v0.12.2>`_ (May 08, 2025)
======================================================================================

.. tab-set::


  .. tab-item:: Added

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - closed system
          - `#1021 <https://github.com/ansys/pyansys-heart/pull/1021>`_

        * - technical review
          - `#1037 <https://github.com/ansys/pyansys-heart/pull/1037>`_

        * - improve-EP-default-conduction
          - `#1069 <https://github.com/ansys/pyansys-heart/pull/1069>`_

        * - allow passing additional keyword arguments to launch_fluent
          - `#1095 <https://github.com/ansys/pyansys-heart/pull/1095>`_


  .. tab-item:: Fixed

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - modify beam mesh doc strings
          - `#1030 <https://github.com/ansys/pyansys-heart/pull/1030>`_

        * - remove wheelhouse from doc/source/_static dir
          - `#1032 <https://github.com/ansys/pyansys-heart/pull/1032>`_

        * - follow ansys.health namespace
          - `#1036 <https://github.com/ansys/pyansys-heart/pull/1036>`_

        * - documentation build
          - `#1039 <https://github.com/ansys/pyansys-heart/pull/1039>`_

        * - run examples in pipelines
          - `#1040 <https://github.com/ansys/pyansys-heart/pull/1040>`_

        * - fall back to mpiexec when mpirun is not found
          - `#1050 <https://github.com/ansys/pyansys-heart/pull/1050>`_

        * - avoid pyvista 0.45
          - `#1060 <https://github.com/ansys/pyansys-heart/pull/1060>`_

        * - adding EM_CONTROL_TIMESTEP to fiber generation decks
          - `#1071 <https://github.com/ansys/pyansys-heart/pull/1071>`_

        * - mutable default
          - `#1072 <https://github.com/ansys/pyansys-heart/pull/1072>`_

        * - force update node mesh ID in laplacewriter
          - `#1079 <https://github.com/ansys/pyansys-heart/pull/1079>`_

        * - convert int64 data to int32 for visualization
          - `#1097 <https://github.com/ansys/pyansys-heart/pull/1097>`_

        * - changelog actions version in release ci
          - `#1099 <https://github.com/ansys/pyansys-heart/pull/1099>`_


  .. tab-item:: Documentation

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - overall review
          - `#1043 <https://github.com/ansys/pyansys-heart/pull/1043>`_

        * - link to documentation is broken
          - `#1046 <https://github.com/ansys/pyansys-heart/pull/1046>`_

        * - execute all examples nightly doc build
          - `#1054 <https://github.com/ansys/pyansys-heart/pull/1054>`_

        * - update atrial fiber example
          - `#1064 <https://github.com/ansys/pyansys-heart/pull/1064>`_

        * - update user guide and expose pre, post and simulator api docs
          - `#1065 <https://github.com/ansys/pyansys-heart/pull/1065>`_

        * - interactive plots in examples
          - `#1073 <https://github.com/ansys/pyansys-heart/pull/1073>`_

        * - edits based on skimming rendered doc
          - `#1075 <https://github.com/ansys/pyansys-heart/pull/1075>`_

        * - add left ventricle mechanical example
          - `#1076 <https://github.com/ansys/pyansys-heart/pull/1076>`_

        * - add basic ep postprocessor example
          - `#1080 <https://github.com/ansys/pyansys-heart/pull/1080>`_

        * - fix interactive plots in doc build
          - `#1086 <https://github.com/ansys/pyansys-heart/pull/1086>`_

        * - cleanup and fixes for examples
          - `#1087 <https://github.com/ansys/pyansys-heart/pull/1087>`_

        * - switch to ReactionEikonal for ep-mechanics example
          - `#1090 <https://github.com/ansys/pyansys-heart/pull/1090>`_

        * - reduce size of vtksz for doc build
          - `#1091 <https://github.com/ansys/pyansys-heart/pull/1091>`_


  .. tab-item:: Dependencies

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - update flit-core requirement from <3.11,>=3.2 to >=3.2,<4
          - `#1025 <https://github.com/ansys/pyansys-heart/pull/1025>`_

        * - bump pytest-cov from 6.0.0 to 6.1.1
          - `#1026 <https://github.com/ansys/pyansys-heart/pull/1026>`_

        * - update numpy requirement from <=2.2.4 to <=2.2.5
          - `#1059 <https://github.com/ansys/pyansys-heart/pull/1059>`_


  .. tab-item:: Maintenance

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - update CHANGELOG for v0.11.0
          - `#1023 <https://github.com/ansys/pyansys-heart/pull/1023>`_

        * - bump version to 0.12.dev0
          - `#1033 <https://github.com/ansys/pyansys-heart/pull/1033>`_

        * - bump the actions group across 1 directory with 4 updates
          - `#1034 <https://github.com/ansys/pyansys-heart/pull/1034>`_

        * - bump ansys/actions from 9.0.0 to 9.0.2 in the actions group
          - `#1048 <https://github.com/ansys/pyansys-heart/pull/1048>`_

        * - use intelmpi on runner for doc build
          - `#1061 <https://github.com/ansys/pyansys-heart/pull/1061>`_

        * - update nightly and release doc builds
          - `#1070 <https://github.com/ansys/pyansys-heart/pull/1070>`_

        * - only run release workflow on tag push
          - `#1098 <https://github.com/ansys/pyansys-heart/pull/1098>`_


  .. tab-item:: Miscellaneous

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - standardize type hints for ``pre``, ``post``, and ``utils`` subpackages
          - `#1018 <https://github.com/ansys/pyansys-heart/pull/1018>`_

        * - remove unused and outdated method
          - `#1035 <https://github.com/ansys/pyansys-heart/pull/1035>`_

        * - improve how conduction paths and their data are managed
          - `#1041 <https://github.com/ansys/pyansys-heart/pull/1041>`_

        * - consolidate _BeamsMesh functionality into Mesh
          - `#1042 <https://github.com/ansys/pyansys-heart/pull/1042>`_

        * - only print LS-DYNA stdout to debug level
          - `#1081 <https://github.com/ansys/pyansys-heart/pull/1081>`_

        * - deprecate update parts
          - `#1089 <https://github.com/ansys/pyansys-heart/pull/1089>`_


`0.11.0 <https://github.com/ansys/pyansys-heart/releases/tag/v0.11.0>`_ (April 07, 2025)
========================================================================================

.. tab-set::


  .. tab-item:: Added

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - add changelog actions and changelog documentation
          - `#908 <https://github.com/ansys/pyansys-heart/pull/908>`_

        * - handle incompressibility consistently
          - `#909 <https://github.com/ansys/pyansys-heart/pull/909>`_

        * - refactor-beam-networks
          - `#932 <https://github.com/ansys/pyansys-heart/pull/932>`_

        * - add D-RBM method for left ventricle model
          - `#933 <https://github.com/ansys/pyansys-heart/pull/933>`_

        * - compute ventricle thickening
          - `#945 <https://github.com/ansys/pyansys-heart/pull/945>`_

        * - set stiffness damping
          - `#980 <https://github.com/ansys/pyansys-heart/pull/980>`_

        * - add module for custom exceptions
          - `#990 <https://github.com/ansys/pyansys-heart/pull/990>`_

        * - Append user k files
          - `#992 <https://github.com/ansys/pyansys-heart/pull/992>`_


  .. tab-item:: Fixed

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - add-EMCONTROLTIMESTEP-in-ep
          - `#922 <https://github.com/ansys/pyansys-heart/pull/922>`_

        * - fix cap types and cap type check
          - `#935 <https://github.com/ansys/pyansys-heart/pull/935>`_

        * - refactor part id assignment post wrap
          - `#946 <https://github.com/ansys/pyansys-heart/pull/946>`_

        * - syntax error
          - `#950 <https://github.com/ansys/pyansys-heart/pull/950>`_

        * - tox file correction and improvement
          - `#956 <https://github.com/ansys/pyansys-heart/pull/956>`_

        * - `test_ep_postprocessor` tests on Github runner
          - `#971 <https://github.com/ansys/pyansys-heart/pull/971>`_

        * - reassign part ids when no orphan cells are found
          - `#983 <https://github.com/ansys/pyansys-heart/pull/983>`_

        * - shutil.which for wsl
          - `#995 <https://github.com/ansys/pyansys-heart/pull/995>`_

        * - pinned versions for direct dependencies
          - `#996 <https://github.com/ansys/pyansys-heart/pull/996>`_


  .. tab-item:: Documentation

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - Cleanup
          - `#923 <https://github.com/ansys/pyansys-heart/pull/923>`_

        * - add the landing page
          - `#949 <https://github.com/ansys/pyansys-heart/pull/949>`_

        * - refactor user guide and getting started
          - `#955 <https://github.com/ansys/pyansys-heart/pull/955>`_

        * - contributing guide improvement
          - `#961 <https://github.com/ansys/pyansys-heart/pull/961>`_

        * - update docstrings and standardize periods
          - `#991 <https://github.com/ansys/pyansys-heart/pull/991>`_


  .. tab-item:: Dependencies

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - bump tox from 4.24.1 to 4.24.2
          - `#910 <https://github.com/ansys/pyansys-heart/pull/910>`_

        * - bump ansys-dpf-core from 0.13.4 to 0.13.6
          - `#912 <https://github.com/ansys/pyansys-heart/pull/912>`_

        * - cleanup dependencies list
          - `#913 <https://github.com/ansys/pyansys-heart/pull/913>`_

        * - bump ansys-fluent-core from 0.29.0 to 0.30.0
          - `#940 <https://github.com/ansys/pyansys-heart/pull/940>`_

        * - update numpy requirement from <=2.2.3 to <=2.2.4
          - `#941 <https://github.com/ansys/pyansys-heart/pull/941>`_

        * - bump the docs-deps group across 1 directory with 2 updates
          - `#954 <https://github.com/ansys/pyansys-heart/pull/954>`_


  .. tab-item:: Maintenance

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - self hosted runner
          - `#904 <https://github.com/ansys/pyansys-heart/pull/904>`_

        * - workflow improvements
          - `#951 <https://github.com/ansys/pyansys-heart/pull/951>`_

        * - mark and cleanup tests that require dpf
          - `#981 <https://github.com/ansys/pyansys-heart/pull/981>`_

        * - release to private pypi
          - `#1019 <https://github.com/ansys/pyansys-heart/pull/1019>`_


  .. tab-item:: Miscellaneous

    .. list-table::
        :header-rows: 0
        :widths: auto

        * - clean up deprecated dump model
          - `#914 <https://github.com/ansys/pyansys-heart/pull/914>`_

        * - volume meshing and mesher module
          - `#915 <https://github.com/ansys/pyansys-heart/pull/915>`_

        * - name of Material 295
          - `#918 <https://github.com/ansys/pyansys-heart/pull/918>`_

        * - cleanup and introduce new environment variables to manage automation
          - `#919 <https://github.com/ansys/pyansys-heart/pull/919>`_

        * - volume meshing and mesher module (#915)
          - `#921 <https://github.com/ansys/pyansys-heart/pull/921>`_

        * - create misc module
          - `#924 <https://github.com/ansys/pyansys-heart/pull/924>`_

        * - rename landmarks module to landmark_utils
          - `#927 <https://github.com/ansys/pyansys-heart/pull/927>`_

        * - move slerp methods to misc
          - `#930 <https://github.com/ansys/pyansys-heart/pull/930>`_

        * - download module
          - `#934 <https://github.com/ansys/pyansys-heart/pull/934>`_

        * - rename custom keywords and keywords_module
          - `#936 <https://github.com/ansys/pyansys-heart/pull/936>`_

        * - uhcwriter
          - `#937 <https://github.com/ansys/pyansys-heart/pull/937>`_

        * - rename vtkmethods to vtk_utils
          - `#938 <https://github.com/ansys/pyansys-heart/pull/938>`_

        * - cleanup paths in examples
          - `#943 <https://github.com/ansys/pyansys-heart/pull/943>`_

        * - mecha writer clean up
          - `#944 <https://github.com/ansys/pyansys-heart/pull/944>`_

        * - add method to get fluent ui-mode
          - `#957 <https://github.com/ansys/pyansys-heart/pull/957>`_

        * - move symbols to dpf utils and cleanup
          - `#960 <https://github.com/ansys/pyansys-heart/pull/960>`_

        * - replace wget by httpx
          - `#962 <https://github.com/ansys/pyansys-heart/pull/962>`_

        * - cleanup and refactor preprocessor module
          - `#969 <https://github.com/ansys/pyansys-heart/pull/969>`_

        * - rename helpers subpackage and downloader module
          - `#970 <https://github.com/ansys/pyansys-heart/pull/970>`_

        * - dynain file in mechanical simulator
          - `#977 <https://github.com/ansys/pyansys-heart/pull/977>`_

        * - boundary type and anatomy axis exception
          - `#988 <https://github.com/ansys/pyansys-heart/pull/988>`_

        * - remove deprecated arguments and methods
          - `#998 <https://github.com/ansys/pyansys-heart/pull/998>`_

        * - move packages to core
          - `#1014 <https://github.com/ansys/pyansys-heart/pull/1014>`_

        * - change structure of tests
          - `#1017 <https://github.com/ansys/pyansys-heart/pull/1017>`_


.. vale on