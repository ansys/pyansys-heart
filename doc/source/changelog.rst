.. _ref_release_notes:

Release notes
#############

.. vale off

.. towncrier release notes start

`refs/heads/release/0.11 <https://github.com/ansys/pyansys-heart/releases/tag/vrefs/heads/release/0.11>`_ (April 07, 2025)
==========================================================================================================================

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