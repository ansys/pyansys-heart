Contribute as a developer
#########################

.. grid:: 1 2 3 3

    .. grid-item-card:: :fa:`code-fork` Fork the repository
        :padding: 2 2 2 2
        :link: fork-the-repository
        :link-type: ref

        Fork the project to create a copy.

    .. grid-item-card:: :fa:`download` Clone the repository
        :padding: 2 2 2 2
        :link: clone-the-repository
        :link-type: ref

        Clone the repository to download the copy to your local machine.

    .. grid-item-card:: :fa:`download` Install for developers
        :padding: 2 2 2 2
        :link: install-for-developers
        :link-type: ref

        Install the project in editable mode.

    .. grid-item-card:: :fa:`vial-circle-check` Run the tests
        :padding: 2 2 2 2
        :link: run-tests
        :link-type: ref

        Verify your changes to the project by running tests.

    .. grid-item-card:: :fa:`computer` Code style compliance
        :padding: 2 2 2 2
        :link: code-style
        :link-type: ref

        Adhere to code style.

    .. grid-item-card:: :fa:`arrows-spin` Run the CI/CD pipelines
        :padding: 2 2 2 2
        :link: run-pipelines
        :link-type: ref

        Understand the different CI/CD pipelines that are executed
        automatically.


.. _fork-the-repository:

Fork the repository
===================

Forking the repository is the first step to contributing to the project. This
allows you to have your own copy of the project so that you can make changes without
affecting the main project. Once you have made your changes, you can submit a
pull request to the main project to have your changes reviewed and merged.

.. button-link:: https://github.com/ansys/pyansys-heart/fork
    :color: primary
    :align: left

    :fa:`code-fork` Fork this project

.. note::

    If you are an Ansys employee, you can skip this step.

.. _clone-the-repository:

Clone the repository
====================

Clone the repository in development mode:

.. code-block:: bash

    git clone https://github.com/pyansys/pyansys-heart

.. note::

    If you are not an Ansys employee, you must :ref:`fork the repository <fork-the-repository>` and
    replace ``ansys`` with your GitHub user name in the ``git clone`` command.

.. _install-for-developers:

Install for developers
======================

Installing PyAnsys Heart in development mode lets you change the code
and see these changes reflected in your environment without having to reinstall
the library every time you make a change.

Set up a virtual environment
----------------------------

#. Navigate to the project's root directory :

.. code-block::

       cd pyansys-heart

#. Create a virtual environment named ``.venv`` to isolate your Python environment:

.. code-block:: text

    python -m venv .venv

#. Activate the virtual environment:

.. tab-set::

    .. tab-item:: Windows

        .. tab-set::

            .. tab-item:: CMD

                .. code-block:: text

                    .venv\Scripts\activate.bat

            .. tab-item:: PowerShell

                .. code-block:: text

                    .venv\Scripts\Activate.ps1

    .. tab-item:: macOS/Linux/UNIX

        .. code-block:: text

            source .venv/bin/activate

Install in development mode
---------------------------

#. Install PyAnsys Heart in editable mode:

   .. code-block:: text

       python -m pip install --editable .

#. Verify the installation by checking the version of the library:

.. code-block:: python

    from ansys.health.heart import __version__

    print(f"PyAnsys Heart version is {__version__}")

.. jinja::

    .. code-block:: text

       >>> PyAnsys Heart version is {{ PYANSYS_HEART_VERSION }}.

Install Tox
-----------

Once the project is installed, you can install `Tox`_. This is a cross-platform
automation tool. The main advantage of Tox is that it eases routine tasks like project
testing, documentation generation, and wheel building in separate and isolated Python
virtual environments.

#. Install Tox:

   .. code-block:: text

       python -m pip install tox

#. Verify the installation by listing all the different environments
(automation rules) for PyAnsys Heart:

   .. code-block:: text

       python -m tox list

   .. jinja:: toxenvs

       .. dropdown:: Default Tox environments
           :animate: fade-in
           :icon: three-bars

           .. list-table::
               :header-rows: 1
               :widths: auto

               * - Environment
                 - Description
               {% for environment in envs %}
               {% set name, description  = environment.split("->") %}
               * - {{ name }}
                 - {{ description }}
               {% endfor %}

.. _run-tests:

Run the tests
=============

Once you have made your changes, you can run the tests to verify that your
changes did not break the project. PyAnsys Heart tests support different markers
to allow testing with or without coverage (and against specific Python versions).
These markers are associated with dedicated Tox environments.

.. jinja:: toxenvs

    .. dropdown:: Testing environments
        :animate: fade-in
        :icon: three-bars

        .. list-table::
            :header-rows: 1
            :widths: auto

            * - Environment
              - Command
            {% for environment in envs %}
            {% set name, description  = environment.split("->") %}
            {% if name.startswith("tests")%}
            * - {{ name }}
              - python -m tox -e {{ name }}
            {% endif %}
            {% endfor %}

.. Note::

    The preceding test code runs all tests, including those that require Fluent (which take longer). For more
    selective testing, append ``-- -vv -m "not requires_fluent or (not extract_models)"`` or ``-- -vv -m "requires_fluent"``
    to Tox testing commands:

    .. code:: bash

      # run quick tests
      python -m tox -e tests312-coverage -- -vv -m "not requires_fluent or (not extract_models)"
      # run tests requiring Fluent
      python -m tox -e tests312-coverage -- -vv -m "requires_fluent"

.. _code-style:

Check code style
================

PyAnsys Heart follows the PEP 8 standard as described in
`PEP 8 <https://dev.docs.pyansys.com/coding-style/pep8.html>`_ in
the *PyAnsys developer's guide* and implements style checking using
`pre-commit <https://pre-commit.com/>`_.

To ensure your code meets minimum code styling standards, run the following Tox environment:

.. jinja:: toxenvs

    .. dropdown:: Code style environment
        :animate: fade-in
        :icon: three-bars

        .. list-table::
            :header-rows: 1
            :widths: auto

            * - Environment
              - Command
            {% for environment in envs %}
            {% set name, description  = environment.split("->") %}
            {% if name.startswith("code-")%}
            * - {{ name }}
              - python -m tox -e {{ name }}
            {% endif %}
            {% endfor %}

This way, it's not possible for you to push code that fails the style checks::

  $ git commit -am "added my cool feature"
  black....................................................................Passed
  blacken-docs.............................................................Passed
  isort....................................................................Passed
  flake8...................................................................Passed
  codespell................................................................Passed
  pydocstyle...............................................................Passed
  check for merge conflicts................................................Passed
  debug statements (python)................................................Passed
  check yaml...............................................................Passed
  trim trailing whitespace.................................................Passed
  Validate GitHub Workflows................................................Passed

.. _run-pipelines:

Run CI/CD pipelines
===================

PyAnsys Heart has a set of CI/CD pipelines that are executed automatically when certain
events are detected in the repository. Some of these events include opening a
pull request, labeling a pull-request, and tagging a commit.

You can label a pull request to skip certain jobs in the pipeline. Supported
labels are listed on the `PyAnsys Heart labels`_ page.

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Label
      - Description
    * - ``test:skip``
      - Skip the model generation tests
