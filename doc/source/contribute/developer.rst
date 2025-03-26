Contributing as a developer
###########################

.. grid:: 1 2 3 3

    .. grid-item-card:: :fa:`code-fork` Fork the repository
        :padding: 2 2 2 2
        :link: fork-the-repository
        :link-type: ref

        Learn how to fork the project and get your own copy.

    .. grid-item-card:: :fa:`download` Clone the repository
        :padding: 2 2 2 2
        :link: clone-the-repository
        :link-type: ref

        Download your own copy in your local machine.

    .. grid-item-card:: :fa:`download` Install for developers
        :padding: 2 2 2 2
        :link: install-for-developers
        :link-type: ref

        Install the project in editable mode.

    .. grid-item-card:: :fa:`vial-circle-check` Run the tests
        :padding: 2 2 2 2
        :link: run-tests
        :link-type: ref

        Verify your changes by testing the project.

    .. grid-item-card:: :fa:`computer` Code style compliance
        :padding: 2 2 2 2
        :link: code-style
        :link-type: ref

        Adhere to code style

    .. grid-item-card:: :fa:`arrows-spin` Run the CI/CD pipelines
        :padding: 2 2 2 2
        :link: run-pipelines
        :link-type: ref

        Understand the different CI/CD pipelines.


.. _fork-the-repository:

Fork the repository
===================

Forking the repository is the first step to contributing to the project. This
allows you to have your own copy of the project so you can make changes without
affecting the main project. Once you have made your changes, you can submit a
pull-request to the main project to have your changes reviewed and merged.

.. button-link:: https://github.com/ansys/pyansys-heart/fork
    :color: primary
    :align: center

    :fa:`code-fork` Fork this project

.. note::

    If you are an Ansys employee, you can skip this step.

.. _clone-the-repository:

Clone the repository
====================

Clone the latest version of PyAnsys Heart in development mode by running this code:

.. code-block:: bash

    git clone https://github.com/pyansys/pyansys-heart

.. note::

    If you are not an Ansys employee, you need to :ref:`fork the repository <fork-the-repository>` and
    replace ``ansys`` with your GitHub user name in the ``git clone``
    command.

.. _install-for-developers:

Install for developers
======================

Installing PyAnsys Heart in development mode allows you to perform changes to the code
and see the changes reflected in your environment without having to reinstall
the library every time you make a change.

Virtual environment
-------------------

Start by navigating to the project's root directory by running:

.. code-block::

    cd pyansys-heart

Then, create a new virtual environment named ``.venv`` to isolate your system's
Python environment by running:

.. code-block:: text

    python -m venv .venv

Finally, activate this environment by running:

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

Development mode
----------------

Now, install PyAnsys Heart in editable mode by running:

.. code-block:: text

    python -m pip install --editable .

Verify the installation by checking the version of the library:


.. code-block:: python

    from ansys.heart import __version__


    print(f"PyAnsys Heart version is {__version__}")

.. jinja::

    .. code-block:: text

       >>> PyAnsys Heart version is {{ PYANSYS_HEART_VERSION }}

Install tox
-----------

Once the project is installed, you can install `tox`_. This is a cross-platform
automation tool. The main advantage of Tox is that it eases routine tasks like project
testing, documentation generation, and wheel building in separate and isolated Python
virtual environments. To install Tox, run:

.. code-block:: text

    python -m pip install tox

Finally, verify the installation by listing all the different environments
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
modifications did not break the project. PyAnsys Heart tests support different markers
to allow testing with/without coverage (and against specific python versions).
These markers are associated with dedicated `Tox`_ environments.

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

    The preceding test commands run all tests, including those that require Fluent (which take longer). For more
    selective testing, ``-- -vv -m "not requires_fluent or (not extract_models)"`` or ``-- -vv -m "requires_fluent"`` can be
    appended to tox testing commands.

    .. code:: bash

      # run quick tests
      python -m tox -e tests312-coverage -- -vv -m "not requires_fluent or (not extract_models)"
      # run tests requiring Fluent
      python -m tox -e tests312-coverage -- -vv -m "requires_fluent"

.. _code-style:

Check code style
================

PyAnsys Heart follows the PEP8 standard as outlined in
`PEP 8 <https://dev.docs.pyansys.com/coding-style/pep8.html>`_ in
the *PyAnsys Developer's Guide* and implements style checking using
`pre-commit <https://pre-commit.com/>`_.

To ensure your code meets minimum code styling standards, run the following tox environment:

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
pull-request, labelling a pull-request, and tagging a commit.

You can label a pull-request to skip certain jobs in the pipeline. Supported
labels are listed in the `PyAnsys Heart labels`_ page.

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Label
      - Description
    * - ``test:skip``
      - Skip the model generation tests