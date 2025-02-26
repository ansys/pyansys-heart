============
Contributing
============

Overall guidance on contributing to a PyAnsys repository appears in the *PyAnsys Developer's Guide*, see:
`Contribute <https://dev.docs.pyansys.com/>`_. Ensure that you are thoroughly familiar
with this guide before attempting to contribute to `PyAnsys Heart <https://github.com/ansys/pyansys-heart>`_.

The following contribution information is specific to `PyAnsys Heart <https://github.com/ansys/pyansys-heart>`_.

Clone the repository
--------------------
Clone and install the latest version of PyAnsys Heart in
development mode by running this code:

.. code:: bash

    git clone https://github.com/pyansys/pyansys-heart
    cd pyansys-heart
    pip install -e .


Install additional requirements, such as dependencies to build documentation or run (unit)tests:

.. code:: bash

  # dependencies for local doc building
  python -m pip install -e .[doc]
  # dependencies needed for (unit) testing
  python -m pip install -e .[tests]

Run tests to verify your development version.

.. code:: bash

  # run quick tests
  python -m pytest -v -m "not requires_fluent or (not extract_models)"
  # run tests requiring Fluent
  python -m pytest -v -m requires_fluent
  # run all tests
  pytest tests -v


Post issues
-----------
Use the `PyAnsys Heart Issues <https://github.com/ansys/pyansys-heart/issues>`_
page to submit questions, report bugs, and request new features. When possible, you
should use these issue templates:

* Bug, problem, error: For filing a bug report
* Documentation error: For requesting modifications to the documentation
* Adding an example: For proposing a new example
* New feature: For requesting enhancements to the code

If your issue does not fit into one of these template categories, you can click
the link for opening a blank issue.

To reach the project support team, email `pyansys.core@ansys.com <pyansys.core@ansys.com>`_.

View documentation
------------------
Documentation for the latest stable release of PyAnsys Heart is hosted at
https://heart.docs.pyansys.com/.

In the upper right corner of the documentation's title bar, there is an option
for switching from viewing the documentation for the latest stable release
to viewing the documentation for the development version or previously
released versions.

Adhere to code style
--------------------

PyAnsys Heart follows the PEP8 standard as outlined in
`PEP 8 <https://dev.docs.pyansys.com/coding-style/pep8.html>`_ in
the *PyAnsys Developer's Guide* and implements style checking using
`pre-commit <https://pre-commit.com/>`_.

To ensure your code meets minimum code styling standards, run these commands::

  pip install pre-commit
  pre-commit run --all-files

You can also install this as a pre-commit hook by running this command::

  pre-commit install

This way, it's not possible for you to push code that fails the style checks::

  $ pre-commit install
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

Install tox
-----------
Once the project is installed, you can install `Tox <https://tox.wiki/en/stable/>`_ to run tests in an isolated environment.
Tox is a generic virtual environment management and test command line tool that can be used to test your project
in isolated python environments.

.. code:: bash

  python -m pip install tox

To verify that your project is working as expected, run the following command

.. code:: bash

  tox list

This command will list all the environments that are available for testing.

.. jinja:: toxenvs

    .. dropdown:: Default Tox environments
        :animate: fade-in
        :icon: three-bars

        .. list-table::
            :header-rows: 1
            :widths: auto

            * - Environment
              - Description
              - usage
            {% for environment in envs %}
            {% set name, description  = environment.split("->") %}
            * - {{ name }}
              - {{ description }}
              - python -m tox -e {{ name }}
            {% endfor %}

Run CI/CD pipelines
-------------------
You can label a pull-request to skip certain jobs in the pipeline as follows:

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Label
      - Description
    * - ``test:skip``
      - Skip the model generation tests