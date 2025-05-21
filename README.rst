PyAnsys Heart
#############

|pyansys| |python| |pypi| |GH-CI| |MIT| |ruff| |pre-commit|

.. |pyansys| image:: https://img.shields.io/badge/Py-Ansys-ffc107.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAIAAACQkWg2AAABDklEQVQ4jWNgoDfg5mD8vE7q/3bpVyskbW0sMRUwofHD7Dh5OBkZGBgW7/3W2tZpa2tLQEOyOzeEsfumlK2tbVpaGj4N6jIs1lpsDAwMJ278sveMY2BgCA0NFRISwqkhyQ1q/Nyd3zg4OBgYGNjZ2ePi4rB5loGBhZnhxTLJ/9ulv26Q4uVk1NXV/f///////69du4Zdg78lx//t0v+3S88rFISInD59GqIH2esIJ8G9O2/XVwhjzpw5EAam1xkkBJn/bJX+v1365hxxuCAfH9+3b9/+////48cPuNehNsS7cDEzMTAwMMzb+Q2u4dOnT2vWrMHu9ZtzxP9vl/69RVpCkBlZ3N7enoDXBwEAAA+YYitOilMVAAAAAElFTkSuQmCC
   :target: https://docs.pyansys.com/
   :alt: PyAnsys

.. |python| image:: https://img.shields.io/pypi/pyversions/ansys-health-heart?logo=pypi
   :target: https://pypi.org/project/ansys-health-heart/
   :alt: Python

.. |pypi| image:: https://img.shields.io/pypi/v/ansys-health-heart.svg?logo=python&logoColor=white&label=PyPI
   :target: https://pypi.org/project/ansys-health-heart
   :alt: PyPI

.. |GH-CI| image:: https://github.com/ansys/pyansys-heart/actions/workflows/ci_cd_night.yml/badge.svg
   :target: https://github.com/ansys/pyansys-heart/actions/workflows/ci_cd_night.yml
   :alt: GH-CI

.. |MIT| image:: https://img.shields.io/badge/license-MIT-yellow
   :target: https://opensource.org/blog/license/mit
   :alt: MIT

.. |ruff| image:: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
   :target: https://github.com/astral-sh/ruff
   :alt: Ruff

.. |pre-commit| image:: https://results.pre-commit.ci/badge/github/ansys/pyansys-heart/main.svg
   :target: https://results.pre-commit.ci/latest/github/ansys/pyansys-heart/main
   :alt: pre-commit.ci

About
=====

PyAnsys Heart is a Python library developed by the Ansys Healthcare Research
team to model the human heart using Ansys LS-DYNA software. Designed to support
advancements in cardiovascular research, this tool enables detailed,
patient-specific simulations of heart function, capturing complex interactions
like fluid-structure-electrophysiology interaction (FSEI). These simulations
replicate hard-to-measure features such as structural stress, electrical
activity, muscle fiber orientation, and even account for anatomical elements
like the pericardium, valves, and atria.

Installation
============

Ensure you have all the necessary `prerequisites`_ before installing this
software. Then, refer to the `installation guidelines`_ for detailed instructions
on how to install PyAnsys Heart in your system.

Documentation
=============

Documentation for the latest stable release of PyAnsys Heart is hosted at
`PyAnsys Heart documentation`_.

In the upper right corner of the documentation's title bar, there is an option
for switching from viewing the documentation for the latest stable release
to viewing the documentation for the development version or previously
released versions.

Descriptions follow for each documentation section:

- `Getting started`_: Provides an overview of key techniques in cardiac modeling,
  package prerequisites, and installation information.

- `User guide`_: Provides an overview of the capabilities of PyAnsys Heart,
  explaining the key concept of preprocessor, simulator, and postprocessor.

- `API reference`_: Describes PyAnsys Heart endpoints, their capabilities,
  and how to interact with them programmatically.

- `Examples`_: Shows how to use the Preprocessor, Postprocessor, and Simulator
  modules to preprocess, consume, and analyze heart models.

- `Contribute`_: Provides guidelines and instructions on how to contribute
  based on your role in the project.

- `Release notes`_: Provides a summary of notable changes for each version of
  PyAnsys Heart. It helps you keep track of updates, bug fixes, new features, and
  improvements made to the project over time.

On the `PyAnsys Heart Issues <https://github.com/ansys/pyansys-heart/issues>`_ page,
you can create issues to report bugs and request new features. On the
`PyAnsys Heart Discussions <https://github.com/ansys/pyansys-heart/discussions>`_ page
or the `Discussions <https://discuss.ansys.com/>`_ page on the Ansys Developer
portal, you can post questions, share ideas, and get community feedback.
To reach the project support team, email `pyansys.core@ansys.com <mailto:pyansys.core@ansys.com>`_.


.. Documentation links
.. _prerequisites: https://heart.health.docs.pyansys.com/version/stable/getting-started/prerequisites.html
.. _installation guidelines: https://heart.health.docs.pyansys.com/version/stable/getting-started/installation.html
.. _getting started: https://heart.health.docs.pyansys.com/version/stable/getting-started.html
.. _user guide: https://heart.health.docs.pyansys.com/version/stable/user-guide.html
.. _API reference: https://heart.health.docs.pyansys.com/version/stable/api/index.html
.. _examples: https://heart.health.docs.pyansys.com/version/stable/examples/index.html
.. _contribute: https://heart.health.docs.pyansys.com/version/stable/contributing.html
.. _LICENSE: https://github.com/ansys/pyansys-heart/blob/main/LICENSE
.. _release notes: https://heart.health.docs.pyansys.com/version/stable/changelog.html
.. _PyAnsys Heart documentation: https://heart.health.docs.pyansys.com/version/stable/

