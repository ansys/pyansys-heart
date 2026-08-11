# PyAnsys Heart

> [!WARNING]
>
> **This repository is no longer maintained as of 2026-08-12**
>
> This public repository is deprecated and will not receive further public
> updates. Ongoing development continues within Synopsis as part of
> next-generation cardiac simulation solutions.
>
> **What this means**
>
> - No further updates, bug fixes, or pull request reviews
> - Issues will be closed
> - The repository will be archived
>
> For questions, contact [pyansys-core@synopsis.com](mailto:pyansys-core@synopsis.com).
>
> Thank you to everyone who contributed, used, or supported this project! For further questions please reach out to pyansys-core@synopsis.com.

[![PyAnsys](https://img.shields.io/badge/Py-Ansys-ffc107.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAIAAACQkWg2AAABDklEQVQ4jWNgoDfg5mD8vE7q/3bpVyskbW0sMRUwofHD7Dh5OBkZGBgW7/3W2tZpa2tLQEOyOzeEsfumlK2tbVpaGj4N6jIs1lpsDAwMJ278sveMY2BgCA0NFRISwqkhyQ1q/Nyd3zg4OBgYGNjZ2ePi4rB5loGBhZnhxTLJ/9ulv26Q4uVk1NXV/f///////69du4Zdg78lx//t0v+3S88rFISInD59GqIH2esIJ8G9O2/XVwhjzpw5EAam1xkkBJn/bJX+v1365hxxuCAfH9+3b9/+////48cPuNehNsS7cDEzMTAwMMzb+Q2u4dOnT2vWrMHu9ZtzxP9vl/69RVpCkBlZ3N7enoDXBwEAAA+YYitOilMVAAAAAElFTkSuQmCC)](https://docs.pyansys.com/)
[![Python](https://img.shields.io/pypi/pyversions/ansys-health-heart?logo=pypi)](https://pypi.org/project/ansys-health-heart/)
[![PyPI](https://img.shields.io/pypi/v/ansys-health-heart.svg?logo=python&logoColor=white&label=PyPI)](https://pypi.org/project/ansys-health-heart)
[![GH-CI](https://github.com/ansys/pyansys-heart/actions/workflows/ci_cd_night.yml/badge.svg)](https://github.com/ansys/pyansys-heart/actions/workflows/ci_cd_night.yml)
[![MIT](https://img.shields.io/badge/license-MIT-yellow)](https://opensource.org/blog/license/mit)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![pre-commit.ci](https://results.pre-commit.ci/badge/github/ansys/pyansys-heart/main.svg)](https://results.pre-commit.ci/latest/github/ansys/pyansys-heart/main)

## About

PyAnsys Heart is a Python library developed by the Ansys Healthcare Research team to model the human heart using Ansys LS-DYNA software. Designed to support advancements in cardiovascular research, this tool enables detailed, patient-specific simulations of heart function, capturing complex interactions like fluid-structure-electrophysiology interaction (FSEI). These simulations replicate hard-to-measure features such as structural stress, electrical activity, muscle fiber orientation, and even account for anatomical elements like the pericardium, valves, and atria.

## Installation

Ensure you have all the necessary [prerequisites](https://heart.health.docs.pyansys.com/version/stable/getting-started/prerequisites.html) before installing this software. Then, refer to the [installation guidelines](https://heart.health.docs.pyansys.com/version/stable/getting-started/installation.html) for detailed instructions on how to install PyAnsys Heart in your system.

## Documentation

Documentation for the latest stable release of PyAnsys Heart is hosted at [PyAnsys Heart documentation](https://heart.health.docs.pyansys.com/version/stable/).

In the upper right corner of the documentation's title bar, there is an option for switching from viewing the documentation for the latest stable release to viewing the documentation for the development version or previously released versions.

Descriptions follow for each documentation section:

- [Getting started](https://heart.health.docs.pyansys.com/version/stable/getting-started.html): Provides an overview of key techniques in cardiac modeling, package prerequisites, and installation information.
- [User guide](https://heart.health.docs.pyansys.com/version/stable/user-guide.html): Provides an overview of the capabilities of PyAnsys Heart, explaining the key concept of preprocessor, simulator, and postprocessor.
- [API reference](https://heart.health.docs.pyansys.com/version/stable/api/index.html): Describes PyAnsys Heart endpoints, their capabilities, and how to interact with them programmatically.
- [Examples](https://heart.health.docs.pyansys.com/version/stable/examples/index.html): Shows how to use the Preprocessor, Postprocessor, and Simulator modules to preprocess, consume, and analyze heart models.
- [Contribute](https://heart.health.docs.pyansys.com/version/stable/contributing.html): Provides guidelines and instructions on how to contribute based on your role in the project.
- [Release notes](https://heart.health.docs.pyansys.com/version/stable/changelog.html): Provides a summary of notable changes for each version of PyAnsys Heart. It helps you keep track of updates, bug fixes, new features, and improvements made to the project over time.

On the [PyAnsys Heart Issues](https://github.com/ansys/pyansys-heart/issues) page, you can create issues to report bugs and request new features. On the [PyAnsys Heart Discussions](https://github.com/ansys/pyansys-heart/discussions) page or the [Discussions](https://discuss.ansys.com/) page on the Ansys Developer portal, you can post questions, share ideas, and get community feedback. To reach the project support team, email [pyansys-core@synopsis.com](mailto:pyansys-core@synopsis.com).
