# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import os
import sys

import pytest

from ansys.health.heart.settings.settings import DynaSettings

if os.getenv("GITHUB_ACTIONS"):
    is_gh_action = True
else:
    is_gh_action = False
if "win" in sys.platform:
    is_windows = True
else:
    is_windows = False


@pytest.mark.parametrize(
    "dynatype",
    [
        "smp",
        "intelmpi",
        "platformmpi",
        "msmpi",
    ],
)
@pytest.mark.parametrize(
    "platform",
    [
        "windows",
        "linux",
        pytest.param(
            "wsl",
            marks=pytest.mark.xfail(not is_windows, reason="WSL Only valid argument on Windows"),
        ),
    ],
)
def test_get_dyna_commands_001(dynatype, platform):
    """Test if get commands returns right command line if no additional options are given."""
    if dynatype == "msmpi" and platform != "windows":
        pytest.skip("MSMPI and %s are not compatible and does not make sense to test." % platform)

    # define mock data
    settings = DynaSettings(
        lsdyna_path="my-dyna-path.exe",
        dynatype=dynatype,
        num_cpus=2,
        platform=platform,
    )

    commands = settings.get_commands("path-to-input.k")

    if platform == "wsl":
        expected = ["powershell", "-Command", "wsl", "-e", "bash", "-lic", "./run_lsdyna.sh"]
    else:
        if dynatype == "smp":
            expected = ["my-dyna-path.exe", "i=path-to-input.k", "ncpu=2"]
        elif dynatype in ["intelmpi", "platformmpi"]:
            expected = ["mpirun", "-np", "2", "my-dyna-path.exe", "i=path-to-input.k"]
        elif dynatype == "msmpi":
            expected = ["mpiexec", "-np", "2", "my-dyna-path.exe", "i=path-to-input.k"]

    assert commands == expected


@pytest.mark.parametrize(
    "dynatype",
    [
        "smp",
        "intelmpi",
        "platformmpi",
        "msmpi",
    ],
)
@pytest.mark.parametrize(
    "platform",
    [
        "windows",
        "linux",
        pytest.param(
            "wsl",
            marks=pytest.mark.xfail(not is_windows, reason="WSL Only valid argument on Windows"),
        ),
    ],
)
def test_get_dyna_commands_002(dynatype, platform):
    """Test if get commands returns right command line arguments if dyna options are given."""
    if dynatype == "msmpi" and platform != "windows":
        pytest.skip("MSMPI and %s are not compatible and does not make sense to test." % platform)

    # define mock data
    settings = DynaSettings(
        lsdyna_path="my-dyna-path.exe",
        dynatype=dynatype,
        num_cpus=2,
        platform=platform,
        dyna_options="memory=1000",
    )

    commands = settings.get_commands("path-to-input.k")

    if platform == "wsl":
        expected = ["powershell", "-Command", "wsl", "-e", "bash", "-lic", "./run_lsdyna.sh"]
    else:
        if dynatype == "smp":
            expected = ["my-dyna-path.exe", "i=path-to-input.k", "ncpu=2", "memory=1000"]
        elif dynatype in ["intelmpi", "platformmpi"]:
            expected = [
                "mpirun",
                "-np",
                "2",
                "my-dyna-path.exe",
                "i=path-to-input.k",
                "memory=1000",
            ]
        elif dynatype == "msmpi":
            expected = [
                "mpiexec",
                "-np",
                "2",
                "my-dyna-path.exe",
                "i=path-to-input.k",
                "memory=1000",
            ]

    assert commands == expected


@pytest.mark.parametrize(
    "dynatype",
    [
        "smp",
        "intelmpi",
        "platformmpi",
        "msmpi",
    ],
)
@pytest.mark.parametrize(
    "platform",
    [
        "windows",
        "linux",
        pytest.param(
            "wsl",
            marks=pytest.mark.xfail(not is_windows, reason="WSL Only valid argument on Windows"),
        ),
    ],
)
def test_get_dyna_commands_003(dynatype, platform):
    """Test if get commands returns right command line arguments if dyna and mpi options present."""
    if dynatype == "msmpi" and platform != "windows":
        pytest.skip("MSMPI and %s are not compatible and does not make sense to test." % platform)

    # define mock data
    settings = DynaSettings(
        lsdyna_path="my-dyna-path.exe",
        dynatype=dynatype,
        num_cpus=2,
        platform=platform,
        dyna_options="memory=1000",
        mpi_options="-hostfile myhostfile",
    )

    commands = settings.get_commands("path-to-input.k")

    if platform == "wsl":
        expected = ["powershell", "-Command", "wsl", "-e", "bash", "-lic", "./run_lsdyna.sh"]
    else:
        if dynatype == "smp":
            expected = ["my-dyna-path.exe", "i=path-to-input.k", "ncpu=2", "memory=1000"]
        elif dynatype in ["intelmpi", "platformmpi"]:
            expected = [
                "mpirun",
                "-hostfile myhostfile",
                "-np",
                "2",
                "my-dyna-path.exe",
                "i=path-to-input.k",
                "memory=1000",
            ]
        elif dynatype == "msmpi":
            expected = [
                "mpiexec",
                "-hostfile myhostfile",
                "-np",
                "2",
                "my-dyna-path.exe",
                "i=path-to-input.k",
                "memory=1000",
            ]

    assert commands == expected


@pytest.mark.parametrize(
    "dynatype",
    [
        "smp",
        "intelmpi",
        "platformmpi",
        "msmpi",
    ],
)
@pytest.mark.parametrize(
    "platform",
    [
        "windows",
        "linux",
        pytest.param(
            "wsl",
            marks=pytest.mark.xfail(not is_windows, reason="WSL Only valid argument on Windows"),
        ),
    ],
)
def test_get_dyna_commands_004(dynatype, platform):
    """Test if get commands returns right command line arguments if env variable is used."""
    if dynatype == "msmpi" and platform != "windows":
        pytest.skip("MSMPI and %s are not compatible and does not make sense to test." % platform)

    # define some mock data to test expand vars
    os.environ["TMP"] = "/some/tmp/directory/"
    mpi_options = "-hostfile $TMP"

    # define mock data
    settings = DynaSettings(
        lsdyna_path="my-dyna-path.exe",
        dynatype=dynatype,
        num_cpus=2,
        platform=platform,
        dyna_options="memory=1000",
        mpi_options=mpi_options,
    )

    commands = settings.get_commands("path-to-input.k")

    if platform == "wsl":
        expected = ["powershell", "-Command", "wsl", "-e", "bash", "-lic", "./run_lsdyna.sh"]
    else:
        if dynatype == "smp":
            expected = ["my-dyna-path.exe", "i=path-to-input.k", "ncpu=2", "memory=1000"]
        elif dynatype in ["intelmpi", "platformmpi"]:
            expected = [
                "mpirun",
                "-hostfile /some/tmp/directory/",
                "-np",
                "2",
                "my-dyna-path.exe",
                "i=path-to-input.k",
                "memory=1000",
            ]
        elif dynatype == "msmpi":
            expected = [
                "mpiexec",
                "-hostfile /some/tmp/directory/",
                "-np",
                "2",
                "my-dyna-path.exe",
                "i=path-to-input.k",
                "memory=1000",
            ]

    assert commands == expected


def test_modify_settings_from_env_variables(monkeypatch):
    """Test to ensure proper override from environment variables."""
    settings = DynaSettings(lsdyna_path="my-dyna-path.exe", dynatype="intelmpi", num_cpus=2)

    assert settings.lsdyna_path == "my-dyna-path.exe"
    assert settings.dynatype == "intelmpi"
    assert settings.num_cpus == 2

    # temporarily set the environment variables
    monkeypatch.setenv("PYANSYS_HEART_LSDYNA_PATH", "new-dyna-path.exe")
    monkeypatch.setenv("PYANSYS_HEART_LSDYNA_PLATFORM", "wsl")
    monkeypatch.setenv("PYANSYS_HEART_LSDYNA_TYPE", "msmpi")
    monkeypatch.setenv("PYANSYS_HEART_NUM_CPU", "4")

    settings = DynaSettings()

    assert settings.lsdyna_path == "new-dyna-path.exe"
    assert settings.platform == "wsl"
    assert settings.dynatype == "msmpi"
    assert settings.num_cpus == 4
