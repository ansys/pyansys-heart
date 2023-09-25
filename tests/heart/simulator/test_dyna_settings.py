import os
import sys

from ansys.heart.simulator.settings.settings import DynaSettings
import pytest

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
    """Test if get commands returns right command line arguments if dyna options are given"""
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


@pytest.mark.xfail(is_gh_action, reason="No Ansys installation expected on Github runner.")
def test_set_env_variables():
    """Test setting environment variables."""

    import os

    if os.getenv("MPI_ROOT"):
        del os.environ["MPI_ROOT"]

    settings = DynaSettings(
        lsdyna_path="my-dyna-path.exe",
        dynatype="intelmpi",
        num_cpus=2,
        platform="windows",
    )
    settings._set_env_variables()

    assert os.getenv("MPI_ROOT")

    pass
