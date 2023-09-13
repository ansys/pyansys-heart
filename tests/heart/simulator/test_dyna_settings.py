from ansys.heart.simulator.simulator import DynaSettings
import pytest


@pytest.mark.parametrize(
    "dynatype",
    [
        "smp",
        "intelmpi",
        "platformmpi",
        pytest.param("msmpi", marks=pytest.mark.xfail(reason="MSMPI not yet supported")),
    ],
)
@pytest.mark.parametrize("platform", ["windows", "linux", "wsl"])
def test_get_dyna_commands_001(dynatype, platform):
    """Test if get commands returns right command line if no additional options are given."""
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
        elif dynatype in ["intelmpi", "platformmpi", "msmpi"]:
            expected = ["mpirun", "-np", "2", "my-dyna-path.exe", "i=path-to-input.k"]

    assert commands == expected


@pytest.mark.parametrize(
    "dynatype",
    [
        "smp",
        "intelmpi",
        "platformmpi",
        pytest.param("msmpi", marks=pytest.mark.xfail(reason="MSMPI not yet supported")),
    ],
)
@pytest.mark.parametrize("platform", ["windows", "linux", "wsl"])
def test_get_dyna_commands_002(dynatype, platform):
    """Test if get commands returns right command line arguments if dyna options are given"""
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
        elif dynatype in ["intelmpi", "platformmpi", "msmpi"]:
            expected = [
                "mpirun",
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
        pytest.param("msmpi", marks=pytest.mark.xfail(reason="MSMPI not yet supported")),
    ],
)
@pytest.mark.parametrize("platform", ["windows", "linux", "wsl"])
def test_get_dyna_commands_003(dynatype, platform):
    """Test if get commands returns right command line arguments if dyna and mpi options present."""
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
        elif dynatype in ["intelmpi", "platformmpi", "msmpi"]:
            expected = [
                "mpirun",
                "-hostfile myhostfile",
                "-np",
                "2",
                "my-dyna-path.exe",
                "i=path-to-input.k",
                "memory=1000",
            ]

    assert commands == expected


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
