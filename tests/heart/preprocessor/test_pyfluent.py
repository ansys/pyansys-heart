"""Collection of methods to test pyfluent."""
import sys

import pytest

if not sys.platform.startswith("win"):
    pytest.skip("Skipping windows-only tests", allow_module_level=True)

import ansys.fluent.core as pyfluent


def test_pyfluent():
    """Test launching pyfluent in meshing mode."""
    try:
        session = pyfluent.launch_fluent(
            meshing_mode=True,
            precision="double",
            processor_count=1,
            start_transcript=False,
            show_gui=False,
        )
        assert session.check_health() == "SERVING"
        session.exit()
        assert True
    except:
        assert False, "Failed to launch pyfluent in meshing model."
