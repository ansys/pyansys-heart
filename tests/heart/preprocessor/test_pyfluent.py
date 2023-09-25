"""Collection of methods to test pyfluent."""
import os

import ansys.fluent.core as pyfluent
import pytest

# marks all tests with the 'requires_fluent' tag after this line
pytestmark = pytest.mark.requires_fluent

try:
    os.environ["GITHUB_JOB"]
    is_github_job = True
except KeyError:
    is_github_job = False

# disable Fluent gui for github job
if is_github_job:
    os.environ["SHOW_FLUENT_GUI"] = "0"


def test_launch_fluent():
    """Launch pyfluent in meshing mode and check health."""
    try:
        session = pyfluent.launch_fluent(
            mode="meshing",
            precision="double",
            processor_count=1,
            start_transcript=False,
            show_gui=False,
            product_version="23.1.0",
        )
        assert session.health_check_service.status() == "SERVING"
        # try to initialize workflow
        assert session.workflow.InitializeWorkflow(
            WorkflowType="Watertight Geometry"
        ), "Failed workflow"
        session.exit()
        assert True
    except:
        assert False, "Failed to launch pyfluent in meshing mode."
