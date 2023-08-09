"""Collection of methods to test pyfluent."""

import ansys.fluent.core as pyfluent


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
        assert session.check_health() == "SERVING"
        # try to initialize workflow
        assert session.workflow.InitializeWorkflow(
            WorkflowType="Watertight Geometry"
        ), "Failed workflow"
        session.exit()
        assert True
    except:
        assert False, "Failed to launch pyfluent in meshing mode."
