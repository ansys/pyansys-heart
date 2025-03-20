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

"""Collection of methods to test pyfluent."""

import pytest

import ansys.fluent.core as pyfluent
from ansys.fluent.core.launcher.error_handler import LaunchFluentError

# marks all tests with the 'requires_fluent' tag after this line
pytestmark = pytest.mark.requires_fluent

pyfluent.set_console_logging_level("DEBUG")


def test_launch_fluent():
    """Launch pyfluent in meshing mode and check health."""
    try:
        session = pyfluent.launch_fluent(
            mode="meshing",
            precision="double",
            processor_count=1,
            start_transcript=True,
            ui_mode="no_gui",
        )
        assert session._fluent_connection.check_health() == "SERVING"
        # try to initialize workflow
        assert session.workflow.InitializeWorkflow(WorkflowType="Watertight Geometry"), (
            "Failed workflow"
        )
        session.exit()
        assert True
    except (Exception, LaunchFluentError) as e:
        print(f"exception: {e}")
        assert False, "Failed to launch pyfluent in meshing mode."
