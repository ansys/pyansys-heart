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

"""Test to verify the computation of fiber inclination angles."""

import numpy as np

from ansys.health.heart.utils.misc import signed_angle_between_vectors


def test_compute_signed_angles():
    """Test computing signed angles between vectors."""

    # Local coordinate system
    e_c = np.array([1.0, 0.0, 0.0])  # circumferential direction
    e_t = np.array([0.0, 1.0, 0.0])  # transverse direction
    e_l = np.array([0.0, 0.0, 1.0])  # longitudinal direction

    # Test fiber directions
    fiber = np.array([e_c, -e_c, e_l, 2 * e_l, e_t, -e_t, [1.0, 1.0, 1.0], [1.0, 0.0, -1.0]])
    expected_angles = np.array(
        [
            0.0,  # angle between two equal vectors is 0
            180.0,  # angle between opposite vectors is 180
            90.0,  # angle between circumferential and longitudinal is 90 degrees
            90.0,  # angle between circumferential and longitudinal is 90 degrees
            np.nan,  # undefined due to parallel with normal vector
            np.nan,  # undefined due to parallel with normal vector
            45.0,  # angle should be 45 degrees for the vector [1, 1, 1]
            -45.0,  # angle should be -45 degrees for the vector [1, 0, -1]
        ]
    )

    # Extend input vectors to right shape
    num_fibers = fiber.shape[0]
    e_c = np.tile(e_c, (num_fibers, 1))
    e_t = np.tile(e_t, (num_fibers, 1))
    e_l = np.tile(e_l, (num_fibers, 1))

    computed_angles = signed_angle_between_vectors(x=fiber, y=e_c, n=e_t)

    assert np.allclose(computed_angles, expected_angles, equal_nan=True)
