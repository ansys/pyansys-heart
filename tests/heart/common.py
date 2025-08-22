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

"""Some common functions to test model stats."""


def compare_stats_names(stats: dict, stats_ref: dict):
    """Compare stats names.

    Parameters
    ----------
    stats : dict
        Dictionary with generated stats.
    stats_ref : dict
        Dictionary with reference stats.
    """
    for part_name in stats_ref["PARTS"].keys():
        assert part_name in list(stats["PARTS"].keys()), f"Part: {part_name} missing"

        for surface_name in stats_ref["PARTS"][part_name]["SURFACES"].keys():
            assert surface_name in list(stats["PARTS"][part_name]["SURFACES"].keys()), (
                f"Surface: {surface_name} missing"
            )

        for cap_name in stats_ref["PARTS"][part_name]["CAPS"].keys():
            assert cap_name in list(stats["PARTS"][part_name]["CAPS"].keys()), (
                f"Cap: {cap_name} missing"
            )

    assert sorted(list(stats["CAVITIES"].keys())) == sorted(list(stats_ref["CAVITIES"].keys())), (
        "one or more cavities missing"
    )

    return


def compare_stats_volumes(stats: dict, stats_ref: dict):
    """Compare stats volumes of cavities.

    Parameters
    ----------
    stats : dict
        Dictionary with generated stats.
    stats_ref : dict
        Dictionary with reference stats.
    """
    for cavity_name in stats_ref["CAVITIES"].keys():
        volume = stats["CAVITIES"][cavity_name]["volume"]
        volume_ref = stats_ref["CAVITIES"][cavity_name]["volume"]
        percent_diff = abs(volume - volume_ref) / volume_ref * 100
        assert percent_diff < 1, (
            f"Difference in cavity volumes of {cavity_name} is {percent_diff} percent"
        )

    return


def compare_stats_mesh(stats: dict, stats_ref: dict):
    """Compare mesh stats of the generated model."""
    assert stats["GENERAL"]["total_num_tets"] == stats_ref["GENERAL"]["total_num_tets"], (
        "Total number of tets not the same"
    )
    assert stats["GENERAL"]["total_num_nodes"] == stats_ref["GENERAL"]["total_num_nodes"], (
        "Total number of nodes not the same"
    )

    for part in stats_ref["PARTS"].keys():
        assert stats["PARTS"][part]["num_tets"] == stats_ref["PARTS"][part]["num_tets"], (
            f"Num tets in {part} not the same."
        )

        for surface in stats_ref["PARTS"][part]["SURFACES"].keys():
            ref_num_faces = stats_ref["PARTS"][part]["SURFACES"][surface]["num_faces"]
            num_faces = stats["PARTS"][part]["SURFACES"][surface]["num_faces"]
            assert num_faces == ref_num_faces, f"Number of faces of {surface} do not match."

        for cap in stats_ref["PARTS"][part]["CAPS"].keys():
            ref_num_nodes = stats_ref["PARTS"][part]["CAPS"][cap]["num_nodes"]
            num_nodes = stats["PARTS"][part]["CAPS"][cap]["num_nodes"]
            assert num_nodes == ref_num_nodes, f"Number of nodes of {cap} do not match."

    return
