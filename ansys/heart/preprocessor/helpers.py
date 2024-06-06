# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
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

"""Module containing general helper methods."""

import os

from ansys.heart.core import LOG as LOGGER

heart_version = os.environ["ANSYS_HEART_MODEL_VERSION"]

if not heart_version or heart_version == "v0.1":
    from ansys.heart.preprocessor.models.v0_1.models import HeartModel
elif heart_version == "v0.2":
    from ansys.heart.preprocessor.models.v0_2.models import HeartModel


def model_summary(model: HeartModel, attributes: list = None) -> dict:
    """Generate a dictionary with model information.

    Parameters
    ----------
    model : HeartModel
        HeartModel for which to generate the summary dictionary
    attributes : list
        List of attributes to try to add to the dict.

    Returns
    -------
    dict
        Dictionary with model information.
    """
    sum_dict = {}
    sum_dict["GENERAL"] = {}

    try:
        sum_dict["GENERAL"]["total_num_tets"] = model.mesh.tetrahedrons.shape[0]
        sum_dict["GENERAL"]["total_num_nodes"] = model.mesh.nodes.shape[0]
    except TypeError:
        LOGGER.info("Failed to format General model information.")

    sum_dict["PARTS"] = {}
    sum_dict["CAVITIES"] = {}
    for ii, part in enumerate(model.parts):
        sum_dict["PARTS"][part.name] = {}
        sum_dict["PARTS"][part.name]["num_tets"] = len(part.element_ids)

        sum_dict["PARTS"][part.name]["SURFACES"] = {}
        sum_dict["PARTS"][part.name]["CAPS"] = {}

        for surface in part.surfaces:
            sum_dict["PARTS"][part.name]["SURFACES"][surface.name] = {}
            sum_dict["PARTS"][part.name]["SURFACES"][surface.name]["num_faces"] = (
                surface.triangles.shape[0]
            )

            if attributes:
                for attribute in attributes:
                    try:
                        sum_dict["PARTS"][part.name]["SURFACES"][surface.name][attribute] = getattr(
                            surface.clean(), attribute
                        )
                    except AttributeError:
                        pass

        for cap in part.caps:
            sum_dict["PARTS"][part.name]["CAPS"][cap.name] = {}
            sum_dict["PARTS"][part.name]["CAPS"][cap.name]["num_nodes"] = len(cap.node_ids)

            if attributes:
                for attribute in attributes:
                    try:
                        sum_dict["PARTS"][part.name]["CAPS"][cap.name][attribute] = getattr(
                            cap, attribute
                        )
                    except AttributeError:
                        pass

    for cavity in model.cavities:
        sum_dict["CAVITIES"][cavity.name] = {}
        sum_dict["CAVITIES"][cavity.name]["volume"] = cavity.surface.volume

        if attributes:
            for attribute in attributes:
                try:
                    sum_dict["CAVITIES"][cavity.name][attribute] = getattr(cavity, attribute)
                except AttributeError:
                    pass

    return sum_dict
