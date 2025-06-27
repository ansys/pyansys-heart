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

"""Module for computing heart anatomical landmarks."""

from typing import Literal

from deprecated import deprecated
import numpy as np
from scipy.spatial.transform import Rotation

from ansys.health.heart.models import HeartModel
from ansys.health.heart.objects import CapType


def compute_anatomy_axis(
    mv_center: np.ndarray,
    av_center: np.ndarray,
    apex: np.ndarray,
    first_cut_short_axis: float = 0.2,
) -> tuple[dict, dict, dict]:
    """Compute the long and short axes of the left ventricle.

    Parameters
    ----------
    mv_center : np.ndarray
        Mitral valve center.
    av_center : np.ndarray
        Aortic valve center.
    apex : np.ndarray
        Left ventricle epicardium apex point.
    first_cut_short_axis : float, default: 0.2
        Relative distance between the mitral valve center and apex,
        which is used for defining the center of the short axis.

    Returns
    -------
    tuple[dict, dict, dict]
        4CV, 2CV, and short-axis. Each dictionary contains ``center`` and ``normal``.
    """
    # long 4CAV axis: cross apex, mitral and aortic valve centers
    center = np.mean(np.array([av_center, mv_center, apex]), axis=0)
    normal = np.cross(av_center - apex, mv_center - apex)
    l4cv_axis = {"center": center, "normal": normal / np.linalg.norm(normal)}

    # short axis: from mitral valve center to apex
    sh_axis = apex - mv_center
    # center is  highest possible point but avoid to cut aortic valve plane
    center = mv_center + first_cut_short_axis * sh_axis
    short_axis = {"center": center, "normal": sh_axis / np.linalg.norm(sh_axis)}

    # long 2CAV axis: normal to 4cav axe and pass mv center and apex
    center = np.mean(np.array([mv_center, apex]), axis=0)
    p1 = center + 10 * l4cv_axis["normal"]
    p2 = mv_center
    p3 = apex
    normal = np.cross(p1 - p2, p1 - p3)
    l2cv_axis = {"center": center, "normal": normal / np.linalg.norm(normal)}

    return (l4cv_axis, l2cv_axis, short_axis)


def compute_aha17(
    model: HeartModel,
    short_axis: dict,
    l4cv_axis: dict,
    seg: Literal[16, 17] = 17,
    p_junction: np.ndarray = None,
) -> np.ndarray:
    """Compute the AHA17 label for left ventricle elements.

    Parameters
    ----------
    model : HeartModel
        Heart model.
    short_axis : dict
        Short axis.
    l4cv_axis : dict
        Long 4CV axis.
    seg : Literal[16, 17], default: 17
        Compute 16 or 17 segments.
    p_junction : np.ndarray, default: None
        LV and RV junction points. If these points are given, they defines the start of segment 1.
        If they are not given, the start point is defined by rotating 60 degrees from the 4CV axis.

    Returns
    -------
    np.ndarray
        AHA17 IDs. No concerned elements are assigned with ``np.nan``.
    """
    aha_ids = np.full(len(model.mesh.tetrahedrons), np.nan)

    # get lv elements
    try:
        ele_ids = np.hstack(
            (
                model.left_ventricle.get_element_ids(model.mesh),
                model.septum.get_element_ids(model.mesh),
            )
        )
    except AttributeError:
        ele_ids = np.hstack(model.left_ventricle.get_element_ids(model.mesh))

    # element's center
    elem_center = np.mean(model.mesh.points[model.mesh.tetrahedrons[ele_ids]], axis=1)

    # anatomical points
    for cap in model.left_ventricle.caps:
        if cap.type == CapType.MITRAL_VALVE:
            mv_center = cap.centroid
    for apex in model.left_ventricle.apex_points:
        if "endocardium" in apex.name:
            apex_ed = apex.xyz
        elif "epicardium" in apex.name:
            apex_ep = apex.xyz

    # short axis
    short_normal = short_axis["normal"]
    p_highest = short_axis["center"]

    # define reference cut plane
    if p_junction is not None:
        # CASIS definition: LV and RV junction point
        vec = (p_junction - p_highest) / np.linalg.norm(p_junction - p_highest)
        axe_60 = Rotation.from_rotvec(np.radians(90) * short_normal).apply(vec)
    else:
        # default: rotate 60 from long axis
        long_axis = l4cv_axis["normal"]
        axe_60 = Rotation.from_rotvec(np.radians(60) * short_normal).apply(  # noqa:E501
            long_axis
        )

    axe_120 = Rotation.from_rotvec(np.radians(60) * short_normal).apply(axe_60)
    axe_180 = -Rotation.from_rotvec(np.radians(60) * short_normal).apply(axe_120)
    axe_45 = Rotation.from_rotvec(np.radians(-15) * short_normal).apply(axe_60)
    axe_135 = Rotation.from_rotvec(np.radians(90) * short_normal).apply(axe_45)

    p1_3 = 1 / 3 * (apex_ep - p_highest) + p_highest
    p2_3 = 2 / 3 * (apex_ep - p_highest) + p_highest

    # to have a flat segment 17, project endocardical apex point on short axis
    x = apex_ed - apex_ep
    y = p_highest - apex_ep
    apex_ed = y * np.dot(x, y) / np.dot(y, y) + apex_ep

    # aha17 label assignment
    label = np.full(len(elem_center), np.nan)
    for i, n in enumerate(elem_center):
        # This part contains valves, do not considered by AHA17
        if np.dot(n - p_highest, mv_center - p_highest) > 0:
            continue
        # Basal: segment 1 2 3 4 5 6
        elif np.dot(n - p1_3, mv_center - p1_3) >= 0:
            if np.dot(n - p1_3, axe_60) >= 0:
                if np.dot(n - p1_3, axe_120) >= 0:
                    if np.dot(n - p1_3, axe_180) >= 0:
                        label[i] = 5
                    else:
                        label[i] = 6
                else:
                    label[i] = 4
            else:
                if np.dot(n - p1_3, axe_180) <= 0:
                    if np.dot(n - p1_3, axe_120) >= 0:
                        label[i] = 1
                    else:
                        label[i] = 2
                else:
                    label[i] = 3
        # Mid cavity: segment 7 8 9 10 11 12
        elif np.dot(n - p2_3, mv_center - p2_3) >= 0:
            if np.dot(n - p1_3, axe_60) >= 0:
                if np.dot(n - p1_3, axe_120) >= 0:
                    if np.dot(n - p1_3, axe_180) >= 0:
                        label[i] = 11
                    else:
                        label[i] = 12
                else:
                    label[i] = 10
            else:
                if np.dot(n - p1_3, axe_180) <= 0:
                    if np.dot(n - p1_3, axe_120) >= 0:
                        label[i] = 7
                    else:
                        label[i] = 8
                else:
                    label[i] = 9
        # Apical
        else:
            if seg == 17:
                if np.dot(n - apex_ed, apex_ep - apex_ed) >= 0:
                    label[i] = 17
                else:
                    if np.dot(n - p1_3, axe_45) >= 0:
                        if np.dot(n - p1_3, axe_135) >= 0:
                            label[i] = 16
                        else:
                            label[i] = 15
                    else:
                        if np.dot(n - p1_3, axe_135) >= 0:
                            label[i] = 13
                        else:
                            label[i] = 14

            elif seg == 16:
                if np.dot(n - p1_3, axe_45) >= 0:
                    if np.dot(n - p1_3, axe_135) >= 0:
                        label[i] = 16
                    else:
                        label[i] = 15
                else:
                    if np.dot(n - p1_3, axe_135) >= 0:
                        label[i] = 13
                    else:
                        label[i] = 14

    aha_ids[ele_ids] = label
    return aha_ids


@deprecated(reason="Using gradient from UVC to get better results.")
def compute_element_cs(
    model: HeartModel, short_axis: dict, aha_element: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute elemental coordinate system for AHA elements.

    Parameters
    ----------
    model : HeartModel
        Heart model.
    short_axis : dict
        Short axis.
    aha_element : np.ndarray
        Elements with AHA labels. Compute only on these elements.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        Longitudinal, radial, and circumferential vectors of each AHA element.
    """
    elems = model.mesh.tetrahedrons[aha_element]
    elem_center = np.mean(model.mesh.points[elems], axis=1)

    # compute longitudinal direction, i.e. short axis
    e_l = np.tile(short_axis["normal"], (len(aha_element), 1))

    # compute radial direction
    center_offset = elem_center - model.left_ventricle.apex_points[1].xyz
    e_r = center_offset - (np.sum(e_l * center_offset, axis=1) * e_l.T).T
    # normalize each row
    e_r /= np.linalg.norm(e_r, axis=1)[:, np.newaxis]

    # compute circumferential direction
    e_c = np.cross(e_l, e_r)

    return e_l, e_r, e_c
