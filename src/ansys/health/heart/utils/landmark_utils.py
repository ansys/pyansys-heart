# Copyright (C) 2023 - 2026 ANSYS, Inc. and/or its affiliates.
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

"""Module for computing heart anatomical landmarks.

This module provides utilities for computing heart anatomical landmarks including:

- Long and short axes of the left ventricle
- AHA17 segment landmarks and segmentation
- Coordinate systems for strain calculations
"""

from typing import Literal

from deprecated import deprecated
import numpy as np
import pyvista as pv
from scipy.spatial.transform import Rotation

from ansys.health.heart import LOG as LOGGER
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
    # Compute 4-chamber view (4CV) axis from cardiac valves and apex
    # This is the long axis passing through mitral and aortic valves
    center = np.mean(np.array([av_center, mv_center, apex]), axis=0)
    normal = np.cross(av_center - apex, mv_center - apex)
    l4cv_axis = {"center": center, "normal": normal / np.linalg.norm(normal)}

    # Compute short axis perpendicular to the long axis
    # The axis runs from the mitral valve center toward the apex
    sh_axis = apex - mv_center
    # Place the center at a fraction of the distance to avoid cutting through valve plane
    center = mv_center + first_cut_short_axis * sh_axis
    short_axis = {"center": center, "normal": sh_axis / np.linalg.norm(sh_axis)}

    # Compute 2-chamber view (2CV) axis perpendicular to both 4CV and short axes
    # This axis passes through mitral valve center and apex
    center = np.mean(np.array([mv_center, apex]), axis=0)
    p1 = center + 10 * l4cv_axis["normal"]
    p2 = mv_center
    p3 = apex
    normal = np.cross(p1 - p2, p1 - p3)
    l2cv_axis = {"center": center, "normal": normal / np.linalg.norm(normal)}

    return (l4cv_axis, l2cv_axis, short_axis)


def compute_aha17_points(
    model: HeartModel, short_axis: dict, long_axis: dict, surface: pv.PolyData = None
) -> list[np.ndarray]:
    """Compute 26 AHA17 landmark points on the left ventricle surface.

    The AHA17 landmark points define the geometry for the 17-segment model of the left ventricle.
    These 26 points are distributed across three levels (basal, mid, apical) and are located
    at specific angular positions around the circumference.

    Segment layout diagram::

        P1---H1---P2---H2---P3----H3--P4---H4---P5---H5---P6---H6---P1
        |          |         |         |         |         |         |
        |          |         |         |         |         |         |
        V1   S4   V2    S5  V3    S6  V4    S1  V5   S2   V6   S3   V1
        |          |         |         |         |         |         |
        |          |         |         |         |         |         |
        P7---H7---P8---H8---P9----H9--P10--H10--P11--H11--P12--H12--P7
        |          |         |         |         |         |         |
        |          |         |         |         |         |         |
        V7  S10   V8   S11  V9   S12  V10  S7   V11  S8   V12  S9   V7
        |          |         |         |         |         |         |
        |          |         |         |         |         |         |
        P13--H13--P14--H14--P15--H15--P16--H16--P17--H17--P18--H18--P13
              P19---H19---P20---H20---P21---H21---P22---H22---P19
               |           |           |           |           |
               |           |           |           |           |
               V13  S15   V14   S16   V15   S13   V16   S14   V13
               |           |           |           |           |
               |           |           |           |           |
              P23---H23---P24---H24---P25---H25---P26---H26---P23

    Parameters
    ----------
    model : HeartModel
        Heart model containing the left ventricle geometry.
    short_axis : dict
        Short axis definition with 'center' and 'normal' keys.
    long_axis : dict
        Long axis (4CV) definition with 'center' and 'normal' keys.
    surface : pv.PolyData, optional
        Endocardial surface on which to compute landmarks. If None, the endocardium from
        the heart model is used.

    Returns
    -------
    list[np.ndarray]
        List of 26 landmark point coordinates on the endocardial surface.
    """
    # Note: Epicardium is not supported because it may be incomplete for certain model types
    if surface is None:
        surface: pv.PolyData = model.left_ventricle.endocardium

    # Calculate longitudinal levels (basal, mid, apical) and apex points along the short axis
    p_basal, p_mid, p_apical, apex_endo, apex_epi = _calculate_longitudinal_points(
        model, short_axis
    )
    # Calculate rotation axes for circumferential segmentation boundaries
    axe_60, axe_120, axe_180, axe_45, axe_135 = _calculate_rotation_axis(short_axis, long_axis)

    # Compute landmark points for each level
    points = []

    # Basal and mid-cavity: 6 points each at different angular positions (60° spacing)
    for center in [p_basal, p_mid, p_apical]:
        for normal in [-axe_60, -axe_120, axe_180, axe_60, axe_120, -axe_180]:
            # Ray trace from the level center in each angular direction to find surface intersection
            point = surface.copy().ray_trace(center, center + 1e3 * normal, first_point=True)[0]
            if point.size == 0:
                raise ValueError("Cannot find point on surface from basal, mid, or apical.")
            else:
                points.append(point)

    # Apical region: 4 points at each level (45° and 135° angles)
    for center in [p_apical, apex_endo]:
        for normal in [-axe_45, -axe_135, axe_45, axe_135]:
            # Ray trace to find apical landmark points
            point = surface.copy().ray_trace(center, center + 1e3 * normal, first_point=True)[0]
            if point.size == 0:
                raise ValueError("Cannot find point on surface from apical or apex.")
            else:
                points.append(point)

    return points


def compute_aha17_lines(
    surface: pv.PolyData, points: list[np.ndarray]
) -> tuple[list[pv.PolyData], list[pv.PolyData]]:
    """Compute AHA 17 segment horizontal and vertical lines.

    Parameters
    ----------
    surface : pv.PolyData
        The endocardial surface to project the lines onto.
    points : list[np.ndarray]
        The list of AHA17 landmark points defining the lines.

    Returns
    -------
    tuple[list[pv.PolyData], list[pv.PolyData]]
        Tuple of horizontal and vertical line segments projected onto the surface.
    """
    # Calculate segment centers for line projection (mean of points in each group)
    p_basal = np.mean(np.array(points[0:6]), axis=0)
    p_mid = np.mean(np.array(points[6:12]), axis=0)
    p_apical = np.mean(np.array(points[12:22]), axis=0)
    apex_endo = np.mean(np.array(points[22:25]), axis=0)

    # Create horizontal (circumferential) lines at each level
    # Basal level horizontal lines
    hline = [
        _project_line_segment(surface, p_basal, points[0], points[1]),
        _project_line_segment(surface, p_basal, points[1], points[2]),
        _project_line_segment(surface, p_basal, points[2], points[3]),
        _project_line_segment(surface, p_basal, points[3], points[4]),
        _project_line_segment(surface, p_basal, points[4], points[5]),
        _project_line_segment(surface, p_basal, points[5], points[0]),
        _project_line_segment(surface, p_mid, points[6], points[7]),
        _project_line_segment(surface, p_mid, points[7], points[8]),
        _project_line_segment(surface, p_mid, points[8], points[9]),
        _project_line_segment(surface, p_mid, points[9], points[10]),
        _project_line_segment(surface, p_mid, points[10], points[11]),
        _project_line_segment(surface, p_mid, points[11], points[6]),
        _project_line_segment(surface, p_apical, points[12], points[13]),
        _project_line_segment(surface, p_apical, points[13], points[14]),
        _project_line_segment(surface, p_apical, points[14], points[15]),
        _project_line_segment(surface, p_apical, points[15], points[16]),
        _project_line_segment(surface, p_apical, points[16], points[17]),
        _project_line_segment(surface, p_apical, points[17], points[12]),
        _project_line_segment(surface, p_apical, points[18], points[19]),
        _project_line_segment(surface, p_apical, points[19], points[20]),
        _project_line_segment(surface, p_apical, points[20], points[21]),
        _project_line_segment(surface, p_apical, points[21], points[18]),
        _project_line_segment(surface, apex_endo, points[22], points[23]),
        _project_line_segment(surface, apex_endo, points[23], points[24]),
        _project_line_segment(surface, apex_endo, points[24], points[25]),
        _project_line_segment(surface, apex_endo, points[25], points[22]),
    ]

    # Create vertical (longitudinal) lines connecting levels
    # These lines span from basal through mid to apical levels
    vline = [
        _project_line_segment(surface, p_basal, points[0], points[6]),
        _project_line_segment(surface, p_basal, points[1], points[7]),
        _project_line_segment(surface, p_basal, points[2], points[8]),
        _project_line_segment(surface, p_basal, points[3], points[9]),
        _project_line_segment(surface, p_basal, points[4], points[10]),
        _project_line_segment(surface, p_basal, points[5], points[11]),
        _project_line_segment(surface, p_mid, points[6], points[12]),
        _project_line_segment(surface, p_mid, points[7], points[13]),
        _project_line_segment(surface, p_mid, points[8], points[14]),
        _project_line_segment(surface, p_mid, points[9], points[15]),
        _project_line_segment(surface, p_mid, points[10], points[16]),
        _project_line_segment(surface, p_mid, points[11], points[17]),
        _project_line_segment(surface, p_apical, points[18], points[22]),
        _project_line_segment(surface, p_apical, points[19], points[23]),
        _project_line_segment(surface, p_apical, points[20], points[24]),
        _project_line_segment(surface, p_apical, points[21], points[25]),
    ]

    for id, line in enumerate(hline):
        line.cell_data["id"] = id + 1

    for id, line in enumerate(vline):
        line.cell_data["id"] = id + 1

    return hline, vline


def _project_line_segment(
    surf: pv.PolyData, center: np.ndarray, p1: np.ndarray, p2: np.ndarray
) -> pv.PolyData:
    """Project a line segment onto a surface using ray tracing.

    Parameters
    ----------
    surf : pv.PolyData
        The surface onto which the line segment is projected.
    center : np.ndarray
        The center point used as the ray origin for tracing.
    p1 : np.ndarray
        The starting point of the line segment.
    p2 : np.ndarray
        The ending point of the line segment.

    Returns
    -------
    pv.PolyData
        A spline curve representing the projected line segment on the surface.
    """
    segment_points = np.linspace(p1, p2, num=10)

    # Project each point onto the surface using ray tracing from the center point
    projected_points = []

    for pt in segment_points:
        # Ray trace from center point through each segment point to find surface intersection
        start = center
        end = center + 1.0e3 * (pt - center)
        intersection = surf.copy().ray_trace(start, end, first_point=True)[0]
        if intersection.size == 0:
            # pv.PolyData([p1, p2, center,pt]).save("debug_points.vtp")
            # pv.Line(start,end).save("debug_points2.vtp")
            # surf.save("debug_endo.vtp")
            # raise ValueError(f"Cannot find intersection for segment.")
            LOGGER.warning("Cannot find intersection on the surface.")
        else:
            projected_points.append(intersection)

    projected_points = np.array(projected_points)

    # Create a smooth curve through the projected points
    projected_line = pv.Spline(projected_points, n_points=len(projected_points))

    return projected_line


def _calculate_longitudinal_points(model, short_axis):
    """Calculate basal, mid, and apical landmarks along the short axis.

    Parameters
    ----------
    model : HeartModel
        Heart model containing the left ventricle geometry.
    short_axis : dict
        Short axis information with 'center' and 'normal' keys defining the axis orientation.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        Basal center, mid-cavity center, apical center, endocardial apex, and epicardial apex.
    """
    # Extract endocardial and epicardial apex coordinates from the model
    for apex in model.left_ventricle.apex_points:
        if "endocardium" in apex.name:
            surface = model.left_ventricle.endocardium
            apex_ed = surface.points[np.where(surface["_global-point-ids"] == apex.node_id)[0][0]]
        elif "epicardium" in apex.name:
            surface = model.left_ventricle.epicardium
            apex_ep = surface.points[np.where(surface["_global-point-ids"] == apex.node_id)[0][0]]

    # Define three longitudinal levels along the short axis from basal to apical
    p_basal = short_axis["center"]
    p_mid = 1 / 3 * (apex_ep - p_basal) + p_basal
    p_apical = 2 / 3 * (apex_ep - p_basal) + p_basal

    # Project endocardial apex onto the short axis to create a flat apical segment (segment 17)
    # This ensures the endocardial apex point lies on the short axis plane
    x = apex_ed - apex_ep
    y = p_basal - apex_ep
    apex_ed = 1.5 * y * np.dot(x, y) / np.dot(y, y) + apex_ep
    return p_basal, p_mid, p_apical, apex_ed, apex_ep


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

    if not isinstance(model, HeartModel):
        raise TypeError("Model must be an instance of HeartModel.")

    # get elements of the left ventricle
    try:
        ele_ids = np.hstack(
            [
                model.mesh._global_tetrahedron_ids[
                    np.isin(
                        model.mesh._global_tetrahedron_ids,
                        model.left_ventricle.get_element_ids(model.mesh),
                    )
                ],
                model.mesh._global_tetrahedron_ids[
                    np.isin(
                        model.mesh._global_tetrahedron_ids,
                        model.septum.get_element_ids(model.mesh),
                    )
                ],
            ]
        )
        tetrahedrons = np.vstack(
            [
                model.left_ventricle._get_tetrahedrons(model.mesh),
                model.septum._get_tetrahedrons(model.mesh),
            ]
        )
    except AttributeError:
        tetrahedrons = model.left_ventricle._get_tetrahedrons(model.mesh)
        ele_ids = model.mesh._global_tetrahedron_ids[
            np.isin(
                model.mesh._global_tetrahedron_ids,
                model.left_ventricle.get_element_ids(model.mesh),
            )
        ]

    # Extract element center coordinates for segmentation
    elem_center = np.mean(model.mesh.points[tetrahedrons], axis=1)

    # Extract key anatomical landmarks and planes used for segmentation
    # Mitral valve center
    for cap in model.left_ventricle.caps:
        if cap.type == CapType.MITRAL_VALVE:
            mv_center = cap.centroid
    # Apical points on endocardium and epicardium
    for apex in model.left_ventricle.apex_points:
        if "endocardium" in apex.name:
            apex_ed = apex.xyz
        elif "epicardium" in apex.name:
            apex_ep = apex.xyz

    # Define reference cutting planes using short and long axes
    short_normal = short_axis["normal"]
    p_highest = short_axis["center"]

    # Create segmentation reference plane by rotating the long axis around the short axis
    # This defines the starting point for the 6 basal/mid segments and 4 apical segments
    if p_junction is not None:
        # LV-RV junction points provided, use them to define segment 1 boundary (CASIS definition)
        vec = (p_junction - p_highest) / np.linalg.norm(p_junction - p_highest)
        axe_60 = Rotation.from_rotvec(np.radians(90) * short_normal).apply(vec)
    else:
        # Default: rotate 60 degrees from the 4CV (long) axis
        long_axis = l4cv_axis["normal"]
        axe_60 = Rotation.from_rotvec(np.radians(60) * short_normal).apply(  # noqa:E501
            long_axis
        )

    axe_120 = Rotation.from_rotvec(np.radians(60) * short_normal).apply(axe_60)
    axe_180 = -Rotation.from_rotvec(np.radians(60) * short_normal).apply(axe_120)
    axe_45 = Rotation.from_rotvec(np.radians(-15) * short_normal).apply(axe_60)
    axe_135 = Rotation.from_rotvec(np.radians(90) * short_normal).apply(axe_45)

    # Calculate three longitudinal levels for segmentation (1/3 and 2/3 along the apex direction)
    p1_3 = 1 / 3 * (apex_ep - p_highest) + p_highest
    p2_3 = 2 / 3 * (apex_ep - p_highest) + p_highest

    # Project endocardial apex onto the short axis plane for a flat apical cap (segment 17)
    # This ensures segment 17 is perpendicular to the long axis
    x = apex_ed - apex_ep
    y = p_highest - apex_ep
    apex_ed = y * np.dot(x, y) / np.dot(y, y) + apex_ep

    # Assign AHA17 segment labels to each element based on its position
    label = np.full(len(elem_center), np.nan)
    for i, n in enumerate(elem_center):
        # Skip elements basal to the reference plane (typically containing valves)
        if np.dot(n - p_highest, mv_center - p_highest) > 0:
            continue
        # Basal level: segments 1, 2, 3, 4, 5, 6
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
        # Mid cavity level: segments 7, 8, 9, 10, 11, 12
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
        # Apical level: segments 13, 14, 15, 16, and optionally 17
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


def _calculate_rotation_axis(short: dict, long: dict):
    """Calculate rotation axes for AHA segment boundaries.

    Parameters
    ----------
    short : dict
        Short axis information with 'normal' key defining the rotation axis.
    long : dict
        Long axis (4CV) information with 'normal' key defining the reference direction.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        Five rotation axes at 60, 120, 180, 45, and 135 degree angles used for segmentation.
    """
    short_normal = short["normal"]

    # Calculate reference axes by rotating the long axis around the short axis normal
    # These define the boundaries between the 6 basal/mid and 4 apical segments
    long_axis = long["normal"]
    axe_60 = Rotation.from_rotvec(np.radians(60) * short_normal).apply(  # noqa:E501
        long_axis
    )

    axe_120 = Rotation.from_rotvec(np.radians(60) * short_normal).apply(axe_60)
    axe_180 = -Rotation.from_rotvec(np.radians(60) * short_normal).apply(axe_120)
    axe_45 = Rotation.from_rotvec(np.radians(-15) * short_normal).apply(axe_60)
    axe_135 = Rotation.from_rotvec(np.radians(90) * short_normal).apply(axe_45)

    return axe_60, axe_120, axe_180, axe_45, axe_135


@deprecated(reason="Use gradient from UVC to get better results.")
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
