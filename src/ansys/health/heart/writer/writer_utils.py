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

"""Module for useful methods to help format LS-DYNA keywords."""

import numpy as np
import pandas as pd

from ansys.dyna.core import Deck
from ansys.dyna.core.keywords import keywords


def create_node_keyword(nodes: np.ndarray, offset: int = 0) -> keywords.Node:
    """Create node keyword from a NumPy array of nodes.

    Parameters
    ----------
    nodes : np.ndarray
        NumPy array containing the node coordinates.

    Returns
    -------
    keywords.Node
        Formatted node keyword.
    """
    # create array with node ids
    nids = np.arange(0, nodes.shape[0], 1) + 1

    kw = keywords.Node()
    df = pd.DataFrame(data=np.vstack([nids, nodes.T]).T, columns=kw.nodes.columns[0:4])

    kw.nodes = df

    return kw


def add_nodes_to_kw(nodes: np.ndarray, node_kw: keywords.Node, offset: int = 0) -> keywords.Node:
    """Add nodes to an existing *NODE keyword.

    Notes
    -----
    If nodes are already defined, this method adds both the nodes in the previous
    keyword and the specified array of nodes. It automatically computes
    the index offset in case ``node_kw.nodes`` is not empty.

    Parameters
    ----------
    nodes : np.ndarray
        NumPy array of node coordinates to add.
        If (n,3), the node ID is continuous by offset.
        If (n,4), the first column is the node ID.
    node_kw : keywords.Node
        Node keyword.
    offset : int
        Node ID offset.
    """
    if nodes.shape[1] == 4:
        df = pd.DataFrame(data=nodes, columns=node_kw.nodes.columns[0:4])
    elif nodes.shape[1] == 3:
        # get node id of last node:
        if not node_kw.nodes.empty and offset == 0:
            last_nid = node_kw.nodes.iloc[-1, 0]
            offset = last_nid

        # create array with node IDs
        nids = np.arange(0, nodes.shape[0], 1) + offset + 1

        # create dataframe
        df = pd.DataFrame(data=np.vstack([nids, nodes.T]).T, columns=node_kw.nodes.columns[0:4])

    # concatenate old and new dataframe
    df1 = pd.concat([node_kw.nodes, df], axis=0, ignore_index=True, join="outer")

    node_kw = keywords.Node()
    node_kw.nodes = df1

    return node_kw


def add_beams_to_kw(
    beams: np.ndarray, beam_kw: keywords.ElementBeam, pid: int, offset: int = 0
) -> keywords.ElementBeam:
    """Add beams to an existing *ELEMENT_BEAM keyword.

    Notes
    -----
    If beams are already defined, this method adds both the beams in the previous
    keyword and the specified array of beams. It automatically computes
    the index offset in case ``beam_kw.elements`` is not empty.

    Parameters
    ----------
    beams : np.ndarray
        NumPy array of beam coordinates to add.
    beam_kw : keywords.ElementBeam
        Beam keyword.
    offset : int
        Beam ID offset.
    """
    # get beam id of last beam:
    if not beam_kw.elements.empty and offset == 0:
        last_eid = beam_kw.elements.iloc[-1, 1]
        offset = last_eid

    # create array with beam ids
    eids = np.arange(0, beams.shape[0], 1) + offset + 1
    pid = np.zeros(eids.shape) + pid
    # create dataframe
    df = pd.DataFrame(
        data=np.vstack([eids, pid, beams.T]).T,
        columns=beam_kw.elements.columns[0:4],
        dtype=int,
    )

    # concatenate old and new dataframe
    df1 = pd.concat([beam_kw.elements, df], axis=0, ignore_index=True, join="outer")

    beam_kw = keywords.ElementBeam()
    beam_kw.elements = df1

    return beam_kw


def create_segment_set_keyword(
    segments: np.ndarray, segid: int = 1, title: str = ""
) -> keywords.SetSegment:
    """Create a segment set keyword from an array with the segment set definition.

    Parameters
    ----------
    segments : np.ndarray
        Array of node indices that make up the segment. If three columns are provided,
        it is assumed that the segments are triangular
    segid : int, default: 1
        Segment set ID.
    title : str, default: ""
        Title of the segment set.

    Returns
    -------
    keywords.SetSegment
        Formatted segment set keyword.
    """
    if segments.shape[1] < 3 or segments.shape[1] > 4:
        raise ValueError("Expecting segments to have 3 or 4 columns.")

    if segments.shape[1] == 3:
        # segtype = "triangle"
        segments = np.vstack([segments.T, segments[:, -1]]).T

    kw = keywords.SetSegment(sid=segid)

    if title != "":
        kw.options["TITLE"].active = True
        kw.title = title

    # prepare dataframe
    df = pd.DataFrame(data=segments, columns=kw.segments.columns[0:4])

    kw.segments = df

    return kw


def create_node_set_keyword(
    node_ids: np.ndarray | list[int] | int,
    node_set_id: int = 1,
    title: str = "nodeset-title",
) -> keywords.SetNodeList:
    """Create a nodeset.

    Parameters
    ----------
    node_ids : np.ndarray | list[int] | int
        List of node IDs to include in the nodeset.
    node_set_id : int, default: 1
        ID of the nodeset.
    title : str, default: ``'nodeset-title'``
        Title of the nodeset

    Returns
    -------
    keywords.SetNodeList
        Formatted nodeset.
    """
    if not isinstance(node_ids, (np.ndarray, int, np.int32, np.int64, list)):
        raise TypeError("Expecting node IDs to be array or list of integers or a single integer.")
    if isinstance(node_ids, (int, np.int32, np.int64)):
        node_ids = [node_ids]

    kw = keywords.SetNodeList(sid=node_set_id)

    kw.options["TITLE"].active = True
    kw.title = title
    kw.nodes._data = node_ids
    kw._cards[0].set_value("its", None)  # remove its parameter.

    return kw


def create_element_shell_keyword(
    shells: np.ndarray, part_id: int = 1, id_offset: int = 0
) -> keywords.ElementShell:
    """Create an element shell keyword.

    Notes
    -----
    This method creates an element shell keyword from a NumPy array of elements.
    Each row corresponds to an element.

    """
    num_shells = shells.shape[0]

    kw = keywords.ElementShell()

    if shells.shape[1] == 3:
        # element_type = "triangle"
        columns = kw.elements.columns[0:5]
    elif shells.shape[1] == 4:
        # element_type = "quad"
        columns = kw.elements.columns[0:6]
    else:
        raise ValueError("Type is unknown. Check size of shell array.")

    # create element id array
    element_ids = np.arange(0, num_shells, 1) + 1 + id_offset
    part_ids = np.ones(element_ids.shape) * part_id
    data = np.vstack((element_ids, part_ids, shells.T)).T
    # create pandas dataframe
    df = pd.DataFrame(data=data, columns=columns)
    kw.elements = df
    return kw


def create_element_solid_keyword(
    elements: np.ndarray,
    e_id: np.ndarray,
    part_id: np.ndarray,
    element_type: str = "tetra",
) -> keywords.ElementSolid:
    """Format the *ELEMENT_SOLID keyword with the provided input.

    Parameters
    ----------
    elements : np.ndarray
        NumPy array of integers with element definition.
    part_id : np.ndarray
        Part IDs of each element.
    e_id : np.ndarray
        Element ID.
    element_type : str, default: ``'tetra'``
        Type of element to write

    Returns
    -------
    keywords.ElementSolid
        Formatted *ELEMENT_SOLID keyword.
    """
    kw = keywords.ElementSolid()
    df = pd.DataFrame(columns=kw.elements)

    # prepare element data for dataframe
    df["eid"] = e_id
    df["pid"] = part_id
    if element_type == "tetra":
        df["n1"] = elements[:, 0]
        df["n2"] = elements[:, 1]
        df["n3"] = elements[:, 2]
        df["n4"] = elements[:, 3]
        df["n5"] = elements[:, 3]
        df["n6"] = elements[:, 3]
        df["n7"] = elements[:, 3]
        df["n8"] = elements[:, 3]

        kw.elements = df

        return kw


def create_element_solid_ortho_keyword(
    elements: np.ndarray,
    a_vec: np.ndarray,
    d_vec: np.ndarray,
    e_id: np.ndarray,
    part_id: np.ndarray,
    element_type: str = "tetra",
) -> keywords.ElementSolidOrtho:
    """Format the *ELEMENT_SOLID_ORTHO keyword with the provided input.

    Parameters
    ----------
    elements : np.ndarray
        NumPy array of integers with element definition
    a_vec : np.ndarray
        Vector specifying the A direction.
    d_vec : np.ndarray
        Vector specifying the D direction.
    part_id : np.ndarray
        Part IDs of each element.
    e_id : np.ndarray
        Element ID.
    element_type : str, default: ``'tetra'``
        Type of element to write.

    Returns
    -------
    keywords.ElementSolidOrtho
        Formatted *ELEMENT_SOLID_ORTHO keyword.
    """
    kw = keywords.ElementSolidOrtho()

    df = pd.DataFrame(columns=kw.elements)

    # prepare element data for dataframe
    df["eid"] = e_id
    df["pid"] = part_id

    if element_type == "tetra":
        df["n1"] = elements[:, 0]
        df["n2"] = elements[:, 1]
        df["n3"] = elements[:, 2]
        df["n4"] = elements[:, 3]
        df["n5"] = elements[:, 3]
        df["n6"] = elements[:, 3]
        df["n7"] = elements[:, 3]
        df["n8"] = elements[:, 3]

    df["a1"] = a_vec[:, 0]
    df["a2"] = a_vec[:, 1]
    df["a3"] = a_vec[:, 2]

    df["d1"] = d_vec[:, 0]
    df["d2"] = d_vec[:, 1]
    df["d3"] = d_vec[:, 2]

    kw.elements = df

    return kw


def create_define_curve_kw(
    x: np.ndarray,
    y: np.ndarray,
    curve_name: str = "my-title",
    curve_id: int = 1,
    lcint: int = 15000,
) -> keywords.DefineCurve:
    """Create define curve from x and y values."""
    kw = keywords.DefineCurve()
    kw.options["TITLE"].active = True
    kw.title = curve_name
    kw.lcid = curve_id
    kw.lcint = lcint
    columns = kw.curves.columns
    data = np.array((x, y), dtype=float).T
    df = pd.DataFrame(data=data, columns=columns)
    kw.curves = df

    return kw


def create_define_sd_orientation_kw(
    vectors: np.ndarray, vector_id_offset: int = 0, iop: int = 0
) -> keywords.DefineSdOrientation:
    """Create define SD orientation keyword.

    Parameters
    ----------
    vectors : np.ndarray
        Array of shape Nx3 with the defined vector.
    vector_id_offset : int, default: 0
        Offset for the vector ID.
    iop : int, default: 0
        Option.
    """
    kw = keywords.DefineSdOrientation()
    if len(vectors.shape) == 2:
        num_vectors = vectors.shape[0]
    elif len(vectors.shape) == 1:
        num_vectors = 1
        vectors = vectors[None, :]

    vector_ids = np.arange(0, num_vectors, 1) + vector_id_offset + 1
    iops = np.ones(num_vectors) * iop
    columns = kw.vectors.columns
    data = np.hstack([vector_ids[:, None], iops[:, None], vectors], dtype=np.float32)
    df = pd.DataFrame(data=data, columns=columns[0:5])

    kw.vectors = df

    return kw


def create_discrete_elements_kw(
    nodes: np.ndarray,
    part_id: int,
    vector_ids: np.ndarray | int,
    scale_factor: np.ndarray | float,
    element_id_offset: int = 0,
) -> keywords.ElementDiscrete:
    """Create discrete elements based on input arguments.

    Parameters
    ----------
    nodes : np.ndarray
        Nx2 array with node IDs used for the discrete element.
    part_id : int
        Part ID of the discrete elements given.
    vector_ids : np.ndarray | int
        Orientation IDs (vector IDs) that the spring acts on.
        You can provide either an array of length N or a scalar integer.
    scale_factor : np.ndarray | float
        Scale factor on forces. You can provide either an array of length N
        or a scalar value.
    element_id_offset : int, default: 0
        Offset value for the element IDs.
    init_offset : float, default: 0.0
        Initial offset, which is the initial displacement or rotation at t=0.
    """
    num_elements = nodes.shape[0]

    if isinstance(vector_ids, int):
        vector_ids = np.ones((num_elements), dtype=float) * vector_ids

    if isinstance(scale_factor, float):
        scaling = np.ones((num_elements, 1), dtype=float) * scale_factor
    elif isinstance(scale_factor, np.ndarray):
        scaling = scale_factor[:, None]

    kw = keywords.ElementDiscrete()
    columns = kw.elements.columns

    element_ids = np.arange(0, num_elements, 1) + 1 + element_id_offset
    element_ids = element_ids[:, None]
    part_ids = np.ones((num_elements, 1), dtype=float) * part_id

    data = np.hstack([element_ids, part_ids, nodes, vector_ids[:, None], scaling])

    # set up dataframe
    df = pd.DataFrame(data=data, columns=columns[0:6])

    kw.elements = df

    return kw


def get_list_of_used_ids(keyword_db: Deck, keyword_str: str) -> np.ndarray:
    """Get array of used IDs in the database.

    Notes
    -----
    For example, for *SECTION, you would get *PART and *MAT IDs.

    Parameters
    ----------
    database : Deck
        Database of keywords.
    keyword : str
        Keyword to find.

    Returns
    -------
    np.ndarray
        Array of IDs (integers) that are already used
    """
    ids = np.empty(0, dtype=int)

    valid_kws = [
        "SECTION",
        "PART",
        "MAT",
        "SET_SEGMENT",
        "SET_NODE",
        "DEFINE_CURVE",
        "SET_PART",
    ]

    if keyword_str not in valid_kws:
        raise ValueError("Expecting one of: {0}".format(valid_kws))

    if keyword_str == valid_kws[0]:
        for kw in keyword_db.get_kwds_by_type(valid_kws[0]):
            ids = np.append(ids, kw.secid)

    if keyword_str == valid_kws[1]:
        for kw in keyword_db.get_kwds_by_type(valid_kws[1]):
            ids = np.append(ids, kw.parts["pid"].to_numpy())

    if keyword_str == valid_kws[2]:
        for kw in keyword_db.get_kwds_by_type(valid_kws[2]):
            ids = np.append(ids, kw.mid)

    # special treatment for segment sets:
    if keyword_str == valid_kws[3]:
        for kw in keyword_db.get_kwds_by_type("SET"):
            if "SEGMENT" in kw.subkeyword:
                ids = np.append(ids, kw.sid)

    # special treatment for nodesets
    if keyword_str == valid_kws[4]:
        for kw in keyword_db.get_kwds_by_type("SET"):
            if "NODE" in kw.subkeyword:
                ids = np.append(ids, kw.sid)

    # special treatment for nodesets
    if keyword_str == valid_kws[6]:
        for kw in keyword_db.get_kwds_by_type("SET"):
            if "PART" in kw.subkeyword:
                ids = np.append(ids, kw.sid)

    if keyword_str == valid_kws[5]:
        for kw in keyword_db.get_kwds_by_type("DEFINE"):
            if "CURVE" in kw.subkeyword:
                ids = np.append(ids, kw.lcid)
            elif "FUNCTION" in kw.subkeyword:
                ids = np.append(ids, kw.fid)

    return ids


def fast_element_writer(
    element_kw: keywords.ElementSolidOrtho | keywords.ElementSolid, filename: str
) -> None:
    """Fast implementation of the element writer.

    Notes
    -----
    Use this method as an alternative to the PyDYNA keywords module.

    """
    # TODO: generalize this writer

    if element_kw.subkeyword == "SOLID":
        writer = "solid_writer"
    elif element_kw.subkeyword == "SOLID_ORTHO":
        writer = "solid_ortho_writer"

    # remove columns n9 and n10 if they exist
    elements = element_kw.elements
    try:
        elements = elements.drop("n9", axis=1)
        elements = elements.drop("n10", axis=1)
    except KeyError:
        pass

    elements_np = elements.to_numpy()
    # headers = list(element_kw.elements.columns)

    if writer == "solid_ortho_writer":
        # explicitly cast to ints and floats
        elements = np.array(elements_np[:, 0:10], dtype=int)
        elements_advectors = np.array(elements_np[:, 10:], dtype=float)

        # create list of formatted strings
        list_formatted_strings = []
        line_format = (
            "{:8d}" * 2  # element ID and part ID
            + "\n"
            + "{:8d}" * 8  # node IDs
            + "\n"
            + "{:16e}" * 3  # fiber vector
            + "\n"
            + "{:16e}" * 3  # sheet vector
            + "\n"
        )
        for ii, element in enumerate(elements):
            line_format_str = line_format.format(
                element[0],
                element[1],
                element[2],
                element[3],
                element[4],
                element[5],
                element[6],
                element[7],
                element[8],
                element[9],
                elements_advectors[ii][0],
                elements_advectors[ii][1],
                elements_advectors[ii][2],
                elements_advectors[ii][3],
                elements_advectors[ii][4],
                elements_advectors[ii][5],
            )
            list_formatted_strings.append(line_format_str)

        fid = open(filename, "a")
        fid.write("*ELEMENT_SOLID_ORTHO\n")
        for line in list_formatted_strings:
            fid.write(line)
        fid.close()

    elif writer == "solid_writer":
        list_formatted_strings = []
        # explicitly cast to ints
        elements = np.array(elements_np[:, 0:10], dtype=int)
        line_format = (
            "{:8d}" * 2 + "\n" + "{:8d}" * 8 + "\n"
        )  # element ID and part ID. n1, n2, n3, ...
        for element in elements:
            line_format_str = line_format.format(
                element[0],
                element[1],
                element[2],
                element[3],
                element[4],
                element[5],
                element[6],
                element[7],
                element[8],
                element[9],
            )
            list_formatted_strings.append(line_format_str)
        fid = open(filename, "a")
        fid.write("*ELEMENT_SOLID\n")
        for line in list_formatted_strings:
            fid.write(line)
        fid.close()

    return
