"""Module that contains some useful methods to help format the keywords
"""

from multiprocessing.sharedctypes import Value
import pandas as pd
import numpy as np
from typing import List, Union

from ansys.dyna import keywords
from ansys.dyna.keywords import Deck

# import some custom keywords that avoid buygs in dynalib
from ansys.heart.writer.custom_dynalib_keywords._custom_set_node_list import (
    SetNodeList_custom,
)
from ansys.heart.writer.custom_dynalib_keywords._custom_set_segment_add import (
    SetSegmentAdd_custom,
)




def create_node_keyword(nodes: np.array, offset: int = 0) -> keywords.Node:
    """Creates node keyword from numpy array of nodes

    Parameters
    ----------
    nodes : np.array
        Numpy array containing the node coordinates

    Returns
    -------
    keywords.Node
        Formatted node keyword
    """
    # create array with node ids
    nids = np.arange(0, nodes.shape[0], 1) + 1

    kw = keywords.Node()
    df = pd.DataFrame(
        data=np.vstack([nids, nodes.T]).T, columns=kw.nodes.columns[0:4]
    )

    kw.nodes = df

    return kw


def add_nodes_to_kw(
    nodes: np.array, node_kw: keywords.Node, offset: int = 0
) -> keywords.Node:
    """Adds nodes to an existing node keyword. If nodes are
    already defined, this adds both the nodes in the previous
    keyword and the specified array of nodes. Automatically computes
    the index offset incase node_kw.nodes is not empty

    Parameters
    ----------
    nodes : np.array
        Numpy array of node coordinates to add
    node_kw : keywords.Node
        Node keyword
    offset : int
        Node id offset
    """

    # get node id of last node:
    if not node_kw.nodes.empty and offset == 0:
        last_nid = node_kw.nodes.iloc[-1, 0]
        offset = last_nid

    # create array with node ids
    nids = np.arange(0, nodes.shape[0], 1) + offset + 1

    # create dataframe
    df = pd.DataFrame(
        data=np.vstack([nids, nodes.T]).T, columns=node_kw.nodes.columns[0:4]
    )

    # concatenate old and new dataframe
    df1 = pd.concat(
        [node_kw.nodes, df], axis=0, ignore_index=True, join="outer"
    )

    node_kw = keywords.Node()
    node_kw.nodes = df1

    return node_kw


def create_segment_set_keyword(
    segments: np.array, segid: int = 1, title: str = "seg-title"
) -> keywords.SetSegment:
    """Creates a segment set keyword from an array with the segment set definition

    Parameters
    ----------
    segments : np.array
        Array of node-indices that make up the segment. If three columns are provided 
        it is assumed that the segments are triangular
    segid : int, optional
        Segment set ID, by default 1
    title : str, optional
        Title of the segment set, by default 'seg-title'

    Returns
    -------
    keywords.SetSegment
        Formatted segment set keyword
    """

    if segments.shape[1] < 3 or segments.shape[1] > 4:
        raise ValueError("expecting segments to have 3 or 4 columns")

    if segments.shape[1] == 3:
        segtype = "triangle"
        segments = np.vstack([segments.T, segments[:, -1]]).T

    kw = keywords.SetSegment(sid=segid)

    kw.options["TITLE"].active = True
    kw.title = title

    # prepare dataframe
    df = pd.DataFrame(data=segments, columns=kw.segments.columns[0:4])

    kw.segments = df

    return kw


def create_node_set_keyword(
    node_ids: np.array, node_set_id: int = 1, title: str = "nodeset-title"
) -> SetNodeList_custom:
    """Creates node set

    Parameters
    ----------
    node_ids : np.array
        List of node ids to include in the node set
    node_set_id : int, optional
        Id of the node set, by default 1
    title : str, optional
        Title of the node set, by default 'nodeset-title'

    Returns
    -------
    keywords.SetNodeList
        Formatted node set 
    """

    kw = SetNodeList_custom(sid=node_set_id)

    kw.options["TITLE"].active = True
    kw.title = title

    # rearrange nodeids such that its size is n x 8
    num_nodes = len(node_ids)
    num_cols = int(8)
    num_rows = int(np.ceil(len(node_ids) / 8))

    node_ids2 = np.zeros(int(num_rows * num_cols)) * np.nan

    node_ids2[:num_nodes] = node_ids

    node_ids2 = np.reshape(node_ids2, (num_rows, num_cols))

    # prepare dataframe
    df = pd.DataFrame(data=node_ids2, columns=kw.nodes.columns[0:8])

    kw.nodes = df

    return kw

def create_element_shell_keyword(
    shells: np.array, part_id: int = 1, id_offset: int = 0
) -> keywords.ElementShell:

    """Creates element shell keyword from a numpy array of 
    elements where each row corresponds to an element.
    """
    num_shells = shells.shape[0]

    kw = keywords.ElementShell()

    if shells.shape[1] == 3:
        element_type = "triangle"
        columns = kw.elements.columns[0:5]
    elif shells.shape[1] == 4:
        element_type = "quad"
        columns = kw.elements.columns[0:6]
    else:
        raise ValueError("Uknown type. Check size of shell array")

    # create element id array
    element_ids = np.arange(0, num_shells, 1) + 1 + id_offset
    part_ids = np.ones(element_ids.shape) * part_id
    data = np.vstack((element_ids, part_ids, shells.T)).T
    # create pandas dataframe
    df = pd.DataFrame(data=data, columns=columns)
    kw.elements = df
    return kw


def create_element_solid_ortho_keyword(
    elements: np.array,
    a_vec: np.array,
    d_vec: np.array,
    partid: int = 1,
    id_offset: int = 0,
    element_type="tetra",
) -> keywords.ElementSolidOrtho:
    """Formats the *ELEMENT_SOLID_ORTHO keyword with the provided input

    Parameters
    ----------
    elements : np.array
        Numpy array of ints with element definition
    a_vec : np.array
        Vector specifying the A direction
    d_vec : np.array
        Vector specifying the D direction
    partid : int, optional
        Part id to be used, by default 1
    id_offset : int, optional
        Offset of the element number, by default 0
    element_type : str, optional
        Type of element to write, by default "tetra"

    Returns
    -------
    keywords.ElementSolidOrtho
        Formatted *ELEMENT_SOLID_ORTHO keyword
    """

    kw = keywords.ElementSolidOrtho()

    df = pd.DataFrame(columns=kw.elements)

    # prepare element data for dataframe
    elids = np.arange(1, elements.shape[0] + 1, 1) + id_offset
    partids = np.ones(elements.shape[0]) * partid

    df["eid"] = elids
    df["pid"] = partids

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
    x: np.array,
    y: np.array,
    curve_name: str = "my-title",
    curve_id: int = 1,
    lcint: int = 15000,
) -> keywords.DefineCurve:
    """Creates define curve from x and y values"""
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
    vectors: np.array, vector_id_offset: int = 0, iop: int = 0
) -> keywords.DefineSdOrientation:
    """Creates define SD orientation keyword

    Parameters
    ----------
    vectors : np.array
        Array of shape Nx3 with the defined vector
    vector_id_offset : int, optional
        Offset for the vector id, by default 0
    iop : int, optional
        Option, by default 0
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
    data = np.hstack([vector_ids[:, None], iops[:, None], vectors])
    df = pd.DataFrame(data=data, columns=columns[0:5])

    kw.vectors = df

    return kw


def create_discrete_elements_kw(
    nodes: np.array,
    part_id: int,
    vector_ids: Union[np.array, int],
    scale_factor: Union[np.array, float],
    element_id_offset: int = 0,
) -> keywords.ElementDiscrete:
    """Creates discrete elements based on the input arguments

    Parameters
    ----------
    nodes : np.array
        Nx2 Array with node ids used for the discrete element
    part_id : int
        Part id of the discrete elements given
    vector_ids : Union[np.array, int]
        Orientation ids (vector ids) along which the spring acts. Can be either an array of length N, or a scalar integer
    scale_factor : Union[np.array, float]
        Scale factor on forces, either an array of length N or scalar value
    element_id_offset : int, optional
        Offset value for the element ids, by default 0
    init_offset : float, optional
        Initial offset: initial displacement or rotation at t=0, by default 0.0
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

    data = np.hstack(
        [element_ids, part_ids, nodes, vector_ids[:, None], scaling]
    )

    # set up dataframe
    df = pd.DataFrame(data=data, columns=columns[0:6])

    kw.elements = df

    return kw


def get_list_of_used_ids(keyword_db: Deck, keyword_str: str) -> np.array:
    """Gets array of used ids in the database. E.g. for *SECTION, *PART and *MAT ids

    Parameters
    ----------
    database : Deck
        Database of keywords
    keyword : str
        Keyword which to find

    Returns
    -------
    np.array
        Array of ids which are already used
    """
    ids = np.empty(0, dtype=int)

    valid_kws = ["SECTION", "PART", "MAT", "SET_SEGMENT", "SET_NODE"]

    if keyword_str not in valid_kws:
        raise ValueError("Expecting one of: {0}".format(valid_kws))

    if keyword_str == valid_kws[0]:
        for kw in keyword_db.get_kwds_by_type(valid_kws[0]):
            ids = np.append(ids, kw.secid)

    if keyword_str == valid_kws[1]:
        for kw in keyword_db.get_kwds_by_type(valid_kws[1]):
            ids = np.append(ids, kw.parts["pid"].to_numpy() )

    if keyword_str == valid_kws[2]:
        for kw in keyword_db.get_kwds_by_type(valid_kws[2]):
            ids = np.append(ids, kw.mid)

    # special treatment for segment sets:
    if keyword_str == valid_kws[3]:
        for kw in keyword_db.get_kwds_by_type("SET"):
            if "SEGMENT" in kw.subkeyword:
                ids = np.append(ids, kw.sid)

    # special treatment for node sets
    if keyword_str == valid_kws[4]:
        for kw in keyword_db.get_kwds_by_type("SET"):
            if "NODE" in kw.subkeyword:
                ids = np.append(ids, kw.sid)

    return ids


def fast_element_writer(element_kw: keywords.ElementSolidOrtho, filename: str):
    """Fast implementation of the element writer. Use this as an alternative to
    the dynalib writer
    """
    elements = element_kw.elements.to_numpy()
    headers = list(element_kw.elements.columns)

    # create list of formatted strings
    list_formatted_strings = []
    line_format = (
        "{:>8.0f}" * 10
        + "\n"
        + "{:>16.5e}" * 3
        + "\n"
        + "{:>16.5e}" * 3
        + "\n"
    )
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
            element[10],
            element[11],
            element[12],
            element[13],
            element[14],
            element[15],
        )
        list_formatted_strings.append(line_format_str)

    fid = open(filename, "a")
    fid.write("*ELEMENT_SOLID_ORTHO\n")
    for line in list_formatted_strings:
        fid.write(line)
    fid.close()

    # fid.writelines( list_formatted_strings )

    # for line in list_formatted_strings:

    return


def example_performance():

    import time as time
    import pandas as pd
    import numpy as np
    from ansys.dyna import keywords
    from ansys.dyna.keywords import db as db

    # create some data
    num_elem = 40000

    elements = np.zeros((num_elem, 10), dtype=int)
    elements[:, 0] = np.arange(1, num_elem + 1, 1)
    elements[:, 1:] = [1, 2, 3, 4, 5, 6, 7, 8, 9]

    A = np.zeros((num_elem, 3), dtype=float) + 0.33
    B = np.ones((num_elem, 3), dtype=float) + 0.5

    data = np.hstack([elements, A, B])

    data[1, 7:] = np.nan

    # instantiate keyword
    kw = keywords.ElementSolidOrtho()
    element_df = pd.DataFrame(data=data, columns=kw.elements.columns)
    kw.elements = element_df

    ## faster export method:
    t0 = time.time()
    line_format = ""
    for duplicate_card in kw._cards[0]._cards:

        for field in duplicate_card._fields:
            if field.type == float:
                fmt_append = "{:>" + "{0}.5f".format(field.width) + "}"
            elif field.type == int:
                fmt_append = "{:>" + "{0}.0f".format(field.width) + "}"

            line_format = line_format + fmt_append

        line_format = line_format + "\n"

    # format list of strings from data
    list_formatted = ["*KEYWORD\n", kw.get_title() + "\n"]
    elements = kw.elements.to_numpy()
    for element in elements:
        list_formatted.append(line_format.format(*element))
    list_formatted.append("*END")

    # write formatted list of strings
    with open("test_export_elements1.k", "w") as fid:
        fid.writelines(list_formatted)

    t1 = time.time()

    ## original "slow" export for large keywords
    kw_db = Deck()
    kw_db.append(kw)
    kw_db.export_file("test_export_elements.k ")

    t2 = time.time()

    print("Time fast method: {:.3f}s".format(t1 - t0))
    print("Time default method: {:.3f}s".format(t2 - t1))
    print("Speedup: {:.1f}".format((t2 - t1) / (t1 - t0)))

    return


if __name__ == "__main__":

    # keywords for spring b.c.
    kw = keywords.Part()
    kw = keywords.SectionDiscrete()
    kw = keywords.MatSpringElastic()
    kw = keywords.ElementDiscrete()
    kw = keywords.DefineSdOrientation()

    kw = keywords.DampingGlobal()

    example_performance()
    kw = keywords.DefineControlVolume()
    kw = keywords.SetSegmentAdd()
    kw = keywords.ElementShell()

    kw = SetSegmentAdd_custom()

    kw = create_element_shell_keyword(
        np.array([[93, 94, 95], [1, 2, 3]]), part_id=2, id_offset=5
    )

    elements = np.array(
        [[32545, 31655, 28835, 28894], [28099, 28098, 27191, 33154]], dtype=int
    )
    avec = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 1.0]])
    dvec = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 1.0]]) / 2
    element_kw = create_element_solid_ortho_keyword(
        elements=elements, a_vec=avec, d_vec=dvec, element_type="tetra"
    )

    fast_element_writer(element_kw, "test.k")

    nodes = np.array([[0, 0, 0], [0.1, 0.1, 0.1]])

    nodes1 = np.array([[0.2, 0.2, 0.2], [0.3, 0.3, 0.3]])

    kw = add_nodes_to_kw(nodes, kw)
    print(kw)
    kw = add_nodes_to_kw(nodes1, kw)
    print(kw)

    node_kw = keywords.Node()
    nodes = np.array([[-8.0822012582272, 73.5043525327003, 420.0981873906932]])
    kw = add_nodes_to_kw(nodes, node_kw)

    segments = np.array([[1, 2, 3], [2, 3, 1]], dtype=int)
    kw = create_segment_set_keyword(segments, 2, "test")

    # test section ID
    kw_db = Deck()
    kw_db.append(keywords.SectionSolid(secid=1))
    kw_db.append(keywords.SectionShell(secid=2))

    kw_db.append(keywords.MatNull(mid=1))
    kw_db.append(keywords.MatNull(mid=2))

    kw_db = Deck()

    kw_db.append(SetNodeList_custom(sid=1))
    kw_db.append(SetNodeList_custom(sid=2))
    kw_db.append(SetNodeList_custom(sid=5))
    kw_db.append(SetNodeList_custom(sid=6))

    # kw = keywords.SetNodeList_custom( sid = 1)
    kw = create_node_set_keyword(np.arange(0, 10, 1) + 1, 1, "test")

    kw = keywords.SetNodeListSmooth(sid=1)

    kw = keywords.SetSegment()

    # loops over all sections
    ids = get_list_of_used_ids(kw_db, "SET_SEGMENT")
    print(ids)

    ids = get_list_of_used_ids(kw_db, "SET_NODE")
    print(ids)

    print("Protected")
