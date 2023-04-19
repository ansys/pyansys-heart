"""Module containing miscellaneous methods that don't fit other modules."""


def add_solid_name_to_stl(filename, solid_name, file_type: str = "ascii") -> None:
    """Add name of solid to stl file.
    Note
    ----
    Supports only single block.
    """
    if file_type == "ascii":
        start_str = "solid"
        end_str = "endsolid"
        f = open(filename, "r")
        list_of_lines = f.readlines()
        f.close()
        list_of_lines[0] = "{0} {1}\n".format(start_str, solid_name)
        list_of_lines[-1] = "{0} {1}\n".format(end_str, solid_name)

        f = open(filename, "w")
        f.writelines(list_of_lines)
        f.close()
    # replace part name in binary file
    elif file_type == "binary":
        fid = open(filename, "r+b")
        fid.seek(0)
        data = fid.read(40)
        fid.seek(0)
        string_replace = "{:<40}".format(solid_name).encode()
        fid.write(string_replace)
        fid.close()
    return
