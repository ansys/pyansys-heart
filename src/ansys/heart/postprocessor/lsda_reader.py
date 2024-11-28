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

"""Reader classes for LS-DYNA binary output files in Lsda format."""

import glob
import os
from pathlib import Path
import re
import struct
from typing import ByteString, Dict, Iterable, List, Literal, Self, Tuple


class LsdaTypeSizeError(Exception):
    """Raised in case the data type sizes are not what I expect."""

    pass


class FrameNotInFileError(Exception):
    """Raised in case the file does not contain a certain frame."""

    pass


class VarNotInFrameError(Exception):
    """Raised in case the file does not contain a certain variable."""

    pass


class _Diskfile:
    """Handles all the low level file I/O. Nothing here should be called directly by a user."""

    packsize = [0, "b", "h", 0, "i", 0, 0, 0, "q"]
    packtype = [0, "b", "h", "i", "q", "B", "H", "I", "Q", "f", "d", "s"]
    sizeof = [0, 1, 2, 4, 8, 1, 2, 4, 8, 4, 8, 1]

    def __init__(self, name: str):
        """Init Diskfile object for the given file in read mode.

        Parameters
        ----------
        name : str
            str to binary file
        """
        self.name = name  # file name
        self.ateof = 0  # 1 if the file pointer is at EOF
        self.fp = open(name, "rb")
        s = self.fp.read(8)
        header = struct.unpack("BBBBBBBB", s)

        if header[0] > 8:
            self.fp.seek(header[0])

        self.lengthsize = header[1]
        self.offsetsize = header[2]
        self.commandsize = header[3]
        self.typesize = header[4]
        if header[5] == 0:
            self.ordercode = ">"
        else:
            self.ordercode = "<"
        self.ounpack = self.ordercode + _Diskfile.packsize[self.offsetsize]
        self.lunpack = self.ordercode + _Diskfile.packsize[self.lengthsize]
        self.lcunpack = (
            self.ordercode
            + _Diskfile.packsize[self.lengthsize]
            + _Diskfile.packsize[self.commandsize]
        )
        self.tolunpack = (
            self.ordercode
            + _Diskfile.packsize[self.typesize]
            + _Diskfile.packsize[self.offsetsize]
            + _Diskfile.packsize[self.lengthsize]
        )
        self.comp1 = self.typesize + self.offsetsize + self.lengthsize
        self.comp2 = self.lengthsize + self.commandsize + self.typesize + 1

    def readcommand(self):
        """Read a LENGTH,COMMAND pair from the file at the current location."""
        s = self.fp.read(self.lengthsize + self.commandsize)
        return struct.unpack(self.lcunpack, s)

    def readoffset(self):
        """Read an OFFSET from the file at the current location."""
        s = self.fp.read(self.offsetsize)
        return struct.unpack(self.ounpack, s)[0]


class Symbol:
    """A directory tree structure. A Symbol can be a directory (type==0) or data."""

    def __init__(self, name: str = "", parent: str = None):
        """Init Symbol with the given name and parent.

        Parameters
        ----------
        name : str, optional
            Name of the Symbol in the internal file structure, by default ""
        parent : str, optional
            str to the parent in the internal file structure, by default None
        """
        self.name = name  # name of var or directory
        self.type = 0  # data type
        self.offset = 0  # offset of DATA record in file
        self.length = 0  # number of data entries, or # of children
        self.file = 0  # which file the data is in
        self.children = {}  # directory contents
        self.parent = parent  # directory that holds me
        if parent:
            parent.children[name] = self
            parent.length = len(parent.children)

    def path(self):
        """Return absolute path for this Symbol."""
        if not self.parent:
            return "/"
        sym = self
        ret = "/" + sym.name
        while sym.parent and sym.parent.name != "/":
            sym = sym.parent
            ret = "/" + sym.name + ret
        return ret

    def get(self, name: str) -> Self:
        """Return the Symbol with the indicated name.

        The name can be prefixed with a relative or absolute path.

        Parameters
        ----------
        name : str
            Name of the Symbol to recover in the current directory.

        Returns
        -------
        Self
            The requested Symbol or None if it does not exist.
        """
        # If I am just a variable, let my parent handle this
        if self.type != 0:
            return self.parent.get(name)
        # If I have this variable, return it
        if name in self.children.keys():
            return self.children[name]
        # If name has a path component, then look for it there
        if name[0] == "/":  # absolute path
            parts = name.split("/")[1:]
            sym = self
            while sym.parent:
                sym = sym.parent
            for i in range(len(parts)):
                if sym.children.has_key(parts[i]):
                    sym = sym.children[parts[i]]
                else:
                    return None
            return sym
        if name[0] == ".":  # relative path
            parts = name.split("/")[1:]
            # Throw out any "." in the path -- those are just useless....
            parts = filter(lambda p: p != ".", parts)
            if len(parts) == 0:
                return self
            sym = self
            for i in range(parts):
                if parts[i] == "..":
                    if sym.parent:
                        sym = sym.parent
                elif sym.has_key(parts[i]):
                    sym = sym.children[parts[i]]
                else:
                    return None
            return sym
        # Not found
        return None

    def lread(self, start: int = 0, end: int = 2000000000) -> Tuple:
        """Read data from the file. This routine does NOT follow links.

        Parameters
        ----------
        start : int, optional
            start index, by default 0
        end : int, optional
            end index, by default 2000000000

        Returns
        -------
        Tuple
            If this symbol is a DIRECTORY, this returns a sorted list of the
            contents of the directory, and "start" and "end" are ignored.
            Otherwise, read and return data[start:end] (including start but
            not including end -- standard Python slice behavior).
        """
        if self.type == 0:  # directory -- return listing
            return sorted(self.children.keys())
        if end > self.length:
            end = self.length
        if end < 0:
            end = self.length + end
        if start > self.length:
            return ()
        if start < 0:
            start = self.length + start
        if start >= end:
            return ()
        size = _Diskfile.sizeof[self.type]
        pos = self.offset + self.file.comp2 + len(self.name) + start * size
        self.file.fp.seek(pos)
        self.file.ateof = 0
        #    format = self.file.ordercode + _Diskfile.packtype[self.type]*(end-start)
        #    return struct.unpack(format,self.file.fp.read(size*(end-start)))
        format = "%c%d%c" % (
            self.file.ordercode,
            (end - start),
            _Diskfile.packtype[self.type],
        )
        if self.type == Lsda.LINK:
            return struct.unpack(format, self.file.fp.read(size * (end - start)))[0]
        else:
            return struct.unpack(format, self.file.fp.read(size * (end - start)))

    def read(self, start: int = 0, end: int = 2000000000) -> Tuple:
        """Read data from the file.  Same as lread, but follows links.

        Parameters
        ----------
        start : int, optional
            start index, by default 0
        end : int, optional
            end index, by default 2000000000

        Returns
        -------
        Tuple
            If this symbol is a DIRECTORY, this returns a sorted list of the
            contents of the directory, and "start" and "end" are ignored.
            Otherwise, read and return data[start:end] (including start but
            not including end -- standard Python slice behavior).
        """
        return self._resolve_link(self).lread(start, end)

    def read_raw(self, start: int = 0, end: int = 2000000000) -> ByteString | List:
        """Read data from the file and return as bytestring.

        Parameters
        ----------
        start : int, optional
            start index, by default 0
        end : int, optional
            end index, by default 2000000000

        Returns
        -------
        ByteString
            Contents of the variable if self is a variable.
        List
            Directory listing if self is a directory.
        """
        if self.type == 0:  # directory -- return listing
            return sorted(self.children.keys())
        if end > self.length:
            end = self.length
        if end < 0:
            end = self.length + end
        if start > self.length:
            return ()
        if start < 0:
            start = self.length + start
        if start >= end:
            return ()
        size = _Diskfile.sizeof[self.type]
        pos = self.offset + self.file.comp2 + len(self.name) + start * size
        self.file.fp.seek(pos)
        self.file.ateof = 0
        size = size * (end - start)
        return self.file.fp.read(size)

    def _resolve_link(self, var: Self) -> Self:
        """Follow a link to find what it finally resolves to.

        Parameters
        ----------
        var : Self
            Input Symbol.

        Returns
        -------
        Self
            Symbol at the end of link chain.
        """
        ret = var
        while ret.type == Lsda.LINK:
            ret = ret.get(ret.lread())
        return ret


class Lsda:
    """Main class.

    Holds all the Symbols for an LSDA file.
    Has methods for reading data from the file.
    """

    CD = 2
    DATA = 3
    VARIABLE = 4
    BEGINSYMBOLTABLE = 5
    ENDSYMBOLTABLE = 6
    SYMBOLTABLEOFFSET = 7
    I1 = 1
    I2 = 2
    I4 = 3
    I8 = 4
    U1 = 5
    U2 = 6
    U4 = 7
    U8 = 8
    R4 = 9
    R8 = 10
    LINK = 11

    def __init__(self, files: Tuple[str]):
        """Create the LSDA structure, open the file and read the SYMBOLTABLE.

        Parameters
        ----------
        files : Tuple[str]
            Tuple of file names to open. Any missing %XXX continuation files
            will be automatically included.

        Raises
        ------
        LsdaError
            Error raised if the system types don't have the expected byte lengths.
        """
        #
        # If they only input a single name, put it in a tuple, so I can
        # accept input of either kind
        #
        if not types_ok:
            raise LsdaTypeSizeError
        # if type(files) != type((1,)) and type(files) != type([1]):
        if files is not Tuple and files is not List:
            files = (files,)
        self.files = []
        self.fw = None
        # Open all the files in the list that is input, and anything
        # that looks like a continuation of one of them.
        #
        nameset = set()
        for name in files:
            nameset.add(name)
            nameset = nameset.union(set(glob.glob(name + "%[0-9][0-9]*")))
        #
        # Convert to a list and sort
        namelist = list(nameset)
        namelist.sort()
        for file in namelist:
            self.files.append(_Diskfile(file))
        self.root = Symbol("/")
        for f in self.files:
            #
            # We are already positioned to read the SYMBOLTABLEOFFSET record
            #
            (clen, cmd) = f.readcommand()
            self.cwd = self.root
            if cmd == Lsda.SYMBOLTABLEOFFSET:
                self._readsymboltable(f)

        self.cwd = self.root
        self.dirty_symbols = set()
        self.lastpath = None
        #
        # writing will always be to the last one of the files
        #
        self.fw = None
        self.make_dirs = 0

    def cd(self, path: str, create: int = 2) -> Symbol:  # change CWD
        """Change the current working directory in the file.

        Parameters
        ----------
        path : str
            str to the target directory
        create : int, optional
            Flag that pilots on-the-fly directory creation, by default 2

        Returns
        -------
        Symbol
            The Symbol corresponding to the current working directory after
            changing it to the target path. Returns the original cwd if the
            change was unsuccessful.
        """
        if path == "/":
            self.cwd = self.root
            return self.root
        if path[-1] == "/":  # remove trailing /
            path = path[:-1]
        if path[0] == "/":  # absolute path
            path = path[1:]
            self.cwd = self.root
        # path = string.split(path, "/")
        path = path.split("/")
        for part in path:
            if part == "..":
                if self.cwd.parent:
                    self.cwd = self.cwd.parent
            else:
                if part in self.cwd.children.keys():
                    self.cwd = self.cwd.children[part]
                    if self.cwd.type != 0:  # component is a variable, not a directory!
                        self.cwd = self.cwd.parent
                        break
                elif create == 1 or (create == 2 and self.make_dirs == 1):
                    self.cwd = Symbol(part, self.cwd)  # Create directory on the fly
                else:  # component in path is missing
                    break
        return self.cwd

    def get(self, path: str) -> Symbol:
        """Return the Symbol with the indicated name.

        The name can be prefixed with a relative or absolute path.

        Parameters
        ----------
        path : str
            str of the target Symbol.

        Returns
        -------
        Symbol
            Requested Symbol.
        """
        return self.cwd.get(path)

    def _readsymboltable(self, f: _Diskfile):
        """Read all the SYMBOLTABLEs in the current file.

        Users should never call this.

        Parameters
        ----------
        f : _Diskfile
            Target _Diskfile.
        """
        #
        f.ateof = 0
        while 1:
            f.lastoffset = f.fp.tell()
            offset = f.readoffset()
            if offset == 0:
                return
            f.fp.seek(offset)
            (clen, cmd) = f.readcommand()
            if cmd != Lsda.BEGINSYMBOLTABLE:
                return
            while 1:
                (clen, cmd) = f.readcommand()
                clen = clen - f.commandsize - f.lengthsize
                if cmd == Lsda.CD:
                    path = f.fp.read(clen)
                    path = path.decode("utf-8")
                    self.cd(path, 1)
                elif cmd == Lsda.VARIABLE:
                    self._readentry(f, clen, self.cwd)
                else:  # is end of symbol table...get next part if there is one
                    break

    def _readentry(self, f: _Diskfile, reclen: int, parent: Symbol):
        """Read a VARIABLE record from the file, and construct the proper Symbol.

        Users should never call this.

        Parameters
        ----------
        f : _Diskfile
            Target _Diskfile
        reclen : int
            Length to be read
        parent : Symbol
            Parent of the current symbol
        """
        s = f.fp.read(reclen)
        n = reclen - f.comp1
        name = s[:n]
        # If parent already has a symbol by this name, orphan it....
        name = name.decode("utf-8")
        if name in parent.children.keys():
            var = parent.children[name]
        else:
            var = Symbol(name, parent)
        (var.type, var.offset, var.length) = struct.unpack(f.tolunpack, s[n:])
        var.file = f


class LsdaReader:
    """Reader: Wraps Lsda class to extract data from Lsda binary files (e.g. em_binXXXX files).

    Contains methods to extract and rewrite such data.
    """

    def __init__(self, file: Path):
        """Init Reader class for the given file.

        Parameters
        ----------
        file : Path
            Path of the file to be read.
        """
        self.handle = Lsda(file)
        self.frames = self._list_frames()
        self.nids = self._list_nids()
        self.var_dict = self._build_var_dict()

    def _list_frames(self) -> List[str]:
        """List frames contained in the file.

        Returns
        -------
        List[str]
            Names of the frames.
        """
        root = self.handle.root
        items = root.read()
        frames = [item for item in items if re.fullmatch("frame\d{6}", item)]

        return frames

    def _list_nids(self) -> Tuple[int]:
        """Return a tuple of the node ids where the file data was collected.

        Returns
        -------
        Tuple[int]
            Node ids.
        """
        self.handle.cd("nids")
        nids = self.handle.cwd.get("nids").read()
        self.handle.cwd = self.handle.cd("..")

        return nids

    def _build_var_dict(self) -> Dict[str, List[str]]:
        """Generate a dictionary containing available variables for each frame in the file.

        Returns
        -------
        Dict[str:List[str]]
            Variable dictionary. Keys are frame names. Values are variable name lists.
        """
        res = {}

        for frame in self.frames:
            res[frame] = self._list_vars(frame)

        return res

    def _list_vars(self, frame: str) -> List[str]:
        """Return a list with available variables for a given frame.

        Parameters
        ----------
        frame : str
            Target frame name.

        Returns
        -------
        List[str]
            Variable name list for the target frame.
        """
        self.check_frame(frame)
        self.handle.cd(frame)
        vars = self.handle.cwd.read()
        self.handle.cwd = self.handle.cd("..")

        return vars

    def get_var(
        self, frame: str, var: Literal["TMP", "ECP", "ICP", "AT", "RT", "CA2", "time"]
    ) -> Tuple[float]:
        """Return the variable values for a given frame.

        Parameters
        ----------
        frame : str
            Target frame name.
        var : Literal['TMP', 'ECP', 'ICP', 'AT', 'RT', 'CA2', 'time']
            Target variable name.
            TMP: Transmembrane Potential. For emsol=14, this contains the first activation time.
            ECP: Extracellular Potential. For emsol=14, this contains the second activation time.
            ICP: Intracellular Potential. For emsol=14, this contains the third activation time.
            AT: Activation Time. For emsol=14, this contains the fourth activation time.
            RT: Repolarization Time. For emsol=14, this contains the fifth activation time.
            CA2: Calcium Ion Concentration. For emsol=14, this contains the sixth activation time.
            time: Instant of the frame export.

        Returns
        -------
        Tuple[float]
            Requested variable values.
        """
        self.check_frame(frame)
        self.check_var(frame, var)
        self.handle.cd(frame)
        var = self.handle.cwd.get(var).read()
        self.handle.cwd = self.handle.cd("..")

        return var

    def write_frame_to_snapshot(self, path: Path, frame: str, var: str):
        """Write the variable values in the frame into a binary file in Ansys ROM builder format.

        Parameters
        ----------
        path : Path
            Path for writing the file.
        frame : str
            Target frame name.
        var : str
            Target variable name.
        """
        vec = self.get_var(frame, var)
        self._write_binary(path, vec)

    def write_frames_to_snapshots(
        self, dir: Path, var: str, prefix: str = "snapshot", extension: str = ".bin"
    ):
        """Write the requested variable values in for all available frames to binary files.

        One binary snapshot is created per frame.
        Snapshots are written in the given directory in Ansys ROM builder format.
        Snapshot file names are as follows: <prefix><frame_id><extension>.


        Parameters
        ----------
        dir : Path
            Directory for writing snapshots.
        var : str
            Target variable.
        prefix : str, optional
            File output prefix, by default "snapshot"
        extension : str, optional
            File output extension, by default ".bin"
        """
        for frame in self.frames:
            frame_id = re.sub("frame0{,5}", "", frame)
            current_path = os.path.join(dir, f"{prefix}{frame_id}{extension}")
            self.write_frame_to_snapshot(current_path, frame, var)

    def check_frame(self, frame: str):
        """Check if the given frame was loaded into self.frames.

        Parameters
        ----------
        frame : str
            Frame name

        Raises
        ------
        FrameNotInFileError
            Placeholder exception for frames that don't exist.
        """
        if frame not in self.frames:
            raise FrameNotInFileError("{frame} does not exist in this file.")

    def check_var(self, frame: str, var: str):
        """Check if the given var exists in the given frame using self.var_dict.

        Parameters
        ----------
        frame : str
            Frame name.
        var : str
            Variable name.

        Raises
        ------
        VarNotInFrameError
            Placeholder exception for variables that don't exist in a given frame.
        """
        if var not in self.var_dict[frame]:
            raise VarNotInFrameError("{var} does not exist in {frame}.")

    def _write_binary(self, fp: Path, vector: Iterable):
        """Write vector to binary snapshot format usable in Ansys ROM builders.

        Parameters
        ----------
        fp : Path
            Path of the file to write.
        vector : Iterable
            Vector to be contained in the file.
        """
        with open(fp, mode="wb") as f:
            f.write(struct.pack("Q", len(vector)))

            for item in vector:
                f.write(struct.pack("d", item))


# Do sanity check of type lengths.  types_ok will be checked whenever
# a new file is opened, and if things don't match then an exception will
# be raised -- this should prevent unknown errors due to type size problems
#
types = [("b", 1), ("h", 2), ("i", 4), ("q", 8), ("f", 4), ("d", 8)]
x = 17
types_ok = 1
for a, b in types:
    s = struct.pack(a, x)
    if len(s) != b:
        print("LSDA: initialization error")
        print("Data type %s has length %d instead of %d" % (a, len(s), b))
        types_ok = 0

if __name__ == "__main__":
    pass
