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

import os
import warnings

if os.name == "nt":
    supported_versions = ["222"]
    awp_roots = {ver: os.environ.get(f"AWP_ROOT{ver}", "") for ver in supported_versions}
    installed_versions = {
        ver: path for ver, path in awp_roots.items() if path and os.path.isdir(path)
    }

elif os.name == "posix":
    # assert False, "Posix not supported yet"
    warnings.warn("Checks for posix disabled", Warning)
    pass

if os.name == "posix":
    warnings.warn("Skipping product installation checks", Warning)

elif os.name == "nt" and installed_versions:
    # use latest installed version that is supported and has Spaceclaim and Fluent
    SC_EXE = ""
    FLUENT_EXE = ""
    for installed_version, ansys_dir in sorted(installed_versions.items(), reverse=True):
        sc_exe = os.path.join(ansys_dir, "scdm", "SpaceClaim.exe")
        fluent_exe = os.path.join(ansys_dir, "fluent", "ntbin", "win64", "fluent.exe")
        if os.path.isfile(sc_exe) and os.path.isfile(fluent_exe):
            SC_EXE = sc_exe
            FLUENT_EXE = fluent_exe
            break

    if not os.path.isfile(SC_EXE):
        warnings.warn("Spaceclaim not found", Warning)
    if not os.path.isfile(FLUENT_EXE):
        warnings.warn("Fluent not found", Warning)

else:
    warnings.warn(
        "No valid Ansys installations found. Valid Ansys installations include: Ansys %s"
        % supported_versions,
        Warning,
    )

import os
