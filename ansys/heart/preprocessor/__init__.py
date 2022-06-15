import os

if os.name == "nt":
    supported_versions = ["212", "221"]
    awp_roots = {
        ver: os.environ.get(f"AWP_ROOT{ver}", "") for ver in supported_versions
    }
    installed_versions = {
        ver: path for ver, path in awp_roots.items() if path and os.path.isdir(path)
    }    

elif os.name == "posix":
    assert False, "Posix not supported yet"

if installed_versions:    
    # get latest installed version:
    ansys_dir = installed_versions[ sorted(installed_versions.keys() )[-1] ]

    SC_EXE = os.path.join(ansys_dir, "scdm", "SpaceClaim.exe")
    FLUENT_EXE = os.path.join(ansys_dir, "fluent", "ntbin", "win64", "fluent.exe")
    assert os.path.isfile( SC_EXE ), "Spaceclaim not found"
    assert os.path.isfile( FLUENT_EXE ), "Fluent not found"
else:
    assert False, "No valid installations found. Valid Ansys installations include: Ansys %s" % supported_versions
    exit()


