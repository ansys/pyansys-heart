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
    # use latest installed version that is supported and has Spaceclaim and Fluent
    SC_EXE = ""
    FLUENT_EXE = ""
    for installed_version, ansys_dir in sorted( installed_versions.items(), reverse = True ):
        sc_exe = os.path.join(ansys_dir, "scdm", "SpaceClaim.exe")
        fluent_exe = os.path.join(ansys_dir, "fluent", "ntbin", "win64", "fluent.exe")        
        if os.path.isfile(sc_exe) and os.path.isfile(sc_exe):
            SC_EXE     = sc_exe
            FLUENT_EXE = fluent_exe
            break

    assert os.path.isfile( SC_EXE ), "Spaceclaim not found"
    assert os.path.isfile( FLUENT_EXE ), "Fluent not found"
else:
    assert False, "No valid installations found. Valid Ansys installations include: Ansys %s" % supported_versions
    exit()


