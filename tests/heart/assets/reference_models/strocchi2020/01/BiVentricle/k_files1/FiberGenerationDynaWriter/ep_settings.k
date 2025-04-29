*KEYWORD
*EM_CONTROL
$#   emsol     numls   macrodt   dimtype    nperio    unused   ncylfem   ncylbem
        11         4       1.0                                    5000          
*EM_CONTROL_TIMESTEP
$#  tstype    dtcons      lcid    factor     tsmin     tsmas     rlcsf    mecats
         1       1.0                 1.0                            25         0
*EM_CONTROL_EP
$#olvetype  numsplit   actusig  ionsolvr
         4         1 100000000          
*EM_SOLVER_FEM
$#  reltol    maxite     stype    precon   uselast  ncyclfem
     1e-06     10000         1         2         1      5000
*EM_OUTPUT
$#    mats      matf      sols      solf      mesh       mem    timing
         1         1         1         1         0         0         0
*END