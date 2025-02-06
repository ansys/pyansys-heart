*KEYWORD
*SET_PART_LIST
$#     sid       da1       da2       da3       da4    solver
         1       0.0       0.0       0.0       0.0MECH      
         1         2         3
*EM_CONTROL_EP
$#olvetype  numsplit   actusig  ionsolvr
         4         1 100000000          
*EM_CONTROL
$#   emsol     numls   macrodt   dimtype    nperio    unused   ncylfem   ncylbem
        11         4       1.0                                    5000          
*EM_EP_ISOCH
$# idisoch    idepol    dplthr    irepol    rplthr cyclenmin    apdmin
         1         1     -20.0         1     -40.0    2000.0     100.0
*EM_SOLVER_FEM
$#  reltol    maxite     stype    precon   uselast  ncyclfem
     1e-06     10000         1         2         1      5000
*EM_OUTPUT
$#    mats      matf      sols      solf      mesh       mem    timing
         1         1         1         1         0         0         0
*SET_NODE_LIST_TITLE
$#                                                                         title
Stim nodes                                                                      
$#     sid       da1       da2       da3       da4    solver       its    unused
        10       0.0       0.0       0.0       0.0MECH                          
$#   nodes     nodes
     20396     27595
*EM_EP_TENTUSSCHER_STIMULUS
$#  stimid   settype     setid
        10         2        10
$#stimstrt     stimt   stimdur   stimamp
       0.0     800.0       2.0      50.0
$ EM-MECA coupling control
*EM_CONTROL_TIMESTEP
$   TSTYPE   DTCONST      LCID    FACTOR     DTMIN     DTMAX
         1       1.0
*EM_CONTROL_COUPLING
$    THCPL     SMCPL    THLCID    SMLCID
         1         0

*END