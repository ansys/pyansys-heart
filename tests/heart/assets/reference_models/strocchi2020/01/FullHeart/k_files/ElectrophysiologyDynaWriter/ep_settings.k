*KEYWORD
*EM_CONTROL
$#   emsol     numls   macrodt   dimtype    nperio    unused   ncylfem   ncylbem
        11         4       1.0                                    5000          
*EM_EP_ISOCH
$# idisoch    idepol    dplthr    irepol    rplthr
         1         1     -20.0         1     -40.0
*EM_CONTROL_EP
$#olvetype  numsplit   actusig
         4         1 100000000
*EM_SOLVER_FEM
$#  reltol    maxite     stype    precon   uselast  ncyclfem
     1e-06     10000         1         2         1      5000
*EM_OUTPUT
$#    mats      matf      sols      solf      mesh       mem    timing
         1         1         1         1         0         0         0
*SET_NODE_LIST_TITLE
$#                                                                         title
Stim nodes                                                                      
$#     sid       da1       da2       da3       da4    solver
        16       0.0       0.0       0.0       0.0MECH      
$#      n1        n2        n3        n4        n5        n6        n7        n8
      9337     27833                                                            
*EM_EP_TENTUSSCHER_STIMULUS
$#  stimid   settype     setid
         1         2        16
$#stimstrt     stimt   stimdur   stimamp
       0.0     800.0      20.0      50.0
*END