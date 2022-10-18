*KEYWORD
*TITLE
BiVentricle zero-pressure
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas     nosol
       1.0         0       0.0       0.0     1e+08         0
*CONTROL_IMPLICIT_DYNAMICS
$#   imass     gamma      beta    tdybir    tdydth    tdybur     irate     alpha
         0       0.5      0.25       0.0     1e+28     1e+28         0       0.0
*CONTROL_IMPLICIT_AUTO
$#   iauto    iteopt    itewin     dtmin     dtmax     dtexp     kfail    kcycle
         1        11         5      0.01       0.1                              
*CONTROL_IMPLICIT_GENERAL
$#  imflag       dt0    imform      nsbs       igs     cnstn      form    zero_v
         1       0.1         2         1         2         0         0         0
*CONTROL_IMPLICIT_SOLUTION
$#  nsolvr    ilimit    maxref     dctol     ectol     rctol     lstol    abstol
        12        11        15     0.001      0.01     1e+10       0.9     1e-10
$#   dnorm    diverg     istif   nlprint    nlnorm   d3itctl     cpchk
         2         1         1         0       2.0         0         0
$#  arcctl    arcdir    arclen    arcmth    arcdmp    arcpsi    arcalf    arctim
         0         0       0.0         1         2       0.0       0.0       0.0
$#   lsmtd     lsdir      irad      srad      awgt      sred
         4         2       0.0       0.0       0.0       0.0
*CONTROL_IMPLICIT_SOLVER
$#  lsolvr    lprint     negev     order      drcm    drcprm   autospc   autotol
         2         0         2         0         4                   1          
$#  lcpack    mtxdmp    iparm1    rparm1    rparm2                        rparm5
         2         0       500     1e-09     0.001                           0.0
$#  emxdmp    rdcmem    absmem
         0      0.85          
*DATABASE_NODOUT
$#      dt    binary      lcur     ioopt   option1   option2
       0.2         1         0         1       0.0         0
*SET_NODE_GENERAL
$#     sid       da1       da2       da3       da4    solver       its
         1       0.0       0.0       0.0       0.0MECH      1         
$#  option        e1        e2        e3        e4        e5        e6        e7
ALL                                                                             
*DATABASE_HISTORY_NODE_SET
$#     id1       id2       id3       id4       id5       id6       id7       id8
         1                                                                      
*CONTROL_OUTPUT
$#   npopt    neecho    nrefup    iaccop     opifs    ipnint    ikedit    iflush
         1         1         0         0       0.0         0         0         0
$#   iprtf    ierode   tet10s8    msgmax    ipcurv      gmdt   ip1dblt      eocs
         0         0         2        50         0       0.0         0         0
$#   tolev    newleg    frfreq     minfo    solsig    msgflg    cdetol
         2         0         1         0         0         0      10.0
$# phschng    demden   icrfile   spc2bnd    penout    shlsig   hisnout    engout
         0         0         0                   0         0         0         0
$#    insf    isolsf      ibsf      issf    mlkbag
         0         0         0         0         0
*DATABASE_BINARY_D3PLOT
$#      dt      lcdt      beam     npltc    psetid
       0.5                   0                    
$#   ioopt      rate    cutoff    window      type      pset
         0                                       0         0
*DEFINE_CURVE_TITLE
$#                                                                         title
unit load curve                                                                 
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         1         0       1.0       1.0       0.0       0.0         0       100
$#                a1                  o1
                 0.0                 0.0
                 1.0                 1.0
 1.00099999999999989                 0.0
*LOAD_SEGMENT_SET
$#    ssid      lcid        sf        at
         1         1       2.0       0.0
*LOAD_SEGMENT_SET
$#    ssid      lcid        sf        at
         2         1   0.53333       0.0
*CONTROL_REFERENCE_CONFIGURATION
$# maxiter                                                                target
         3nodes.k                                                               
$#  method      step       tol
         2       1.0       5.0
*INCLUDE
$#                                                                      filename
parts.k                                                                         
*INCLUDE
$#                                                                      filename
nodes.k                                                                         
*INCLUDE
$#                                                                      filename
solid_elements.k                                                                
*INCLUDE
$#                                                                      filename
material.k
*INCLUDE
$#                                                                      filename
segment_sets.k                                                                  
*INCLUDE
$#                                                                      filename
node_sets.k                                                                     
*INCLUDE
$#                                                                      filename
boundary_conditions.k                                                           
*INCLUDE
$#                                                                      filename
cap_elements.k                                                                  
*INCLUDE
$#                                                                      filename
pericardium.k                                                                   
*END