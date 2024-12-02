*KEYWORD
*TITLE
FullHeart
$$- Unit system: g-mm-ms-N-MPa-mJ -$$
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas     nosol
     800.0         0       0.0       0.0     1e+08         0
*CONTROL_IMPLICIT_DYNAMICS
$#   imass     gamma      beta    tdybir    tdydth    tdybur     irate     alpha
         1       0.6      0.38       0.0     1e+28     1e+28         0       0.0
$$ Disable auto step due 0D model $$
*CONTROL_IMPLICIT_AUTO
$#   iauto    iteopt    itewin     dtmin     dtmax     dtexp     kfail    kcycle
         0        11         5       5.0       5.0                              
*CONTROL_IMPLICIT_GENERAL
$#  imflag       dt0    imform      nsbs       igs     cnstn      form    zero_v
         1       5.0         2         1         2         0         0         0
*CONTROL_IMPLICIT_SOLUTION
$#  nsolvr    ilimit    maxref     dctol     ectol     rctol     lstol    abstol
        12        11        15      0.02 1000000.0    1000.0      -0.9    -1e-20
$#   dnorm    diverg     istif   nlprint    nlnorm   d3itctl     cpchk
         1         1         1         3       4.0         0         0
$#  arcctl    arcdir    arclen    arcmth    arcdmp    arcpsi    arcalf    arctim
         0         0       0.0         1         2       0.0       0.0       0.0
$#   lsmtd     lsdir      irad      srad      awgt      sred
         5         2       0.0       0.0       0.0       0.0
*CONTROL_IMPLICIT_SOLVER
$#  lsolvr    lprint     negev     order      drcm    drcprm   autospc   autotol
         2         0         2         0         4                   2          
$#  lcpack    mtxdmp    iparm1    rparm1    rparm2                        rparm5
         2         0       500     1e-09     0.001                           0.0
$#  emxdmp    rdcmem    absmem
         0      0.85          
*CONTROL_ACCURACY
$#     osu       inn    pidosu      iacc     exacc
         1         4                   1          
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
*DATABASE_ICVOUT
$#      dt    binary      lcur     ioopt
       5.0         2         0         1
*DATABASE_ABSTAT
$#      dt    binary      lcur     ioopt
       5.0         2         0         1
*DATABASE_GLSTAT
$#      dt    binary      lcur     ioopt
       0.1         2         0         1
*DATABASE_MATSUM
$#      dt    binary      lcur     ioopt
       0.1         2         0         1
*DATABASE_BINARY_D3PLOT
$#      dt      lcdt      beam     npltc    psetid
      20.0                   0                    
$#   ioopt      rate    cutoff    window      type      pset
         0                                       0         0
*DATABASE_EXTENT_BINARY
$#   neiph     neips    maxint    strflg    sigflg    epsflg    rltflg    engflg
        27         0         0         1         1         1         1         1
$#  cmpflg    ieverp    beamip     dcomp      shge     stssz    n3thdt   ialemat
         0         0         0         1         1         1         2         1
$# nintsld   pkp_sen      sclp     hydro     msscl     therm    intout    nodout
         0         0       1.0         0         0         0                    
$#    dtdt    resplt     neipb   quadsld    cubsld   deleres
         0         1                   0         0         0
*DAMPING_GLOBAL
$#    lcid    valdmp       stx       sty       stz       srx       sry       srz
         1       0.0       0.0       0.0       0.0       0.0       0.0       0.0
*DEFINE_CURVE_TITLE
$#                                                                         title
global damping [ms^-1]                                                          
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         1         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1
                 0.00.100000000000000006
               1e+260.100000000000000006
$$ Left ventricle stiffness damping [ms]
*DAMPING_PART_STIFFNESS
$#     pid      coef
         1      -0.2
$$ Right ventricle stiffness damping [ms]
*DAMPING_PART_STIFFNESS
$#     pid      coef
         2      -0.2
$$ Septum stiffness damping [ms]
*DAMPING_PART_STIFFNESS
$#     pid      coef
         7      -0.2
$$ Left atrium stiffness damping [ms]
*DAMPING_PART_STIFFNESS
$#     pid      coef
         3      -0.2
$$ Right atrium stiffness damping [ms]
*DAMPING_PART_STIFFNESS
$#     pid      coef
         4      -0.2
$$ Aorta stiffness damping [ms]
*DAMPING_PART_STIFFNESS
$#     pid      coef
         5      -0.2
$$ Pulmonary artery stiffness damping [ms]
*DAMPING_PART_STIFFNESS
$#     pid      coef
         6      -0.2
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
control_volume.k                                                                
*INCLUDE
$#                                                                      filename
pericardium.k                                                                   
*END