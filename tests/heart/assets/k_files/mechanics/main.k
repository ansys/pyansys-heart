*KEYWORD
*TITLE
BiVentricle
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas     nosol
       5.0         0     0.001       0.0     1e+08         0
*CONTROL_IMPLICIT_DYNAMICS
$#   imass     gamma      beta    tdybir    tdydth    tdybur     irate     alpha
         1       0.6      0.38       1.0     1e+28     1e+28         0       0.0
*DEFINE_CURVE_TITLE
$#                                                                         title
time step control                                                               
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         1         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1
                 0.00.050000000000000003
                 1.00.050000000000000003
 1.01000000000000001               0.001
                 5.0                0.01
*CONTROL_IMPLICIT_AUTO
$#   iauto    iteopt    itewin     dtmin     dtmax     dtexp     kfail    kcycle
         1        11         5     0.001      -1.0                              
*CONTROL_IMPLICIT_GENERAL
$#  imflag       dt0    imform      nsbs       igs     cnstn      form    zero_v
         1      0.01         2         1         2         0         0         0
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
     0.001         2         0         1
*DATABASE_ABSTAT
$#      dt    binary      lcur     ioopt
     0.001         2         0         1
*DATABASE_ELOUT
$#      dt    binary      lcur     ioopt   option1   option2   option3   option4
       0.1         2         0         1         0         0         0         0
*DATABASE_GLSTAT
$#      dt    binary      lcur     ioopt
       0.1         2         0         1
*DATABASE_MATSUM
$#      dt    binary      lcur     ioopt
       0.1         2         0         1
*DEFINE_CURVE_TITLE
$#                                                                         title
d3plot out control                                                              
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         2         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1
                 0.0                 0.5
                 1.0                 0.5
 1.050000000000000040.050000000000000003
                 5.00.050000000000000003
*DATABASE_BINARY_D3PLOT
$#      dt      lcdt      beam     npltc    psetid
      0.05         2         0                    
$#   ioopt      rate    cutoff    window      type      pset
         1                                       0         0
*DATABASE_EXTENT_BINARY
$#   neiph     neips    maxint    strflg    sigflg    epsflg    rltflg    engflg
        27         0         0         1         1         1         1         1
$#  cmpflg    ieverp    beamip     dcomp      shge     stssz    n3thdt   ialemat
         0         0         0         1         1         1         2         1
$# nintsld   pkp_sen      sclp     hydro     msscl     therm    intout    nodout
         0         0       1.0         0         0         0                    
$#    dtdt    resplt     neipb   quadsld    cubsld   deleres
         0         0                   0         0         0
*DAMPING_GLOBAL
$#    lcid    valdmp       stx       sty       stz       srx       sry       srz
         3       0.0       0.0       0.0       0.0       0.0       0.0       0.0
*DEFINE_CURVE_TITLE
$#                                                                         title
damping                                                                         
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         3         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1
                 0.0               500.0
               1e+26               500.0
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