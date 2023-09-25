*KEYWORD
$--------1---------2---------3---------4---------5---------6---------7---------8
$                                                                              $
$                             EM CONTROL                                       $
$                                                                              $
$--------1---------2---------3---------4---------5---------6---------7---------8
*EM_CONTROL
$--------1---------2---------3---------4---------5---------6---------7---------8
$    emsol     numls  emdtinit   emdtmax   emtinit    emtend  ncyclFem  ncyclBem
        11         4       &dt
*EM_SOLVER_FEM
$--------1---------2---------3---------4---------5---------6---------7---------8
$   relTol   maxIter solveType   precond   uselast
     1.e-9     10000         1         2         1
*EM_SOLVER_FEM
$--------1---------2---------3---------4---------5---------6---------7---------8
$   relTol   maxIter solveType   precond   uselast
     1.e-9     10000         1         2         1
$--------1---------2---------3---------4---------5---------6---------7---------8
$                                                                              $
$                             EM MAT                                           $
$                                                                              $
$--------1---------2---------3---------4---------5---------6---------7---------8
*EM_MAT_003
$--------1---------2---------3---------4---------5---------6---------7---------8
$#intramid     mtype   sigmaXX   sigmaYY   sigmaZZ      beta        Cm
         1         2      &s11      &s22      &s33     &beta       &Cm
$--------1---------2---------3---------4---------5---------6---------7---------8
$# sigmaXY   sigmaXZ   sigmaYX   sigmaYZ   sigmaZX   sigmaZY

$--------1---------2---------3---------4---------5---------6---------7---------8
$#    AOPT        XP        YP        ZP        A1        A2        A3      MACF
        2.                             0         0         0         1
$--------1---------2---------3---------4---------5---------6---------7---------8
$#      V1        V2        V3        D1        D2        D3
                                       0        -1         0
*EM_MAT_003
$--------1---------2---------3---------4---------5---------6---------7---------8
$#intramid     mtype   sigmaXX   sigmaYY   sigmaZZ      beta        Cm
         2         2    &s11_2    &s22_2    &s33_2     &beta       &Cm
$--------1---------2---------3---------4---------5---------6---------7---------8
$# sigmaXY   sigmaXZ   sigmaYX   sigmaYZ   sigmaZX   sigmaZY

$--------1---------2---------3---------4---------5---------6---------7---------8
$#    AOPT        XP        YP        ZP        A1        A2        A3      MACF
        2.                             0         0         0         1
$--------1---------2---------3---------4---------5---------6---------7---------8
$#      V1        V2        V3        D1        D2        D3
                                       0        -1         0
$--------1---------2---------3---------4---------5---------6---------7---------8
$                                                                              $
$                             EM TENTUSSCHER STIMULUS                          $
$                                                                              $
$--------1---------2---------3---------4---------5---------6---------7---------8
*EM_EP_TENTUSSCHER
$--------1---------2---------3---------4---------5---------6---------7---------8
$    matid
         1
$ CONSTS
$--------1---------2---------3---------4---------5---------6---------7---------8
$  consts1   consts2   consts3   consts4   consts5  consts52  consts53  consts10
$        R         T         F        Cm        Vc       Vsr       Vss     p_KNa
  8314.472       31096485.3415     0.185  0.016404  0.0010940.00005468      0.03
$--------1---------2---------3---------4---------5---------6---------7---------8
$ consts11  consts12  consts13
$      K_o      Na_o      Ca_o
       5.4      140.        2.
$--------1---------2---------3---------4---------5---------6---------7---------8
$ consts14  consts15  consts16  consts17  consts18  consts19  consts20  consts21
$     g_K1      g_Kr      g_Ks      g_Na     g_bNa     g_CaL     g_bCa      g_to
     5.405     0.153     0.392    14.838   0.0002  0.0000398  0.000592     0.294
$ consts31  consts33
$    g_pCa      g_pK
    0.1238    0.0146
$--------1---------2---------3---------4---------5---------6---------7---------8
$ consts22  consts23  consts24  consts25  consts26  consts27  consts28  consts29
$    P_NaK      K_mK     K_mNa    k_NaCa      ksat     alpha     gamma     Km_Ca
     2.724        1.       40.     1000.       0.1       2.5      0.35      1.38
$ consts30  consts32
$   Km_Nai     K_pCa
      87.5    0.0005
$--------1---------2---------3---------4---------5---------6---------7---------8
$ consts34  consts35  consts36  consts37  consts38  consts39  consts40
$      k1'       k2'        k3        k4        EC    max_sr    min_sr
      0.15     0.045      0.06     0.005       1.5       2.5        1.
$ consts41  consts44  consts42  consts45  consts43
$    V_rel    V_leak    V_xfer   Vmax_up      K_up
     0.102   0.00036    0.0038  0.006375   0.00025
$ consts46  consts47  consts48  consts49  consts50  consts51
$    Buf_c   K_buf_c    Buf_sr  K_buf_sr    Buf_ss  K_buf_ss
       0.2     0.001       10.       0.3       0.4   0.00025
$ INIT STATES
$--------1---------2---------3---------4---------5---------6---------7---------8
$  states1   states2   states3   states4  states11  states18  states19
$        V       K_i      Na_i      Ca_i     Ca_ss     Ca_sr        R'
    -85.23    136.89     8.604  0.000126   0.00036      3.64    0.9073
$--------1---------2---------3---------4---------5---------6---------7---------8
$  states5   states6   states7   states8   states9  states10  states12  states13
$      xr1       xr2        xs         m         h         j         d         f
   0.00621    0.4712    0.0095   0.00172    0.7444    0.7045  3.373e-5    0.7888
$ states14  states15  states16  states17
$       f2     fCass         s         r
    0.9755    0.9953  0.999998   2.42e-8
*EM_EP_TENTUSSCHER
$--------1---------2---------3---------4---------5---------6---------7---------8
$    matid
         2
$ CONSTS
$--------1---------2---------3---------4---------5---------6---------7---------8
$  consts1   consts2   consts3   consts4   consts5  consts52  consts53  consts10
$        R         T         F        Cm        Vc       Vsr       Vss     p_KNa
  8314.472       31096485.3415     0.185  0.016404  0.0010940.00005468      0.03
$--------1---------2---------3---------4---------5---------6---------7---------8
$ consts11  consts12  consts13
$      K_o      Na_o      Ca_o
       5.4      140.        2.
$--------1---------2---------3---------4---------5---------6---------7---------8
$ consts14  consts15  consts16  consts17  consts18  consts19  consts20  consts21
$     g_K1      g_Kr      g_Ks      g_Na     g_bNa     g_CaL     g_bCa      g_to
     5.405     0.153     0.392    14.838   0.0002  0.0000398  0.000592     0.294
$ consts31  consts33
$    g_pCa      g_pK
    0.1238    0.0146
$--------1---------2---------3---------4---------5---------6---------7---------8
$ consts22  consts23  consts24  consts25  consts26  consts27  consts28  consts29
$    P_NaK      K_mK     K_mNa    k_NaCa      ksat     alpha     gamma     Km_Ca
     2.724        1.       40.     1000.       0.1       2.5      0.35      1.38
$ consts30  consts32
$   Km_Nai     K_pCa
      87.5    0.0005
$--------1---------2---------3---------4---------5---------6---------7---------8
$ consts34  consts35  consts36  consts37  consts38  consts39  consts40
$      k1'       k2'        k3        k4        EC    max_sr    min_sr
      0.15     0.045      0.06     0.005       1.5       2.5        1.
$ consts41  consts44  consts42  consts45  consts43
$    V_rel    V_leak    V_xfer   Vmax_up      K_up
     0.102   0.00036    0.0038  0.006375   0.00025
$ consts46  consts47  consts48  consts49  consts50  consts51
$    Buf_c   K_buf_c    Buf_sr  K_buf_sr    Buf_ss  K_buf_ss
       0.2     0.001       10.       0.3       0.4   0.00025
$ INIT STATES
$--------1---------2---------3---------4---------5---------6---------7---------8
$  states1   states2   states3   states4  states11  states18  states19
$        V       K_i      Na_i      Ca_i     Ca_ss     Ca_sr        R'
    -85.23    136.89     8.604  0.000126   0.00036      3.64    0.9073
$--------1---------2---------3---------4---------5---------6---------7---------8
$  states5   states6   states7   states8   states9  states10  states12  states13
$      xr1       xr2        xs         m         h         j         d         f
   0.00621    0.4712    0.0095   0.00172    0.7444    0.7045  3.373e-5    0.7888
$ states14  states15  states16  states17
$       f2     fCass         s         r
    0.9755    0.9953  0.999998   2.42e-8
*EM_EP_TENTUSSCHER_STIMULUS
$--------1---------2---------3---------4---------5---------6---------7---------8
$   stimid   setType      nsid
         1         2         2
$ constants for stimulus
$--------1---------2---------3---------4---------5---------6---------7---------8
$  consts6   consts7   consts8   consts9
$stimStart     stimT   stimDur   stimAmp
        0.    10000.        2.       50.
$--------1---------2---------3---------4---------5---------6---------7---------8
$                                                                              $
$                             EM OUTPUT                                        $
$                                                                              $
$--------1---------2---------3---------4---------5---------6---------7---------8
*EM_OUTPUT
$--------1---------2---------3---------4---------5---------6---------7---------8
$     matS      matF      solS      solF      mesh    memory    timing d3plotAsc
         2         2         2         2         0         0         0         0
$      mf2       gmv d3plotFor  timeHist spoolesTr    intDat
         0         0         0         0         0         0
*END
