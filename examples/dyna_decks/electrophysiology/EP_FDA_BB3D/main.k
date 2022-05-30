$-----------------------------------------------------------------------------
$
$ Example provided by Pierre (LSTC)
$X 1. Run file as is.
$X    Requires LS-DYNA MPP Dev 142100 (or higher) with double precision 
$X
$X-----------------------------------------------------------------------------
$# UNITS: (Dimensionless)
$X-----------------------------------------------------------------------------
*KEYWORD
*TITLE
EM_EP Benchmark FDA BiDomain w/ Bath
$--------1---------2---------3---------4---------5---------6---------7---------8
$                                                                              $
$                             PARAMETER DEFINITION                             $
$                                                                              $
$--------1---------2---------3---------4---------5---------6---------7---------8
*PARAMETER
$--------1---------2---------3---------4---------5---------6---------7---------8
$
$ --- struc ---
$
R       dt    0.0025
R  dt_plot       0.1
$
$ --- conductivity ---
$
R      s11 0.1114533
R      s22 0.1215854
R      s33 0.0303964
R     s11b 0.0464164
R     s22b  0.050361
R     s33b   0.01259
R   s_bath  0.023082
$
$ --- material constants ---
$
R      rho   8928.57
R        E    200.e9
R        v        .3
*INCLUDE
mesh.k
*INCLUDE
struc.k
*INCLUDE
em.k
*END
