$-----------------------------------------------------------------------------
$
$ Example provided by Pierre (LSTC)
$X 1. Run file as is.
$X    Requires LS-DYNA MPP Dev 142100 (or higher) with double precision
$X
$X------------------------------------------------------------------------------
$# UNITS: g/mm/ms
*KEYWORD
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
R       dt       0.1
R  dt_plot         1
$
$ --- conductivity ---
$
R      s11    0.1e-3
R      s22    0.1e-3
R      s33    0.1e-3
R    s11_2  0.166e-3
R    s22_2  0.166e-3
R    s33_2  0.166e-3
R   s_bath     1.e-3
R     beta      0.14
R       Cm      0.01
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
