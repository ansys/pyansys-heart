$-----------------------------------------------------------------------------
$
$ Example provided by Pierre (LSTC)
$X 1. Run file as is.
$X    Requires LS-DYNA MPP Dev 142100 (or higher) with double precision 
$X
$X------------------------------------------------------------------------------
$# UNITS: g/mm/ms
*KEYWORD
*TITLE
EM_EP Benchmark FDA MonoDomain
$---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
$                                                                              $
$                             PARAMETER DEFINITION                             $
$                                                                              $
$---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
*PARAMETER
$--------1---------2---------3---------4---------5---------6---------7---------8
$
$ --- struc ---
$
R       dt         1
R  dt_plot         5     
$
$ --- conductivity ---
$
R      s11   1.33e-4 
R      s22   1.33e-4 
R      s33   1.33e-4
$
$ --- material constants ---
$
R      rho        10
R        E      3.e4
R        v       .45
*INCLUDE
mesh.k
*INCLUDE
struc.k
*INCLUDE
em.k
*END
