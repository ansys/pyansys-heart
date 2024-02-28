*KEYWORD
*CONTROL_SOLUTION
$ thermal only
$     soln       nlq     isnan     lcint     lcacc     ncdcf
         1
*CONTROL_THERMAL_SOLVER
$ steady-state thermal analysis
$    atype     ptype    solver     cgtol       gpt    eqheat     fwork       sbc
         0         0         1
*CONTROL_TERMINATION
$   endtim    endcyc     dtmin    endeng    endmas     nosol
        1.
$
$================================ OUTPUTDATA ===================================
$
*DATABASE_BINARY_D3PLOT
$       dt
       1.0
*DATABASE_GLSTAT
$       dt
       1.0
*DATABASE_MATSUM
$       dt
       1.0
*DATABASE_TPRINT
$       dt
       1.0
$
$============================== PART DEFINITIONS ===============================
$
*PART
$ title

$      pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         1         1         0         0         0         0         0         1
$
$============================= SECTION PROPERTIES ==============================
$
*SECTION_SOLID
$    secid    elform
         1         1
$
$======================== THERMAL MATERIAL PROPERTIES ==========================
$
*MAT_THERMAL_ISOTROPIC
$     tmid       tro     tgrlc    tgmult      tlat      hlat
         1    1.0e-9
$       hc        tc
    1.0000       1.0
$
$======================== NODE/ELEMENT/SET DEFINITIONS ==========================
$
*INCLUDE
nodes.k
solid_elements.k
node_sets.k
*CASE
1,UVC_LONGITUDIANL

1
*CASE
2,UVC_TRANSMURAL

2
*CASE
3,UVC_ROTATIONAL

3
$==================== THERMAL BOUNDARY and INITIAL CONDITIONS ==================
$
*CASE_BEGIN_1
*BOUNDARY_TEMPERATURE_SET
$  NID/SID     TLCID     TMULT       LOC    TDEATH    TBIRTH
         2         0       0.0
         1         0       1.0
*CASE_END_1
$
*CASE_BEGIN_2
*BOUNDARY_TEMPERATURE_SET
$  NID/SID     TLCID     TMULT       LOC    TDEATH    TBIRTH
        20         0       0.0
        30         0       1.0
*CASE_END_2
$
*CASE_BEGIN_3
*BOUNDARY_TEMPERATURE_SET
$  NID/SID     TLCID     TMULT       LOC    TDEATH    TBIRTH
       100         0     -3.14
       200         0     +3.14
       300         0      0.00
*CASE_END_3
$
*END
