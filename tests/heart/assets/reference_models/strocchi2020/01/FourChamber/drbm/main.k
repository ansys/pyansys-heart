*KEYWORD
*CONTROL_SOLUTION
$#    soln       nlq     isnan     lcint     lcacc     ncdcf     nocop
         1                   0       100         0         1         0
*CONTROL_THERMAL_SOLVER
$#   atype     ptype    solver    unused       gpt    eqheat     fwork       sbc
         0         0        11                   8       1.0       1.0       0.0
$#  msglvl    maxitr    abstol    reltol     omega    unused    unused       tsf
         0       500     1e-10     1e-06       0.0                           1.0
$#   mxdmp      dtvf    varden    unused     ncycl
         0       0.0         0                   1
*DATABASE_BINARY_D3PLOT
$#      dt      lcdt      beam     npltc    psetid
       1.0                   0                    
$#   ioopt      rate    cutoff    window      type      pset
         0                                       0         0
*DATABASE_GLSTAT
$#      dt    binary      lcur     ioopt
       1.0         0         0         1
*DATABASE_MATSUM
$#      dt    binary      lcur     ioopt
       1.0         0         0         1
*DATABASE_TPRINT
$#      dt    binary      lcur     ioopt
       1.0         0         0         1
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas     nosol
       1.0         0       1.0       0.0     1e+08         0
*CASE
$#  caseid                                                                 jobid
         1trans                                                                 
$#                                                                      commands
                                                                                
$#   scid1     scid2     scid3     scid4     scid5     scid6     scid7     scid8
         1                                                                      
*CASE_BEGIN_1
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         1         0       2.0         0     1e+20       0.0
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         2         0      -1.0         0     1e+20       0.0
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         3         0       0.0         0     1e+20       0.0
*CASE_END_1
*CASE
$#  caseid                                                                 jobid
         2ab_l                                                                  
$#                                                                      commands
                                                                                
$#   scid1     scid2     scid3     scid4     scid5     scid6     scid7     scid8
         2                                                                      
*CASE_BEGIN_2
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         4         0       1.0         0     1e+20       0.0
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         8         0       0.0         0     1e+20       0.0
*CASE_END_2
*CASE
$#  caseid                                                                 jobid
         3ab_r                                                                  
$#                                                                      commands
                                                                                
$#   scid1     scid2     scid3     scid4     scid5     scid6     scid7     scid8
         3                                                                      
*CASE_BEGIN_3
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         6         0       1.0         0     1e+20       0.0
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         9         0       0.0         0     1e+20       0.0
*CASE_END_3
*CASE
$#  caseid                                                                 jobid
         4ot_l                                                                  
$#                                                                      commands
                                                                                
$#   scid1     scid2     scid3     scid4     scid5     scid6     scid7     scid8
         4                                                                      
*CASE_BEGIN_4
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         5         0       1.0         0     1e+20       0.0
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         8         0       0.0         0     1e+20       0.0
*CASE_END_4
*CASE
$#  caseid                                                                 jobid
         5ot_r                                                                  
$#                                                                      commands
                                                                                
$#   scid1     scid2     scid3     scid4     scid5     scid6     scid7     scid8
         5                                                                      
*CASE_BEGIN_5
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         7         0       1.0         0     1e+20       0.0
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         9         0       0.0         0     1e+20       0.0
*CASE_END_5
*CASE
$#  caseid                                                                 jobid
         6w_l                                                                   
$#                                                                      commands
                                                                                
$#   scid1     scid2     scid3     scid4     scid5     scid6     scid7     scid8
         6                                                                      
*CASE_BEGIN_6
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         4         0       1.0         0     1e+20       0.0
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         8         0       1.0         0     1e+20       0.0
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         5         0       0.0         0     1e+20       0.0
*CASE_END_6
*CASE
$#  caseid                                                                 jobid
         7w_r                                                                   
$#                                                                      commands
                                                                                
$#   scid1     scid2     scid3     scid4     scid5     scid6     scid7     scid8
         7                                                                      
*CASE_BEGIN_7
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         6         0       1.0         0     1e+20       0.0
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         9         0       1.0         0     1e+20       0.0
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         7         0       0.0         0     1e+20       0.0
*CASE_END_7
*CASE
$#  caseid                                                                 jobid
         8lr                                                                    
$#                                                                      commands
                                                                                
$#   scid1     scid2     scid3     scid4     scid5     scid6     scid7     scid8
         8                                                                      
*CASE_BEGIN_8
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         1         0       1.0         0     1e+20       0.0
*BOUNDARY_TEMPERATURE_SET
$#    nsid      lcid     cmult       loc    tdeath    tbirth
         2         0      -1.0         0     1e+20       0.0
*CASE_END_8
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
node_sets.k                                                                     
*END