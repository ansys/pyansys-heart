*KEYWORD
*CONTROL_TIME_STEP
$#  dtinit    tssfac      isdo    tslimt     dt2ms      lctm     erode     ms1st
       1.0       0.0         0       0.0       1.0         0         0         0
$#  dt2msf   dt2mslc     imscl    unused    unused     rmscl     emscl      ihdo
                                                                                
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas     nosol
      10.0         0       0.0       0.0     1e+08         0
*DATABASE_BINARY_D3PLOT
$#      dt      lcdt      beam     npltc    psetid
       1.0                   0                    
$#   ioopt      rate    cutoff    window      type      pset
         0                                       0         0
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
ep_settings.k                                                                   
*INCLUDE
$#                                                                      filename
create_fiber.k                                                                  
*END