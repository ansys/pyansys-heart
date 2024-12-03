*KEYWORD
*MAT_ELASTIC
$#     mid        ro         e        pr        da        db    unused
         1     1e-06       1.0                                        
*MAT_ELASTIC
$#     mid        ro         e        pr        da        db    unused
         2     1e-06       1.0                                        
*MAT_ELASTIC
$#     mid        ro         e        pr        da        db    unused
         3     1e-06       1.0                                        
$$ Left ventricle $$
*EM_MAT_003
$#     mid     mtype   sigma11   sigma22   sigma33      beta        cm
         1         2       0.5       0.1       0.1     140.0      0.01
$# sigma12   sigma13   sigma21   sigma23   sigma31   sigma32      aopt   lambda_
                                                                     2          
$#      xp        yp        zp        a1        a2        a3      macf
                                     0.0       0.0       1.0         1
$#      v1        v2        v3        d1        d2        d3
                                     0.0      -1.0       0.0
$$ Right ventricle $$
*EM_MAT_003
$#     mid     mtype   sigma11   sigma22   sigma33      beta        cm
         2         2       0.5       0.1       0.1     140.0      0.01
$# sigma12   sigma13   sigma21   sigma23   sigma31   sigma32      aopt   lambda_
                                                                     2          
$#      xp        yp        zp        a1        a2        a3      macf
                                     0.0       0.0       1.0         1
$#      v1        v2        v3        d1        d2        d3
                                     0.0      -1.0       0.0
$$ Septum $$
*EM_MAT_003
$#     mid     mtype   sigma11   sigma22   sigma33      beta        cm
         3         2       0.5       0.1       0.1     140.0      0.01
$# sigma12   sigma13   sigma21   sigma23   sigma31   sigma32      aopt   lambda_
                                                                     2          
$#      xp        yp        zp        a1        a2        a3      macf
                                     0.0       0.0       1.0         1
$#      v1        v2        v3        d1        d2        d3
                                     0.0      -1.0       0.0
*END