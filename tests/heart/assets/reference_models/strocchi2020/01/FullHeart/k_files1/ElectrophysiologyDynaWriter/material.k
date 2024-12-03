*KEYWORD
*MAT_ELASTIC
$#     mid        ro         e        pr        da        db    unused
         1     1e-06       1.0                                        
*MAT_ELASTIC
$#     mid        ro         e        pr        da        db    unused
         2     1e-06       1.0                                        
*MAT_ELASTIC
$#     mid        ro         e        pr        da        db    unused
         7     1e-06       1.0                                        
*MAT_ELASTIC
$#     mid        ro         e        pr        da        db    unused
         3     1e-06       1.0                                        
*MAT_ELASTIC
$#     mid        ro         e        pr        da        db    unused
         4     1e-06       1.0                                        
*MAT_ELASTIC
$#     mid        ro         e        pr        da        db    unused
         5     1e-06       1.0                                        
*MAT_ELASTIC
$#     mid        ro         e        pr        da        db    unused
         6     1e-06       1.0                                        
*MAT_ELASTIC
$#     mid        ro         e        pr        da        db    unused
         8     1e-06       1.0                                        
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
         7         2       0.5       0.1       0.1     140.0      0.01
$# sigma12   sigma13   sigma21   sigma23   sigma31   sigma32      aopt   lambda_
                                                                     2          
$#      xp        yp        zp        a1        a2        a3      macf
                                     0.0       0.0       1.0         1
$#      v1        v2        v3        d1        d2        d3
                                     0.0      -1.0       0.0
$$ Left atrium $$
*EM_MAT_001
$#     mid     mtype     sigma     eosid    unused   lambda_    deatht   rdltype
         3         4       0.5                                   1e+28         0
$#   empty      freq     epsrr    eosid2     epsri    eosid3      eps0
                                                                      
$#    beta        cm
     140.0      0.01
$$ Right atrium $$
*EM_MAT_001
$#     mid     mtype     sigma     eosid    unused   lambda_    deatht   rdltype
         4         4       0.5                                   1e+28         0
$#   empty      freq     epsrr    eosid2     epsri    eosid3      eps0
                                                                      
$#    beta        cm
     140.0      0.01
$$ Aorta $$
*EM_MAT_001
$#     mid     mtype     sigma     eosid    unused   lambda_    deatht   rdltype
         5         4       0.5                                   1e+28         0
$#   empty      freq     epsrr    eosid2     epsri    eosid3      eps0
                                                                      
$#    beta        cm
     140.0      0.01
$$ Pulmonary artery $$
*EM_MAT_001
$#     mid     mtype     sigma     eosid    unused   lambda_    deatht   rdltype
         6         4       0.5                                   1e+28         0
$#   empty      freq     epsrr    eosid2     epsri    eosid3      eps0
                                                                      
$#    beta        cm
     140.0      0.01
$$ Atrioventricular isolation $$
*EM_MAT_001
$#     mid     mtype     sigma     eosid    unused   lambda_    deatht   rdltype
         8         1       0.0                                   1e+28         0
$#   empty      freq     epsrr    eosid2     epsri    eosid3      eps0
                                                                      
$#    beta        cm
       0.0       0.0
*END