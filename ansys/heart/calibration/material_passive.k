*KEYWORD
*CONTROL_MPP_IO_NODUMP
*PARAMETER
$ LSOPT inputs
rp1,<<ps1>>
rp2,<<ps2>>
$$$ra,1.05
$$$rb,7.542
$$$raf,3.465
$$$rbf,14.472
ra,2.36e-3
rb,1.75
raf,0.49e-3
rbf,9.01
*PARAMETER_EXPRESSION
ra_ps,a*p1
rb_ps,b*p2
raf_ps,af*p1
rbf_ps,bf*p2
*MAT_295
$#     mid       rho      aopt
         1     1e-03       2.0
$#   title     itype      beta        nu
ISO               -3       0.0     0.499
$#      k1        k2
     &a_ps     &b_ps
$#   title     atype    intype        nf
ANISO             -1         0         1
$#   theta         a         b
       0.0      0.08      0.76
$#   ftype      fcid        k1        k2
         1         0    &af_ps    &bf_ps
$#      xp        yp        zp        a1        a2        a3      macf
                                                                     1
$#      v1        v2        v3        d1        d2        d3      beta       ref
                                                                                
*MAT_295
$#     mid       rho      aopt
         2     1e-03       2.0
$#   title     itype      beta        nu
ISO               -3       0.0     0.499
$#      k1        k2
     &a_ps     &b_ps
$#   title     atype    intype        nf
ANISO             -1         0         1
$#   theta         a         b
       0.0      0.08      0.76
$#   ftype      fcid        k1        k2
         1         0    &af_ps    &bf_ps
$#      xp        yp        zp        a1        a2        a3      macf
                                                                     1
$#      v1        v2        v3        d1        d2        d3      beta       ref
                                                                                
*MAT_295
$#     mid       rho      aopt
         3     1e-03       2.0
$#   title     itype      beta        nu
ISO               -3       0.0     0.499
$#      k1        k2
     &a_ps     &b_ps
$#   title     atype    intype        nf
ANISO             -1         0         1
$#   theta         a         b
       0.0      0.08      0.76
$#   ftype      fcid        k1        k2
         1         0    &af_ps    &bf_ps
$#      xp        yp        zp        a1        a2        a3      macf
                                                                     1
$#      v1        v2        v3        d1        d2        d3      beta       ref
                                                                                
*END