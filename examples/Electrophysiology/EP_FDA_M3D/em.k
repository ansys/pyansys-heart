*KEYWORD
$---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
$                                                                              $
$                             EM CONTROL                                       $
$                                                                              $
$---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
*EM_CONTROL
$--------1---------2---------3---------4---------5---------6---------7---------8
$    emsol     numls  emdtinit   emdtmax   emtinit    emtend  ncyclFem  ncyclBem
        11         1       &dt
*EM_CONTROL_TIMESTEP
$--------1---------2---------3---------4---------5---------6---------7---------8
$   TSTYPE   DTCONST      LCID    FACTOR     DTMIN     DTMAX
         1       &dt                    
*EM_SOLVER_FEM
$--------1---------2---------3---------4---------5---------6---------7---------8
$   relTol   maxIter solveType   precond   uselast
     1.e-6     10000         1         2         1
$---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
$                                                                              $
$                             EM MAT                                           $
$                                                                              $
$---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
*EM_MAT_003
$--------1---------2---------3---------4---------5---------6---------7---------8
$    matid     mtype   sigmaXX   sigmaYY   sigmaZZ      beta        Cm
         1         2      &s11      &s22      &s33         3         2
$--------1---------2---------3---------4---------5---------6---------7---------8
$  sigmaXY   sigmaXZ   sigmaYX   sigmaYZ   sigmaZX   sigmaZY      AOPT
                                                                    2.
$--------1---------2---------3---------4---------5---------6---------7---------8
$       XP        YP        ZP        A1        A2        A3      MACF
                                       1         0         0
$--------1---------2---------3---------4---------5---------6---------7---------8
$       V1        V2        V3        D1        D2        D3
                                       0         1         0
$---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
$                                                                              $
$                             EM CELLMODEL DEFINEFUCTION                       $
$                                                                              $
$---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
*EM_EP_CELLMODEL_DEFINEFUNCTION
$--------1---------2---------3---------4---------5---------6---------7---------8
$    matid  numState funcVsDer potVslion                                  iwrite
         1         3         1         0                                       1
$ define functions for d/dt as a function of (V,u1,u2,u3,...)
$     dVdt     du1dt     du2dt     du3dt     du4dt     du5dt     du6dt     du7dt
       200       201       202       203     
$ define functions for u(t=0,x) as a function of (x,y,z)
$      V_0      u1_0      u2_0      u3_0      u4_0      u5_0      u6_0      u7_0
       300       301       302       303
*EM_EP_CELLMODEL_DEFINEFUNCTION
$--------1---------2---------3---------4---------5---------6---------7---------8
$    matid  numState
         2         3
$ define functions for d/dt as a function of (V,u1,u2,u3,...)
$     dVdt     du1dt     du2dt     du3dt     du4dt     du5dt     du6dt     du7dt
       200       201       202       203     
$ define functions for u(t=0,x) as a function of (x,y,z)
$      V_0      u1_0      u2_0      u3_0      u4_0      u5_0      u6_0      u7_0
       300       301       302       303
$ ****************************************************************************
*DEFINE_FUNCTION
200
float dVdt(
            float ksi, float Cm,
            float V, float u1, float u2, float u3, float u4,
            float u5, float u6, float u7)
{ 
  ksi=3;
  Cm=2;	
  float beta;
$  beta=1.;
  beta=-8.6;
  float ddt;
  ddt=-Cm/2.*(u1+u3-V)*u2**2.*(V-u3)+beta*(V-u3)/ksi;
  ddt=-ddt/Cm;
  return ddt;
}
$ ****************************************************************************
*DEFINE_FUNCTION
201
float du1dt(
            float ksi, float Cm,
            float V, float u1, float u2, float u3, float u4,
            float u5, float u6, float u7)
{
  ksi=3;
  Cm=2;	
  float ddt;
  ddt=(u1+u3-V)**2.*u2**2.;
  ddt=ddt+0.5*(u1+u3-V)*u2**2.*(V-u3);
  return ddt;
}
$ ****************************************************************************
*DEFINE_FUNCTION
202
float du2dt(
            float ksi, float Cm,
            float V, float u1, float u2, float u3, float u4,
            float u5, float u6, float u7)
{ 
  ksi=3;
  Cm=2;	
  float ddt;
  ddt=-(u1+u3-V)*u2**3.;
  return ddt;
}
$ ****************************************************************************
*DEFINE_FUNCTION
203
float du3dt(
            float ksi, float Cm,
            float V, float u1, float u2, float u3, float u4,
            float u5, float u6, float u7)
{
  ksi=3;
  Cm=2;	
  float ddt;
  ddt=0.;
  return ddt;
}
$ ****************************************************************************
*DEFINE_FUNCTION
300
float V_0(
            float time, float x, float y, float z)
{ 
  float init;
  float F;
  float pi;
  pi=3.14159265359;
$  F=pi*z;
  F=cos(pi*x)*cos(2.*pi*y)*cos(3.*pi*z);
  init=sqrt(1.+time)*F;
$ printf("fonctionF %f",F);	
  return F;
}
$ ****************************************************************************
*DEFINE_FUNCTION
301
float u1_0(
            float time, float x, float y, float z)
{ 
  float init;
  float F;
  float G;
  float pi;
  pi=3.14159265359;
  F=cos(pi*x)*cos(2.*pi*y)*cos(3.*pi*z);
  G=1.+x*y**2.*z**3.;
  init=(1.+time)*G+sqrt(1.+time)*F;
  return init;
}
$ ****************************************************************************
*DEFINE_FUNCTION
302
float u2_0(
            float time, float x, float y, float z)
{ 
  float init;
  float F;
  float G;
  float pi;
  pi=3.14159265359;
  F=cos(pi*x)*cos(2.*pi*y)*cos(3.*pi*z);
  G=1.+x*y**2.*z**3.;
  init=1./(1.+time)/sqrt(G);
  return init;
}
$ ****************************************************************************
*DEFINE_FUNCTION
303
float u3_0(
            float time, float x, float y, float z)
{ 
  float init;
  float F;
  float G;
  float pi;
  pi=3.14159265359;
  F=cos(pi*x)*cos(2.*pi*y)*cos(3.*pi*z);
  G=1.+x*y**2.*z**3.;
  init=0.;
  return init;
}
$---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
$                                                                              $
$                             EM OUTPUT                                        $
$                                                                              $
$---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
*EM_OUTPUT
$--------1---------2---------3---------4---------5---------6---------7---------8
$     matS      matF      solS      solF      mesh    memory    timing d3plotAsc
         4         4         4         4         0         0         0         0
$      mf2       gmv d3plotFor  timeHist spoolesTr    intDat
         0         0         0         0         0         0
*END
