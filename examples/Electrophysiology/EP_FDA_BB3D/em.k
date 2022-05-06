*KEYWORD
$--------1---------2---------3---------4---------5---------6---------7---------8
$                                                                              $
$                             EM CONTROL                                       $
$                                                                              $
$--------1---------2---------3---------4---------5---------6---------7---------8
*EM_CONTROL
$--------1---------2---------3---------4---------5---------6---------7---------8
$    emsol     numls  emdtinit   emdtmax   emtinit    emtend  ncyclFem  ncyclBem
        12         2       &dt                                       
*EM_CONTROL_TIMESTEP
$--------1---------2---------3---------4---------5---------6---------7---------8
$   TSTYPE   DTCONST      LCID    FACTOR     DTMIN     DTMAX
         1       &dt                    
*EM_SOLVER_FEM
$--------1---------2---------3---------4---------5---------6---------7---------8
$   relTol   maxIter solveType   precond   uselast
     1.e-6     10000         1         2         1
$--------1---------2---------3---------4---------5---------6---------7---------8
$                                                                              $
$                             EM MAT                                           $
$                                                                              $
$--------1---------2---------3---------4---------5---------6---------7---------8
$ ****************************************************************************
$ tissue
$ ****************************************************************************
*EM_MAT_005
$--------1---------2---------3---------4---------5---------6---------7---------8
$#intramid     mtype   sigmaXX   sigmaYY   sigmaZZ      beta        Cm
         1         2      &s11      &s22      &s33         3         2
$--------1---------2---------3---------4---------5---------6---------7---------8
$# sigmaXY   sigmaXZ   sigmaYX   sigmaYZ   sigmaZX   sigmaZY
                                                                   
$--------1---------2---------3---------4---------5---------6---------7---------8
$#extra               sigmaXXb  sigmaYYb  sigmaZZb
                         &s11b     &s22b     &s33b
$--------1---------2---------3---------4---------5---------6---------7---------8
$#sigmaXYb  sigmaXZb  sigmaYXb  sigmaYZb  sigmaZXb  sigmaZYb
               
$--------1---------2---------3---------4---------5---------6---------7---------8
$#    AOPT        XP        YP        ZP        A1        A2        A3      MACF
        2.                             0         1         0         0
$--------1---------2---------3---------4---------5---------6---------7---------8
$#      V1        V2        V3        D1        D2        D3
                                       0         1         0
$ ****************************************************************************
$ bath
$ ****************************************************************************
*EM_MAT_005
$--------1---------2---------3---------4---------5---------6---------7---------8
$#intramid     mtype   sigmaXX   sigmaYY   sigmaZZ      beta        Cm
         2         4        0.        0.        0.         3         2
$--------1---------2---------3---------4---------5---------6---------7---------8
$# sigmaXY   sigmaXZ   sigmaYX   sigmaYZ   sigmaZX   sigmaZY
                                                                   
$--------1---------2---------3---------4---------5---------6---------7---------8
$#extra               sigmaXXb  sigmaYYb  sigmaZZb
                       &s_bath   &s_bath   &s_bath
$--------1---------2---------3---------4---------5---------6---------7---------8
$#sigmaXYb  sigmaXZb  sigmaYXb  sigmaYZb  sigmaZXb  sigmaZYb
               
$--------1---------2---------3---------4---------5---------6---------7---------8
$#    AOPT        XP        YP        ZP        A1        A2        A3      MACF
        2.                             0         1         0         0
$--------1---------2---------3---------4---------5---------6---------7---------8
$#      V1        V2        V3        D1        D2        D3
                                       0         1         0
$--------1---------2---------3---------4---------5---------6---------7---------8
$                                                                              $
$                             EM CELLMODEL DEFINEFUCTION                       $
$                                                                              $
$--------1---------2---------3---------4---------5---------6---------7---------8
*EM_EP_CELLMODEL_DEFINEFUNCTION
$--------1---------2---------3---------4---------5---------6---------7---------8
$    matid  numState funcVsDer potVsIion                                  iwrite
         1         3         1         0                                       2
$ define functions for d/dt as a function of (V,u1,u2,u3,...)
$     dVdt     du1dt     du2dt     du3dt     du4dt     du5dt     du6dt     du7dt
       200       201       202       203     
$ define functions for u(t=0,x) as a function of (x,y,z)
$      V_0      u1_0      u2_0      u3_0      u4_0      u5_0      u6_0      u7_0
       300       301       302       303
$ ****************************************************************************
*EM_EP_CELLMODEL_DEFINEFUNCTION
$--------1---------2---------3---------4---------5---------6---------7---------8
$    matid  numState funcVsDer potVsIion                                  iwrite
         2         3         1         0                                       2
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
            float time, float timeStep,
            float V, float u1, float u2, float u3, float u4,
            float u5, float u6, float u7)
{ 
  ksi=3;
  Cm=2;
  float beta;
  beta=-0.32218254069;
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
            float time, float timeStep,
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
            float time, float timeStep,
            float V, float u1, float u2, float u3, float u4,
            float u5, float u6, float u7)
{ 
  ksi=3;
  Cm=2;
  float ddt;
  ddt=-(u1+u3-V)*u2**3.;
$ printf("V_value: %f",ddt);
  return ddt;
}
$ ****************************************************************************
*DEFINE_FUNCTION
203
float du3dt(
            float ksi, float Cm,
            float time, float timeStep,
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
  float alpha;
  float Se;
  float k;
  float sigma_e11;
  float sigma11;
  pi=3.14159265359;
  alpha=0.01;
  k=1/sqrt(2);
  sigma11=0.1114533;
  sigma_e11=(1-k)*sigma11/k;
  Se=sigma_e11;
  F=cos(pi*x);
  init=F-alpha/Se*x;
$ printf("fonctionF %f",F);	
  return init;
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
  F=cos(pi*x);
  G=1.+x*y**2.*z**3.;
  init=G+F;
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
  F=cos(pi*x);
  G=1.+x*y**2.*z**3.;
  init=1./sqrt(G);
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
  float alpha;
  float Se;
  float k;
  float sigma_e11;
  float sigma11;
  pi=3.14159265359;
  alpha=0.01;
  k=1/sqrt(2);
  sigma11=0.1114533;
  sigma_e11=(1-k)*sigma11/k;
  Se=sigma_e11;
  F=cos(pi*x);
  G=1.+x*y**2.*z**3.;
  init=-alpha/Se*x;
  return init;
}
$--------1---------2---------3---------4---------5---------6---------7---------8
$                                                                              $
$                             EM BCs                                           $
$                                                                              $
$--------1---------2---------3---------4---------5---------6---------7---------8
*EM_BOUNDARY_PRESCRIBED
$--------1---------2---------3---------4---------5---------6---------7---------8
$     bpid    bptype bpsettype    xsetid     value     xlcid    xprint
         1         3         2         1        0.
*EM_BOUNDARY_PRESCRIBED
$--------1---------2---------3---------4---------5---------6---------7---------8
$     bpid    bptype bpsettype    xsetid     value     xlcid    xprint
         2         4         2         2 +2.267e-5
$--------1---------2---------3---------4---------5---------6---------7---------8
$                                                                              $
$                             EM OUTPUT                                        $
$                                                                              $
$--------1---------2---------3---------4---------5---------6---------7---------8
*EM_OUTPUT
$--------1---------2---------3---------4---------5---------6---------7---------8
$     matS      matF      solS      solF      mesh    memory    timing d3plotAsc
         4         4         4         4         0         0         0         0
$      mf2       gmv d3plotFor  timeHist spoolesTr    intDat
         0         0         0         0         0         0
*END
