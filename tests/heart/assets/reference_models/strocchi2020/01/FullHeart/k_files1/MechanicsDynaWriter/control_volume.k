*KEYWORD
*SECTION_SHELL
$#   secid    elform      shrf       nip     propt   qr/irid     icomp     setyp
         8         4    0.8333         3       1.0         0         0         1
$#      t1        t2        t3        t4      nloc     marea      idof    edgset
       1.0       0.0       0.0       0.0       0.0       0.0       0.0         0
*MAT_NULL
$#     mid        ro        pc        mu     terod     cerod        ym        pr
        14     0.001                                                            
*PART
$#                                                             heading
null flow area                                                        
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
        26         8        14                                                  
*DEFINE_CONTROL_VOLUME
$#      id       sid        p0
         1         1          
*SET_SEGMENT_ADD
$#     sid
        29
$#    sets      sets
        16        17
*DEFINE_CONTROL_VOLUME_FLOW_AREA
$#    FAID     FCIID     FASID   FASTYPE       PID
         1         1        29         2        26
*DEFINE_CONTROL_VOLUME_INTERACTION
$#      id     cvid1     cvid2     lcid      area 
         1         1         0         5          
*DEFINE_FUNCTION
$#     fid                                                               heading
         5afterload_windkessel_left                                             
$#                                                                      function
float afterload_windkessel_left(float t, float dp)
{
$   numerical constant
    int Implicit=1;
$   Only used for Euler Implicit
    float gamma = 0.6;
$   physical constants
    float Rp, Ca;
    float Ra, Rv;
    Rp = 0.00015518725895105997;
    Ca = 12887655.942365238;
    Ra = 4.6556177685317995e-06;
    Rv = 4.6556177685317995e-06;
$   physical variables
    float chi_av, chi_mv;
    float pk, part, pven;
    float qk, qven, qart, qp;
    float vart;
$   only for save data purpose
    float pk2, part2, pven2;
    float qk2, qven2, qart2, qp2;
    float vart2;
$   constant pre load:
    pven = 0.0019998358112249997;
$   time related variables
    int icall=0, is_new_dt=0;
$   t: current time (input)
$   t_last: time of the last call
$   t_old: time of the last step
    float t_last=0.0, t_old=0.0, dt;
$   initialisation at t=0
    if (icall == 0) {
          part = 0.009332567119049997;
$   initial arterial volume
          vart = Ca * part;
          qp = part / Rp;
    }
$   determine if function call is at start of timestep:
    if ( t-t_last > 1e-9 ) {
        is_new_dt = 1;
        dt = t - t_old;
    }
    else if ( t-t_last < 0. ) {
        printf("  ## Warning bisection may not be properly handled ##");
        printf("  ## Warning: dt_old: %f", dt );
        is_new_dt = 0;
        dt = dt - (t_last-t);
        printf("## Warning: dt_new: %f", dt );
$       abort(0);
    } else
    {
        is_new_dt = 0;
    }
    if ( is_new_dt ) {
$   Save system states of last time step (at t_old)
$   The converged pressure value of last time step (at t_old)
        pk2 = dp;
        part2 = part;
        pven2 = pven;
        vart2 = vart;
        qp2 = qp;
        qart2 = qart;
        qven2 = qven;
        qk2 = qk;
$   Update system states for new time step (at t)
        vart = vart + dt * (qart-qp);
        part = vart / Ca;
        qp = part / Rp;
    }
    if (Implicit){
$   LSDYNA will integrate cavity volume implicitly: V^t = V^t_old+dt*Q^t
$   LSDYNAs input dp is interpolated by dp=(1-r)*p^t_old+r*p^t+1_i
$   This is not suitable to check the valve opening (to compute Q at t)
$   We retrieve firstly p^t at this iteration
        pk = (dp -(1-gamma)*pk2)/gamma;
    } else
    {
$   LSDYNA will integrate cavity volume explicitly: V^t = V^t_old+dt*Q^t_old
        pk = pk2;
    }
    
$   Update valve indicator functions
    if (pven >= pk )
    {
        chi_mv = 1;
    } else
    {
        chi_mv = 1.e-16;
    }
    if ( pk >= part )
    {
        chi_av = 1;
    } else {
        chi_av = 1.e-16;
    }
$   compute flow: In - Out
    qven = 0.0;
    qart = chi_av * ( ( pk  - part) / Ra );
    qk  = qven - qart;
$   write data to file
$   Note: we write at the first call of t, write the states for time t_old
    char fn_data[] = "afterload_windkessel_left.csv";
    FILE *f_data;
    if (icall == 0){
        f_data=fopen(fn_data, "w");
        fprintf(f_data, "icall,time,pk,part,pven");
        fprintf(f_data, ",qart,qp,qven,qk,vart\n"); 
        fclose(f_data);
     }
    else if ( is_new_dt ) {
        f_data=fopen(fn_data, "a");
        fprintf(f_data, "%d,%.6e,",icall,t_old);
        fprintf(f_data, "%.6e,%.6e,%.6e,",pk2,part2,pven2);
        fprintf(f_data, "%.6e,%.6e,%.6e,%.6e,",qart2,qp2,qven2,qk2);
        fprintf(f_data, "%.6e\n",vart2);
        fclose(f_data);
$       
        t_old = t;
    }
    
$   Update counters
    t_last = t;
    icall = icall + 1;
    
$   LSDYNA defines outflow as positive
    return -qk;
}
*DEFINE_CONTROL_VOLUME
$#      id       sid        p0
         2         2          
*SET_SEGMENT_ADD
$#     sid
        30
$#    sets      sets
        18        19
*DEFINE_CONTROL_VOLUME_FLOW_AREA
$#    FAID     FCIID     FASID   FASTYPE       PID
         2         2        30         2        26
*DEFINE_CONTROL_VOLUME_INTERACTION
$#      id     cvid1     cvid2     lcid      area 
         2         2         0         6          
*DEFINE_FUNCTION
$#     fid                                                               heading
         6afterload_windkessel_right                                            
$#                                                                      function
float afterload_windkessel_right(float t, float dp)
{
$   numerical constant
    int Implicit=1;
$   Only used for Euler Implicit
    float gamma = 0.6;
$   physical constants
    float Rp, Ca;
    float Ra, Rv;
    Rp = 1.9398407368882496e-05;
    Ca = 57994451.74064357;
    Ra = 1.6294662189861296e-06;
    Rv = 2.3278088842658997e-06;
$   physical variables
    float chi_av, chi_mv;
    float pk, part, pven;
    float qk, qven, qart, qp;
    float vart;
$   only for save data purpose
    float pk2, part2, pven2;
    float qk2, qven2, qart2, qp2;
    float vart2;
$   constant pre load:
    pven = 0.0010665790993199997;
$   time related variables
    int icall=0, is_new_dt=0;
$   t: current time (input)
$   t_last: time of the last call
$   t_old: time of the last step
    float t_last=0.0, t_old=0.0, dt;
$   initialisation at t=0
    if (icall == 0) {
          part = 0.0019998358112249997;
$   initial arterial volume
          vart = Ca * part;
          qp = part / Rp;
    }
$   determine if function call is at start of timestep:
    if ( t-t_last > 1e-9 ) {
        is_new_dt = 1;
        dt = t - t_old;
    }
    else if ( t-t_last < 0. ) {
        printf("  ## Warning bisection may not be properly handled ##");
        printf("  ## Warning: dt_old: %f", dt );
        is_new_dt = 0;
        dt = dt - (t_last-t);
        printf("## Warning: dt_new: %f", dt );
$       abort(0);
    } else
    {
        is_new_dt = 0;
    }
    if ( is_new_dt ) {
$   Save system states of last time step (at t_old)
$   The converged pressure value of last time step (at t_old)
        pk2 = dp;
        part2 = part;
        pven2 = pven;
        vart2 = vart;
        qp2 = qp;
        qart2 = qart;
        qven2 = qven;
        qk2 = qk;
$   Update system states for new time step (at t)
        vart = vart + dt * (qart-qp);
        part = vart / Ca;
        qp = part / Rp;
    }
    if (Implicit){
$   LSDYNA will integrate cavity volume implicitly: V^t = V^t_old+dt*Q^t
$   LSDYNAs input dp is interpolated by dp=(1-r)*p^t_old+r*p^t+1_i
$   This is not suitable to check the valve opening (to compute Q at t)
$   We retrieve firstly p^t at this iteration
        pk = (dp -(1-gamma)*pk2)/gamma;
    } else
    {
$   LSDYNA will integrate cavity volume explicitly: V^t = V^t_old+dt*Q^t_old
        pk = pk2;
    }
    
$   Update valve indicator functions
    if (pven >= pk )
    {
        chi_mv = 1;
    } else
    {
        chi_mv = 1.e-16;
    }
    if ( pk >= part )
    {
        chi_av = 1;
    } else {
        chi_av = 1.e-16;
    }
$   compute flow: In - Out
    qven = 0.0;
    qart = chi_av * ( ( pk  - part) / Ra );
    qk  = qven - qart;
$   write data to file
$   Note: we write at the first call of t, write the states for time t_old
    char fn_data[] = "afterload_windkessel_right.csv";
    FILE *f_data;
    if (icall == 0){
        f_data=fopen(fn_data, "w");
        fprintf(f_data, "icall,time,pk,part,pven");
        fprintf(f_data, ",qart,qp,qven,qk,vart\n"); 
        fclose(f_data);
     }
    else if ( is_new_dt ) {
        f_data=fopen(fn_data, "a");
        fprintf(f_data, "%d,%.6e,",icall,t_old);
        fprintf(f_data, "%.6e,%.6e,%.6e,",pk2,part2,pven2);
        fprintf(f_data, "%.6e,%.6e,%.6e,%.6e,",qart2,qp2,qven2,qk2);
        fprintf(f_data, "%.6e\n",vart2);
        fclose(f_data);
$       
        t_old = t;
    }
    
$   Update counters
    t_last = t;
    icall = icall + 1;
    
$   LSDYNA defines outflow as positive
    return -qk;
}
*DEFINE_CONTROL_VOLUME
$#      id       sid        p0
         3         3          
*SET_SEGMENT_ADD
$#     sid
        31
$#    sets      sets      sets      sets      sets      sets
        20        21        22        23        24        25
*DEFINE_CONTROL_VOLUME_FLOW_AREA
$#    FAID     FCIID     FASID   FASTYPE       PID
         3         3        31         2        26
*DEFINE_CONTROL_VOLUME_INTERACTION
$#      id     cvid1     cvid2     lcid      area 
         3         3         0         7          
*DEFINE_FUNCTION
$#     fid                                                               heading
         7constant_flow_left_atrium                                             
$#                                                                      function
float constant_flow_left_atrium(float t, float dp, float area) 
{
return -83.0;
}
*DEFINE_CONTROL_VOLUME_INTERACTION
$#      id     cvid1     cvid2     lcid      area 
         4         3         1         8          
*DEFINE_FUNCTION
$#     fid                                                               heading
         8valve_mitral                                                          
$#                                                                      function
float valve_mitral(float t, float dp)
{
   float Rv, qven;
   Rv = 1e-06;
   if (dp >= 0.0 )
   {
       qven =  dp  / Rv ;
   }
   else{
       qven=0.0;
   }
   return qven;
}
*DEFINE_CONTROL_VOLUME
$#      id       sid        p0
         4         4          
*SET_SEGMENT_ADD
$#     sid
        32
$#    sets      sets      sets
        26        27        28
*DEFINE_CONTROL_VOLUME_FLOW_AREA
$#    FAID     FCIID     FASID   FASTYPE       PID
         4         5        32         2        26
*DEFINE_CONTROL_VOLUME_INTERACTION
$#      id     cvid1     cvid2     lcid      area 
         5         4         0         9          
*DEFINE_FUNCTION
$#     fid                                                               heading
         9constant_flow_right_atrium                                            
$#                                                                      function
float constant_flow_right_atrium(float t, float dp, float area) 
{
return -83.0;
}
*DEFINE_CONTROL_VOLUME_INTERACTION
$#      id     cvid1     cvid2     lcid      area 
         6         4         2        10          
*DEFINE_FUNCTION
$#     fid                                                               heading
        10valve_tricuspid                                                       
$#                                                                      function
float valve_tricuspid(float t, float dp)
{
   float Rv, qven;
   Rv = 1e-06;
   if (dp >= 0.0 )
   {
       qven =  dp  / Rv ;
   }
   else{
       qven=0.0;
   }
   return qven;
}
*END