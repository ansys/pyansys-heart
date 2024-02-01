*KEYWORD
*DEFINE_CONTROL_VOLUME
$#      id       sid        p0
         1         1          
*DEFINE_CONTROL_VOLUME
$#      id       sid        p0
         2         2          
*DEFINE_CONTROL_VOLUME_INTERACTION
$#      id     cvid1     cvid2     lcid_     area_
         1         1         0        10          
*DEFINE_CONTROL_VOLUME_INTERACTION
$#      id     cvid1     cvid2     lcid_     area_
         2         2         0        11          
*DEFINE_FUNCTION
$#     fid                                                               heading
        10constant_preload_windkessel_afterload_left                            
$#                                                                      function
float constant_preload_windkessel_afterload_left(float t, float dp)
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
    Rv = 9.311235537063599e-07;
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
    qven = chi_mv * ( ( pven - pk ) / Rv );
    qart = chi_av * ( ( pk  - part) / Ra );
    qk  = qven - qart;
$   used to debug
$   Note: we write at every call of t, write states for time t
    char fn_bug[] = "circulation_model_debug.csv";
    FILE *f_bug;
    if (icall == 0){
        f_bug=fopen(fn_bug, "w");
        fprintf(f_bug, "icall,is_new_dt,time,dp,pk,part,pven");
        fprintf(f_bug, ",qart,qp,qven,qk,vart\n");      
        fprintf(f_bug, "%d,%d,%.6e,%.6e,",icall,is_new_dt,t,dp);
        fprintf(f_bug, "%.6e,%.6e,%.6e,",pk,part,pven);
        fprintf(f_bug, "%.6e,%.6e,%.6e,%.6e,",qart,qp,qven,qk);
        fprintf(f_bug, "%.6e\n",vart);
      }
    else {
        f_bug=fopen(fn_bug, "a");
        fprintf(f_bug, "%d,%d,%.6e,%.6e,",icall,is_new_dt,t,dp);
        fprintf(f_bug, "%.6e,%.6e,%.6e,",pk,part,pven);
        fprintf(f_bug, "%.6e,%.6e,%.6e,%.6e,",qart,qp,qven,qk);
        fprintf(f_bug, "%.6e\n",vart);
    }
    fclose(f_bug);
$   write data to file
$   Note: we write at the first call of t, write the states for time t_old
    char fn_data[] = "constant_preload_windkessel_afterload_left.csv";
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
*DEFINE_FUNCTION
$#     fid                                                               heading
        11constant_preload_windkessel_afterload_right                           
$#                                                                      function
float constant_preload_windkessel_afterload_right(float t, float dp)
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
    Rv = 4.6556177685317995e-07;
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
    qven = chi_mv * ( ( pven - pk ) / Rv );
    qart = chi_av * ( ( pk  - part) / Ra );
    qk  = qven - qart;
$   used to debug
$   Note: we write at every call of t, write states for time t
    char fn_bug[] = "circulation_model_debug.csv";
    FILE *f_bug;
    if (icall == 0){
        f_bug=fopen(fn_bug, "w");
        fprintf(f_bug, "icall,is_new_dt,time,dp,pk,part,pven");
        fprintf(f_bug, ",qart,qp,qven,qk,vart\n");      
        fprintf(f_bug, "%d,%d,%.6e,%.6e,",icall,is_new_dt,t,dp);
        fprintf(f_bug, "%.6e,%.6e,%.6e,",pk,part,pven);
        fprintf(f_bug, "%.6e,%.6e,%.6e,%.6e,",qart,qp,qven,qk);
        fprintf(f_bug, "%.6e\n",vart);
      }
    else {
        f_bug=fopen(fn_bug, "a");
        fprintf(f_bug, "%d,%d,%.6e,%.6e,",icall,is_new_dt,t,dp);
        fprintf(f_bug, "%.6e,%.6e,%.6e,",pk,part,pven);
        fprintf(f_bug, "%.6e,%.6e,%.6e,%.6e,",qart,qp,qven,qk);
        fprintf(f_bug, "%.6e\n",vart);
    }
    fclose(f_bug);
$   write data to file
$   Note: we write at the first call of t, write the states for time t_old
    char fn_data[] = "constant_preload_windkessel_afterload_right.csv";
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
*END