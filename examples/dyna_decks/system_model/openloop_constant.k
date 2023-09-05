*DEFINE_FUNCTION
16
float constant_pre_after_load(float t, float dp)
{
$   numerical constant
    int Implicit=1;
$   Only used for Euler Implicit
    float gamma = 0.6;

$   physical constants
    float Ra, Rv;
    Ra = (13.332 / 83333)/10;
    Rv = (0.4 / 83333);

$   constant pre/afterload:
    float pven, part;
    pven = 2;
    part = 8;

$   physical variables
    float chi_av, chi_mv;
    float pk, qk, qven, qart;

$   only for save data purpose
    float pk2, qk2, qven2, qart2;

$   time related variables
    int icall=0, is_new_dt=0;
$   t: current time (input)
$   t_last: time of the last call
$   t_old: time of the last step
    float t_last=0.0, t_old=0.0, dt;

$   determine if function call is at start of timestep:
    if ( t-t_last > 1e-9 ) {
        is_new_dt = 1;
        dt = t - t_old;
    }
    else if ( t-t_last < 0. ) {
        printf("## bisection is not handled now ##");
        is_new_dt = 0;
    } else
    {
        is_new_dt = 0;
    }
    if ( is_new_dt ) {
$   Save system states of last time step (at t_old)
$   The converged pressure value of last time step (at t_old)
        pk2 = dp;

        qart2 = qart;
        qven2 = qven;
        qk2 = qk;
    }

    if (Implicit){
$   LSDYNA will integrate cavity volume implicitly: V^t = V^t_old+dt*Q^t
$   LSDYNA's input dp is interpolated by dp=(1-r)*p^t_old+r*p^t+1_i
$   This is not suitable to check the valve's opening (to compute Q at t)
$   We retreive firstly p^t at this iteration
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

$   write data to file
$   Note: we write at the first call of t, write the states for time t_old
    char fn_data[] = "circulation_model.csv";
    FILE *f_data;
    if (icall == 0){
        f_data=fopen(fn_data, "w");
        fprintf(f_data, "icall,time,pk,part,pven");
        fprintf(f_data, ",qart,qven,qk\n");
        fclose(f_data);
     }
    else if ( is_new_dt ) {
        f_data=fopen(fn_data, "a");
        fprintf(f_data, "%d,%.6e,",icall,t_old);
        fprintf(f_data, "%.6e,%.6e,%.6e,",pk2,part,pven);
        fprintf(f_data, "%.6e,%.6e,%.6e\n",qart2,qven2,qk2);
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