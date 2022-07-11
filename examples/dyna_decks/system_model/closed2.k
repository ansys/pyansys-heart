$ alternative method to update circulation model
$ This method takes the initial volume of cavity(LV) as input
*DEFINE_FUNCTION
16
float lv_closed_load(float t, float dp)
{
$   numerical constant
    int Implicit=1;
$   Only used for Euler Implicit
    float gamma = 0.6;

$   physical constants
    float Ra, Rv;
    float Rp, Ca, Cv;
	float Va0, Vv0, Vlv0, Vtot;
    Ra = 1.0e-5;
    Rv = 5.0e-6;
    Rp = 1.2e-4;
    Ca = 2.5e4;
	Cv = 6.0e5;
	Va0 = 5.0e5;
	Vv0 = 3.0e6;
	Vlv0 = 79234;
	Vtot = Vlv0+Va0+Ca*8+Vv0+Cv*2;
	
$   physical variables
    float chi_av, chi_mv;
    float pk, part, pven;
    float qk, qven, qart, qp;
    float vart, vven, vlv;

$   only for save data purpose
    float pk2, part2, pven2;
    float qk2, qven2, qart2, qp2;
    float vart2, vven2,vlv2;

$   time related variables
    int icall=0, is_new_dt=0;
$   t: current time (input)
$   t_last: time of the last call
$   t_old: time of the last step
    float t_last=0.0, t_old=0.0, dt;

$   initialisation at t=0
    if (icall == 0) {
          part = 8;
$		  pven = 2;
$   initial volume
          vart = Va0 + Ca * part;
          vlv = Vlv0;
          vven = Vtot-vlv-vart;
          pven = (vven-Vv0) / Cv;		  
$   initial flow		  
          qp = part / Rp;
$   initial values for system update 		  
		  vart2 = vart;
		  vlv2 = vlv;
		  part2 = part;
		  pven2 = pven;
    }

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
        part2 = part;
        pven2 = pven;

        vart2 = vart;
        vven2 = vven;
        vlv2 = vlv;

        qp2 = qp;
        qart2 = qart;
        qven2 = qven;
        qk2 = qk;
    }
	
    if (Implicit){
$   Use pressure at t to update flow at t
$   LSDYNA's input dp is interpolated by dp=(1-r)*p^t_old+r*p^t+1_i
$   This is not suitable to check the valve's opening (to compute Q at t)
$   We retreive firstly p^t at this iteration
        pk = (dp -(1-gamma)*pk2)/gamma;
    } else
    {
$   Use pressure at t_old to update flow at t
$   This implies LSDYNA will integrate cavity volume explicitly:
$   V^t = V^t_old+dt*Q^t_old
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

$   Update system states for at t
$   will be save for the first call of next t
	float qf;
$   qf: the same method as LSDYNA to update Volume
	qf = (1-gamma)*(qven2-qart2)+gamma*(qven-qart);
	vlv = vlv2 + dt * qf;
$   Euler implicit to update other Volume
	vart = vart2 + dt * (qart-qp);
	vven = Vtot-vlv-vart;
	part = (vart-Va0) / Ca;
	pven = (vven-Vv0) / Cv;
	qp = (part-pven) / Rp;			 
	
$   used to debug
$   Note: we write at every call of t, write states for time t
    char fn_bug[] = "circulation_model_debug.csv";
    FILE *f_bug;
    if (icall == 0){
        f_bug=fopen(fn_bug, "w");
        fprintf(f_bug, "icall,is_new_dt,time,dp,pk,part,pven");
        fprintf(f_bug, ",qart,qp,qven,qk,vart,vven\n");      
        fprintf(f_bug, "%d,%d,%.6e,%.6e,",icall,is_new_dt,t,dp);
        fprintf(f_bug, "%.6e,%.6e,%.6e,",pk,part,pven);
        fprintf(f_bug, "%.6e,%.6e,%.6e,%.6e,",qart,qp,qven,qk);
        fprintf(f_bug, "%.6e,%.6e\n",vart,vven);
      }
    else {
        f_bug=fopen(fn_bug, "a");
        fprintf(f_bug, "%d,%d,%.6e,%.6e,",icall,is_new_dt,t,dp);
        fprintf(f_bug, "%.6e,%.6e,%.6e,",pk,part,pven);
        fprintf(f_bug, "%.6e,%.6e,%.6e,%.6e,",qart,qp,qven,qk);
        fprintf(f_bug, "%.6e,%.6e\n",vart,vven);
    }
    fclose(f_bug);

$   write data to file
$   Note: we write at the first call of t, write the states for time t_old
    char fn_data[] = "circulation_model.csv";
    FILE *f_data;
    if (icall == 0){
        f_data=fopen(fn_data, "w");
        fprintf(f_data, "icall,time,pk,part,pven");
        fprintf(f_data, ",qart,qp,qven,qk,vart,vven,vlv\n");
        fclose(f_data);
     }
    else if ( is_new_dt ) {
        f_data=fopen(fn_data, "a");
        fprintf(f_data, "%d,%.6e,",icall,t_old);
        fprintf(f_data, "%.6e,%.6e,%.6e,",pk2,part2,pven2);
        fprintf(f_data, "%.6e,%.6e,%.6e,%.6e,",qart2,qp2,qven2,qk2);
        fprintf(f_data, "%.6e,%.6e,%.6e\n",vart2,vven2,vlv2);
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