*DEFINE_FUNCTION
17
float bv_closed_rv(float t, float dp)
{
$   numerical constant
    int Implicit=1;
$   Only used for Euler Implicit
    float gamma = 0.6;
$   Pre filling duration
    t_pre_fill = 1.0;

$   physical constants
    float Ra, Rv;
    float Rp, Ca, Cv;
	float Va0, Vv0;
    Ra = 4.46e-6;
    Rv = 1.10e-6;
    Rp = 1.49e-4;
    Ca = 1.53e4;
	Cv = 4.59e4;
	Va0 = 7.04e5;
	Vv0 = 3.16e6;
$   pulmonary circulation
	float R_pa, R_pv;
    float R_pp, C_pa, C_pv;
	float V_pa0, V_pv0;
    R_pa = 2.48e-6;
    R_pv = 2.18e-6;
    R_pp = 1.78e-5;
    C_pa = 4.59e4;
	C_pv = 1.53e4;
	V_pa0 = 7.83e4;
	V_pv0 = 5.13e5;
$   Only used during the initialisation 
$   (p_art and p_part are not important since we impose constant filling)
$   But they determine how much blood in the Closed loop
$   Vtot = Vlv0+Va0+Ca*part+Vv0+Cv*pven;  ~ 5000 mL
    float part0, pven0, p_part0, p_pven0;
    part0 = 8;
	pven0 = 0.533;
	p_part0 = 0;
	p_pven0 = 2;
	
$   physical variables
$   means different things for lv and rv call
    float chi_out, chi_in;
	float pvv, qvv;

    float part, pven;
    float qven, qart, qp;
    float vart, vven;
	
$   pulmonary circulation
    float p_part, p_pven;
    float q_pven, q_part, q_pp;
    float v_part, v_pven;

$   only for save data purpose
	float pvv2, qvv2;

    float part2, pven2;
    float qven2, qart2, qp2;
    float vart2, vven2;
	float p_part2, p_pven2;
    float q_pven2, q_part2, q_pp2;
    float v_part2, v_pven2;

$   time related variables
    int icall=0, is_new_dt=0;
$   t: current time (input)
$   t_last: time of the last call
$   t_old: time of the last step
    float t_last=0.0, t_old=0.0, dt;

$   temporal variable used to exchange info for two cavities
	double arr_r[2];
	
$   initialisation at t=0
    if (icall == 0) {												   
          part = part0;
		  pven = pven0;
          p_part = p_part0;
		  p_pven = p_pven0;
		  
$   initial volume
          vart = Va0 + Ca * part;
          vven = Vv0 + Cv * pven;
          v_part = V_pa0 + C_pa * p_part;
          v_pven = V_pv0 + C_pv * p_pven;

$   initial flow		  
          qp = (part-pven) / Rp;
          q_pp = (p_part-p_pven) / R_pp;
		  qart = 0;
		  q_pven = p_pven / R_pv;
		  q_part = 0;
		  qven = pven / Rv;
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
$   Read flow from LV cavity (converged value at t_old)	
		FILE *fid2;
		fid2 = fopen("LV_flow.bin", "rb");
		fread( arr_r, sizeof(arr_r), 2, fid2 );
		fclose(fid2);	
		qart = arr_r[0];
		q_pven = arr_r[1];
		
$   The converged pressure value of last time step (at t_old)
        pvv2 = dp;
        part2 = part;
        pven2 = pven;
        p_part2 = p_part;
        p_pven2 = p_pven;

        vart2 = vart;
        vven2 = vven;
        v_part2 = v_part;
        v_pven2 = v_pven;

        qp2 = qp;
        qart2 = qart;
        qven2 = qven;
        q_pp2 = q_pp;
        q_part2 = q_part;
        q_pven2 = q_pven;
        qvv2 = qvv;

$   Update system states for new time step (at t)
        vart = vart + dt * (qart-qp);
        vven = vven + dt * (qp-qven);
        part = (vart-Va0) / Ca;
		if (t>=t_pre_fill){
			pven = (vven-Vv0) / Cv;
		} else{
			pven = pven0;
		}
        qp = (part-pven) / Rp;
$	rv		
		v_part = v_part + dt * (q_part-q_pp);
        v_pven = v_pven + dt * (q_pp-q_pven);
        p_part = (v_part-V_pa0) / C_pa;
		if (t>=t_pre_fill){
			p_pven = (v_pven-V_pv0) / C_pv;
		} else{
			p_pven = p_pven0;
		}
        q_pp = (p_part-p_pven) / R_pp;
    }
    if (Implicit){
$   LSDYNA will integrate cavity volume implicitly: V^t = V^t_old+dt*Q^t
$   LSDYNA's input dp is interpolated by dp=(1-r)*p^t_old+r*p^t+1_i
$   This is not suitable to check the valve's opening (to compute Q at t)
$   We retreive firstly p^t at this iteration
        pvv = (dp -(1-gamma)*pvv2)/gamma;
    } else
    {
$   LSDYNA will integrate cavity volume explicitly: V^t = V^t_old+dt*Q^t_old
        pvv = pvv2;
    }


$   Update valve indicator functions
    if (pven >= pvv )
    {
        chi_in = 1;
    } else
    {
        chi_in = 1.e-16;
    }
    if ( pvv >= p_part )
    {
        chi_out = 1;
    } else {
        chi_out = 1.e-16;
    }

$   compute flow: In - Out
    qven = chi_in * ( ( pven - pvv ) / Rv );
    q_part = chi_out * ( ( pvv  - p_part) / R_pa );
    qvv  = qven - q_part;

$   Write flow for RV cavity (Only the converged value will be read)	
	arr_r[0] = q_part;
	arr_r[1] = qven;
	FILE *fid1;
	fid1 = fopen("RV_flow.bin", "wb");
	fwrite( arr_r, sizeof(arr_r), 2, fid1 );
	fclose(fid1);
		
$   write data to file
$   Note: we write at the first call of t, write the states for time t_old
    char fn_data[] = "circulation_model_r.csv";
    FILE *f_data;
    if (icall == 0){
        f_data=fopen(fn_data, "w");
        fprintf(f_data, "icall,time,pvv,part,pven,p_part,p_pven,");
        fprintf(f_data, "qvv,qart,qp,qven,q_part,q_pp,q_pven,"); 
        fprintf(f_data, "vart,vven,v_part,v_pven\n"); 
        fclose(f_data);
     }
    else if ( is_new_dt ) {
        f_data=fopen(fn_data, "a");
        fprintf(f_data, "%d,%.6e,",icall,t_old);
        fprintf(f_data, "%.6e,%.6e,%.6e,",pvv2,part2,pven2);
        fprintf(f_data, "%.6e,%.6e,",p_part2,p_pven2);
        fprintf(f_data, "%.6e,%.6e,%.6e,%.6e,",qvv2,qart2,qp2,qven2);
        fprintf(f_data, "%.6e,%.6e,%.6e,",q_part2,q_pp2,q_pven2);
        fprintf(f_data, "%.6e,%.6e,",vart2,vven2);
        fprintf(f_data, "%.6e,%.6e\n",v_part2,v_pven2);
        fclose(f_data);
$       
        t_old = t;
    }
    
$   Update counters
    t_last = t;
    icall = icall + 1;
    
$   LSDYNA defines outflow as positive
    return -qvv;
}
*END