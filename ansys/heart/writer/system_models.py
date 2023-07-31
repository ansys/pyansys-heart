"""Module contains system model class.

Note
----
Valid models include:
1. Open loop constant pre-after load
2. Open-loop 3-element Windkessel
3. Closed loop

"""

from ansys.dyna.keywords import keywords


def ed_load_template():
    """
    Define function template to apply ED pressure.

    arg0: define function ID
    arg1: define function name
    arg2: target pressure
    arg3: flow if pressure is not reached
    """
    template = (
        "*DEFINE_FUNCTION\n"
        "{0}\n"
        "float {1}(float t, float dp, float area) \n"
        "{{\n"
        "float flow1;\n"
        "if (dp <= {2})\n"
        "{{\n"
        "flow1= {3};\n"
        "}} else\n"
        "{{\n"
        "flow1= 0.0;\n"
        "}}\n"
        "return flow1;\n"
        "}}"
    )
    return template


def windkessel_template():
    """Windkessel template.

    Note
    ----
    p_pre          p_fem                 p_art
    o--|Rv|--->|---|FEM|-->|-------|Ra| --- o -----+
                                            |      |
                                            |      |
                                          -----    |
                                           Ca     |Rp|
                                          -----    |
                                            |      |
                                p_ven       |      |
                                  o--------- o -----+

    """
    template = (
        "float {0}(float t, float dp)\n"
        "{{\n"
        "$   numerical constant\n"
        "    int Implicit={1};\n"
        "$   Only used for Euler Implicit\n"
        "    float gamma = 0.6;\n"
        "\n"
        "$   physical constants\n"
        "    float Rp, Ca;\n"
        "    float Ra, Rv;\n"
        "    Rp = {2};\n"  # "    Rp = 1.2e-4;\n"
        "    Ca = {3};\n"  # "    Ca = 2.5e4;\n"
        "    Ra = {4};\n"  # "    Ra = 1.0e-5;\n"
        "    Rv = {5};\n"  # "    Rv = 5.0e-6;\n"
        "\n"
        "$   physical variables\n"
        "    float chi_av, chi_mv;\n"
        "    float pk, part, pven;\n"
        "    float qk, qven, qart, qp;\n"
        "    float vart;\n"
        "\n"
        "$   only for save data purpose\n"
        "    float pk2, part2, pven2;\n"
        "    float qk2, qven2, qart2, qp2;\n"
        "    float vart2;\n"
        "\n"
        "$   constant pre load:\n"
        "    pven = {6};\n"
        "\n"
        "$   time related variables\n"
        "    int icall=0, is_new_dt=0;\n"
        "$   t: current time (input)\n"
        "$   t_last: time of the last call\n"
        "$   t_old: time of the last step\n"
        "    float t_last=0.0, t_old=0.0, dt;\n"
        "\n"
        "$   initialisation at t=0\n"
        "    if (icall == 0) {{\n"
        "          part = {7};\n"
        "$   initial arterial volume\n"
        "          vart = Ca * part;\n"
        "          qp = part / Rp;\n"
        "    }}\n"
        "\n"
        "$   determine if function call is at start of timestep:\n"
        "    if ( t-t_last > 1e-9 ) {{\n"
        "        is_new_dt = 1;\n"
        "        dt = t - t_old;\n"
        "    }}\n"
        "    else if ( t-t_last < 0. ) {{\n"
        '        printf("  ## Warning bisection may not be properly handled ##");\n'
        '        printf("  ## Warning: dt_old: %f", dt );\n'
        "        is_new_dt = 0;\n"
        "        dt = dt - (t_last-t);\n"
        '        printf("## Warning: dt_new: %f", dt );\n'
        "$       abort(0);\n"
        "    }} else\n"
        "    {{\n"
        "        is_new_dt = 0;\n"
        "    }}\n"
        "\n"
        "    if ( is_new_dt ) {{\n"
        "$   Save system states of last time step (at t_old)\n"
        "$   The converged pressure value of last time step (at t_old)\n"
        "        pk2 = dp;\n"
        "        part2 = part;\n"
        "        pven2 = pven;\n"
        "\n"
        "        vart2 = vart;\n"
        "\n"
        "        qp2 = qp;\n"
        "        qart2 = qart;\n"
        "        qven2 = qven;\n"
        "        qk2 = qk;\n"
        "\n"
        "$   Update system states for new time step (at t)\n"
        "        vart = vart + dt * (qart-qp);\n"
        "        part = vart / Ca;\n"
        "        qp = part / Rp;\n"
        "    }}\n"
        "    if (Implicit){{\n"
        "$   LSDYNA will integrate cavity volume implicitly: V^t = V^t_old+dt*Q^t\n"
        "$   LSDYNAs input dp is interpolated by dp=(1-r)*p^t_old+r*p^t+1_i\n"
        "$   This is not suitable to check the valve opening (to compute Q at t)\n"
        "$   We retrieve firstly p^t at this iteration\n"
        "        pk = (dp -(1-gamma)*pk2)/gamma;\n"
        "    }} else\n"
        "    {{\n"
        "$   LSDYNA will integrate cavity volume explicitly: V^t = V^t_old+dt*Q^t_old\n"
        "        pk = pk2;\n"
        "    }}\n"
        "    \n"
        "$   Update valve indicator functions\n"
        "    if (pven >= pk )\n"
        "    {{\n"
        "        chi_mv = {9};\n"
        "    }} else\n"
        "    {{\n"
        "        chi_mv = 1.e-16;\n"
        "    }}\n"
        "    if ( pk >= part )\n"
        "    {{\n"
        "        chi_av = {9};\n"
        "    }} else {{\n"
        "        chi_av = 1.e-16;\n"
        "    }}\n"
        "\n"
        "$   compute flow: In - Out\n"
        "    qven = chi_mv * ( ( pven - pk ) / Rv );\n"
        "    qart = chi_av * ( ( pk  - part) / Ra );\n"
        "    qk  = qven - qart;\n"
        "\n"
        "$   used to debug\n"
        "$   Note: we write at every call of t, write states for time t\n"
        '    char fn_bug[] = "circulation_model_debug.csv";\n'
        "    FILE *f_bug;\n"
        "    if (icall == 0){{\n"
        '        f_bug=fopen(fn_bug, "w");\n'
        '        fprintf(f_bug, "icall,is_new_dt,time,dp,pk,part,pven");\n'
        '        fprintf(f_bug, ",qart,qp,qven,qk,vart\\n");      \n'
        '        fprintf(f_bug, "%d,%d,%.6e,%.6e,",icall,is_new_dt,t,dp);\n'
        '        fprintf(f_bug, "%.6e,%.6e,%.6e,",pk,part,pven);\n'
        '        fprintf(f_bug, "%.6e,%.6e,%.6e,%.6e,",qart,qp,qven,qk);\n'
        '        fprintf(f_bug, "%.6e\\n",vart);\n'
        "      }}\n"
        "    else {{\n"
        '        f_bug=fopen(fn_bug, "a");\n'
        '        fprintf(f_bug, "%d,%d,%.6e,%.6e,",icall,is_new_dt,t,dp);\n'
        '        fprintf(f_bug, "%.6e,%.6e,%.6e,",pk,part,pven);\n'
        '        fprintf(f_bug, "%.6e,%.6e,%.6e,%.6e,",qart,qp,qven,qk);\n'
        '        fprintf(f_bug, "%.6e\\n",vart);\n'
        "    }}\n"
        "    fclose(f_bug);\n"
        "\n"
        "$   write data to file\n"
        "$   Note: we write at the first call of t, write the states for time t_old\n"
        '    char fn_data[] = "{8}";\n'
        "    FILE *f_data;\n"
        "    if (icall == 0){{\n"
        '        f_data=fopen(fn_data, "w");\n'
        '        fprintf(f_data, "icall,time,pk,part,pven");\n'
        '        fprintf(f_data, ",qart,qp,qven,qk,vart\\n"); \n'
        "        fclose(f_data);\n"
        "     }}\n"
        "    else if ( is_new_dt ) {{\n"
        '        f_data=fopen(fn_data, "a");\n'
        '        fprintf(f_data, "%d,%.6e,",icall,t_old);\n'
        '        fprintf(f_data, "%.6e,%.6e,%.6e,",pk2,part2,pven2);\n'
        '        fprintf(f_data, "%.6e,%.6e,%.6e,%.6e,",qart2,qp2,qven2,qk2);\n'
        '        fprintf(f_data, "%.6e\\n",vart2);\n'
        "        fclose(f_data);\n"
        "$       \n"
        "        t_old = t;\n"
        "    }}\n"
        "    \n"
        "$   Update counters\n"
        "    t_last = t;\n"
        "    icall = icall + 1;\n"
        "    \n"
        "$   LSDYNA defines outflow as positive\n"
        "    return -qk;\n"
        "}}\n"
    )
    return template


def define_function_windkessel(
    function_id: int = 10,
    function_name: str = "constant_preload_windkessel_after_load",
    implicit: bool = True,
    constants: dict = {"Rv": 5.0e-6, "Ra": 1.0e-5, "Rp": 1.2e-4, "Ca": 2.5e4, "Pven": 2},
    initialvalues: dict = {"part_init": 8},
    ivc=False,
):
    """Generate a Windkessel define function.

    Note
    ----
    Yields a formatted define function with a constant pre-load and WK afterload.

    """
    if implicit:
        implicit_flag = 1
    else:
        implicit_flag = 0

    if ivc is False:
        cf = 1
    else:
        # IVC: flow can not pass anyway
        cf = 1e-16

    filename_output = function_name + ".csv"

    wk_template = windkessel_template()

    define_function_str = wk_template.format(
        function_name,
        implicit_flag,
        constants["Rp"],
        constants["Ca"],
        constants["Ra"],
        constants["Rv"],
        constants["Pven"],
        initialvalues,
        filename_output,
        cf,
    )
    # format keyword:
    define_function_kw = keywords.DefineFunction(fid=function_id, function=define_function_str)
    define_function_kw.heading = function_name

    return define_function_kw


if __name__ == "__main__":
    print()
    define_function_windkessel()
