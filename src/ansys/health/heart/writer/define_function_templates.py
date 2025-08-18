# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Collection of define function templates."""

from ansys.dyna.core.keywords import keywords


def _function_alpha(alpha_endo: float = -60, alpha_epi: float = 60) -> str:
    """Define the alpha angle for the fiber definition."""
    return "\n".join(
        [
            "float alpha(",
            "            float x_ele, float y_ele, float z_ele,",
            "            float phi_len, float phi_thi)",
            "{ ",
            "  float alpha1;",
            "  float pi;",
            "  float alpha_endo;",
            "  float alpha_epi;",
            "  pi=3.14159265359;",
            "  alpha_endo={0:.2f}*pi/180;".format(alpha_endo),
            "  alpha_epi={0:.2f}*pi/180;".format(alpha_epi),
            "  alpha1=alpha_endo*(1-phi_thi)+alpha_epi*phi_thi;",
            "  return alpha1;",
            "}",
        ]
    )


def _function_beta(beta_endo: float = 25, beta_epi: float = -65) -> str:
    """Define the beta angle for the fiber definition in ventricles."""
    return "\n".join(
        [
            "    float beta(",
            "            float x_ele, float y_ele, float z_ele,",
            "            float phi_len, float phi_thi)",
            "{  ",
            "  float beta1;",
            "  float pi;",
            "  float beta_endo;",
            "  float beta_epi;",
            "  pi=3.14159265359;",
            "  beta_endo={0:.2f}*pi/180;".format(beta_endo),
            "  beta_epi={0:.2f}*pi/180;".format(beta_epi),
            "  beta1=beta_endo*(1-phi_thi)+beta_epi*phi_thi;",
            "  return beta1;",
            "}",
        ]
    )


def _function_beta_septum(beta_endo: float = -65, beta_epi: float = 25) -> str:
    """Define the beta angle for the fiber definition in the septum."""
    return "\n".join(
        [
            "    float betaW(",
            "            float x_ele, float y_ele, float z_ele,",
            "            float phi_len, float phi_thi)",
            "{  ",
            "  float beta1;",
            "  float pi;",
            "  float beta_endo;",
            "  float beta_epi;",
            "  pi=3.14159265359;",
            "  beta_endo={0:.2f}*pi/180;".format(beta_endo),
            "  beta_epi={0:.2f}*pi/180;".format(beta_epi),
            "  beta1=beta_endo*(1-phi_thi)+beta_endo*phi_thi;",
            "  return beta1;",
            "}",
        ]
    )


def _constant_flow_template() -> str:
    """Constant flow template."""
    template = "float {0}(float t, float dp, float area) \n{{\nreturn {1};\n}}"
    return template


def _valve_template() -> str:
    """Diode valve model template."""
    template = (
        "float {0}(float t, float dp)\n"
        "{{\n"
        "   float Rv, qven;\n"
        "   Rv = {1};\n"
        "   if (dp >= 0.0 )\n"
        "   {{\n"
        "       qven =  dp  / Rv ;\n"
        "   }}\n"
        "   else{{\n"
        "       qven=0.0;\n"
        "   }}\n"
        "   return qven;\n"
        "}}\n"
    )
    return template


def _afterload_windkessel_template() -> str:
    """Use the same model as in :func:`windkessel_template` but without preload."""
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
        "    Rp = {2:.12e};\n"  # "    Rp = 1.2e-4;\n"
        "    Ca = {3:.12e};\n"  # "    Ca = 2.5e4;\n"
        "    Ra = {4:.12e};\n"  # "    Ra = 1.0e-5;\n"
        "    Rv = {5:.12e};\n"  # "    not used !"
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
        "    pven = {6:.12e};\n"
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
        "          part = {7:.12e};\n"
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
        "        chi_mv = 1;\n"
        "    }} else\n"
        "    {{\n"
        "        chi_mv = 1.e-16;\n"
        "    }}\n"
        "    if ( pk >= part )\n"
        "    {{\n"
        "        chi_av = 1;\n"
        "    }} else {{\n"
        "        chi_av = 1.e-16;\n"
        "    }}\n"
        "\n"
        "$   compute flow: In - Out\n"
        "    qven = 0.0;\n"
        "    qart = chi_av * ( ( pk  - part) / Ra );\n"
        "    qk  = qven - qart;\n"
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


def _ed_load_template() -> str:
    """
    Define function template to apply ED pressure.

    Notes
    -----
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
        "if (dp <= {2:.12e})\n"
        "{{\n"
        "flow1= {3:.12e};\n"
        "}} else\n"
        "{{\n"
        "flow1= 0.0;\n"
        "}}\n"
        "return flow1;\n"
        "}}"
    )
    return template


def _windkessel_template() -> str:
    """Windkessel template.

    Notes
    -----
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
        "    int Implicit={1:d};\n"
        "$   Only used for Euler Implicit\n"
        "    float gamma = 0.6;\n"
        "\n"
        "$   physical constants\n"
        "    float Rp, Ca;\n"
        "    float Ra, Rv;\n"
        "    Rp = {2:.12e};\n"  # "    Rp = 1.2e-4;\n"
        "    Ca = {3:.12e};\n"  # "    Ca = 2.5e4;\n"
        "    Ra = {4:.12e};\n"  # "    Ra = 1.0e-5;\n"
        "    Rv = {5:.12e};\n"  # "    Rv = 5.0e-6;\n"
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
        "    pven = {6:.12e};\n"
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
        "          part = {7:.12e};\n"
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
        "        chi_mv = 1;\n"
        "    }} else\n"
        "    {{\n"
        "        chi_mv = 1.e-16;\n"
        "    }}\n"
        "    if ( pk >= part )\n"
        "    {{\n"
        "        chi_av = 1;\n"
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


def _define_function_0d_system(
    function_id: int,
    function_name: str,
    parameters: dict[str, dict[str, float]],
) -> keywords.DefineFunction:
    """Generate a define function.

    Parameters
    ----------
    function_id : int
        Function ID that defines the interaction between control volumes.
    function_name : str
        Function name.
    parameters : dict[str, dict[str, float]]
        Parameters of the system model.

    Returns
    -------
    keywords.DefineFunction
        Formatted keyword for the define function.
    """
    if "constant_preload_windkessel_afterload" in function_name:
        define_function_str = _windkessel_template().format(
            function_name,
            1,  # implicit flag
            parameters["constants"]["Rp"],
            parameters["constants"]["Ca"],
            parameters["constants"]["Ra"],
            parameters["constants"]["Rv"],
            parameters["constants"]["Pven"],
            parameters["initial_value"]["part"],
            function_name + ".csv",
        )
    elif "afterload_windkessel" in function_name:
        define_function_str = _afterload_windkessel_template().format(
            function_name,
            1,  # implicit flag
            parameters["constants"]["Rp"],
            parameters["constants"]["Ca"],
            parameters["constants"]["Ra"],
            parameters["constants"]["Rv"],  # not used
            parameters["constants"]["Pven"],  # not used
            parameters["initial_value"]["part"],
            function_name + ".csv",
        )
    elif "constant_flow" in function_name:
        define_function_str = _constant_flow_template().format(
            function_name,
            parameters["flow"],
        )
    elif "valve" in function_name:
        define_function_str = _valve_template().format(
            function_name,
            parameters["Rv"],
        )

    # format keyword:
    define_function_kw = keywords.DefineFunction(fid=function_id, function=define_function_str)
    define_function_kw.heading = function_name

    return define_function_kw
