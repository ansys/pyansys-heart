"""Collection of define function strings."""


def function_alpha(alpha_endo: float = -60, alpha_epi: float = 60):
    """Define the alpha angle for fiber definition."""
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


def function_beta(beta_endo: float = 25, beta_epi: float = -65):
    """Define the beta angle for fiber definition in ventricles."""
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


def function_beta_septum(beta_endo: float = -65, beta_epi: float = 25):
    """Define the beta angle for fiber definition in the septum."""
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
