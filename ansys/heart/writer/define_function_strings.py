"""Collection of define function strings"""


function1 = "\n".join( [
    "float alpha(",
    "            float x_ele, float y_ele, float z_ele,",
    "            float phi_len, float phi_thi)",
    "{ ",
    "  float alpha1;",
    "  float pi;",
    "  float alpha_endo;",
    "  float alpha_epi;",
    "  pi=3.14159265359;",
    "  alpha_endo=-60*pi/180;",
    "  alpha_epi=60*pi/180;",
    "  alpha1=alpha_endo*(1-phi_thi)+alpha_epi*phi_thi;",
    "  return alpha1;",
    "}",
] )

function2 = "\n".join( [
    "    float beta(",
    "            float x_ele, float y_ele, float z_ele,",
    "            float phi_len, float phi_thi)",
    "{  ",
    "  float beta1;",
    "  float pi;",
    "  float beta_endo;",
    "  float beta_epi;",
    "  pi=3.14159265359;",
    "  beta_epi=(25)*pi/180;",
    "  beta_endo=(-65)*pi/180;",
    "  beta1=beta_endo*(1-phi_thi)+beta_epi*phi_thi;",
    "  return beta1;",
    "}",
] )

function3 = "\n".join( [
    "    float betaW(",
    "            float x_ele, float y_ele, float z_ele,",
    "            float phi_len, float phi_thi)",
    "{  ",
    "  float beta1;",
    "  float pi;",
    "  float beta_endo;",
    "  float beta_epi;",
    "  pi=3.14159265359;",
    "  beta_epi=(25)*pi/180;",
    "  beta_endo=(-65)*pi/180;",
    "  beta1=beta_endo*(1-phi_thi)+beta_endo*phi_thi;",
    "  return beta1;",
    "}",
] )