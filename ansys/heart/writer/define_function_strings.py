# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
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
