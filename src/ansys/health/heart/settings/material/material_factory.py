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
"""Factory for creating material models."""

from ansys.health.heart.settings.material.curve import constant_ca2
from ansys.health.heart.settings.material.material import (
    ACTIVE,
    ANISO,
    ISO,
    ActiveCurve,
    ActiveModel1,
    ActiveModel3,
    HGOFiber,
    Mat295,
)


def _default_myocardium_material() -> Mat295:
    """Get default Mat295 myocardium material."""
    from ansys.health.heart.settings.defaults.mechanics import material

    return _get_myocardium_material(material["myocardium"])


def _default_passive_material() -> Mat295:
    """Get default Mat295 passive material."""
    from ansys.health.heart.settings.defaults.mechanics import material

    return _get_passive_material(material["passive"])


def _get_myocardium_material(settings: dict, ep_coupled: bool = False) -> Mat295:
    """Get Mat295 myocardium material from settings."""
    if "isotropic" not in settings or "anisotropic" not in settings or "active" not in settings:
        raise ValueError("Incomplete myocardium material settings.")

    rho = settings["isotropic"]["rho"].m

    iso = ISO(
        kappa=settings["isotropic"]["kappa"].m,
        k1=settings["isotropic"]["k1"].m,
        k2=settings["isotropic"]["k2"].m,
        beta=2,
    )

    fibers = [HGOFiber(k1=settings["anisotropic"]["k1f"].m, k2=settings["anisotropic"]["k2f"].m)]

    if "k1s" in settings["anisotropic"]:
        sheet = HGOFiber(k1=settings["anisotropic"]["k1s"].m, k2=settings["anisotropic"]["k2s"].m)
        fibers.append(sheet)

    if "k1fs" in settings["anisotropic"]:
        k1fs, k2fs = settings["anisotropic"]["k1fs"].m, settings["anisotropic"]["k2fs"].m
    else:
        k1fs, k2fs = None, None
    aniso = ANISO(fibers=fibers, k1fs=k1fs, k2fs=k2fs)

    max = settings["active"]["taumax"].m
    bt = settings["active"]["beat_time"].m
    ss = settings["active"]["ss"]
    sn = settings["active"]["sn"]

    if not ep_coupled:
        ac_mdoel = ActiveModel1(taumax=max)  # use default field in Model1 except taumax
        curve = ActiveCurve(func=constant_ca2(tb=bt), threshold=0.1, type="ca2")
        active = ACTIVE(
            ss=ss,
            sn=sn,
            model=ac_mdoel,
            ca2_curve=curve,
        )
    else:
        ac_mdoel = ActiveModel3(
            ca2ion50=0.001,
            n=2,
            f=0.0,
            l=1.9,
            eta=1.45,
            sigmax=max,  # MPa
        )

        active = ACTIVE(
            ss=ss,
            sn=sn,
            acthr=0.0002,
            model=ac_mdoel,
            ca2_curve=None,
        )

    return Mat295(rho=rho, iso=iso, aniso=aniso, active=active)


def _get_passive_material(passive_settings: dict) -> Mat295:
    """Get passive Mat295 from settings."""
    passive = Mat295(
        rho=passive_settings["rho"].m,
        iso=ISO(
            itype=passive_settings["itype"],
            beta=2,
            kappa=passive_settings["kappa"].m,
            mu1=passive_settings["mu1"].m,
            alpha1=passive_settings["alpha1"],
        ),
    )
    return passive
