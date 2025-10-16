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

from ansys.health.heart import LOG as LOGGER
import ansys.health.heart.models as models
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
    MechanicalMaterialModel,
)


def _default_myocardium_material(ep_coupled: bool = False) -> Mat295:
    """Get default Mat295 myocardium material.

    Parameters
    ----------
    ep_coupled : bool, default: False
        Whether to use electro-mechanical coupling in material configuration.

    Returns
    -------
    Mat295
        The default myocardium mechanical material model.

    Raises
    ------
    ValueError
        If material settings are incomplete or invalid.
    RuntimeError
        If material creation fails.

    Examples
    --------
    >>> material = _default_myocardium_material(ep_coupled=True)
    >>> print(type(material))
    <class 'ansys.health.heart.settings.material.material.Mat295'>
    """
    try:
        from ansys.health.heart.settings.defaults.mechanics import material

        return _get_myocardium_material(material["myocardium"], ep_coupled=ep_coupled)

    except Exception as e:
        error_msg = f"Failed to create default myocardium material (EP coupled: {ep_coupled}): {e}"
        LOGGER.error(error_msg)
        raise RuntimeError(error_msg) from e


def _default_passive_material() -> Mat295:
    """Get default Mat295 passive material.

    Returns
    -------
    Mat295
        The default passive mechanical material model.

    Raises
    ------
    ValueError
        If material settings are incomplete or invalid.
    RuntimeError
        If material creation fails.

    Examples
    --------
    >>> material = _default_passive_material()
    >>> print(type(material))
    <class 'ansys.health.heart.settings.material.material.Mat295'>
    """
    try:
        from ansys.health.heart.settings.defaults.mechanics import material

        return _get_passive_material(material["passive"])

    except Exception as e:
        error_msg = f"Failed to create default passive material: {e}"
        LOGGER.error(error_msg)
        raise RuntimeError(error_msg) from e


def _get_myocardium_material(settings: dict, ep_coupled: bool = False) -> Mat295:
    """Get Mat295 myocardium material from settings.

    Parameters
    ----------
    settings : dict
        Material settings dictionary containing isotropic, anisotropic, and active properties.
    ep_coupled : bool, default: False
        Whether to use electro-mechanical coupling in active model configuration.

    Returns
    -------
    Mat295
        Configured myocardium material with isotropic, anisotropic, and active components.

    Raises
    ------
    ValueError
        If required settings are missing or incomplete.

    Examples
    --------
    >>> settings = {"isotropic": {...}, "anisotropic": {...}, "active": {...}}
    >>> material = _get_myocardium_material(settings, ep_coupled=True)
    >>> print(material.active is not None)
    True
    """
    # Validate required settings
    required_sections = ["isotropic", "anisotropic", "active"]
    missing_sections = [section for section in required_sections if section not in settings]
    if missing_sections:
        raise ValueError(
            f"Incomplete myocardium material settings. "
            f"Missing sections: {', '.join(missing_sections)}"
        )

    try:
        rho = settings["isotropic"]["rho"].m

        iso = ISO(
            kappa=settings["isotropic"]["kappa"].m,
            k1=settings["isotropic"]["k1"].m,
            k2=settings["isotropic"]["k2"].m,
            beta=2,
        )

        fibers = [
            HGOFiber(k1=settings["anisotropic"]["k1f"].m, k2=settings["anisotropic"]["k2f"].m)
        ]

        if "k1s" in settings["anisotropic"]:
            sheet = HGOFiber(
                k1=settings["anisotropic"]["k1s"].m, k2=settings["anisotropic"]["k2s"].m
            )
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

    except KeyError as e:
        raise ValueError(f"Missing required material property: {e}") from e
    except Exception as e:
        raise ValueError(f"Failed to construct myocardium material: {e}") from e


def _get_passive_material(passive_settings: dict) -> Mat295:
    """Get passive Mat295 from settings.

    Parameters
    ----------
    passive_settings : dict
        Passive material settings dictionary containing density and isotropic properties.

    Returns
    -------
    Mat295
        Configured passive material with isotropic properties only.

    Raises
    ------
    ValueError
        If required settings are missing or construction fails.

    Examples
    --------
    >>> settings = {"rho": 1.0, "itype": 1, "kappa": 1000.0, "mu1": 10.0, "alpha1": 18.5}
    >>> material = _get_passive_material(settings)
    >>> print(material.active is None)
    True
    """
    try:
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

    except KeyError as e:
        raise KeyError(f"Missing required passive material property: {e}") from e
    except Exception as e:
        raise ValueError(f"Failed to construct passive material: {e}") from e


def assign_default_mechanics_materials(
    model: models.HeartModel
    | models.FullHeart
    | models.BiVentricle
    | models.LeftVentricle
    | models.FourChamber,
    ep_coupled: bool = False,
) -> None:
    """Assign default mechanics materials to heart model parts in-place.

    This function modifies the input model by assigning appropriate mechanical
    materials to parts that do not have materials assigned. Parts with fiber
    orientation receive myocardium materials, while parts without fibers receive
    passive materials.

    Parameters
    ----------
    model : HeartModel or FullHeart or BiVentricle or LeftVentricle or FourChamber
        The heart model to assign materials to. Modified in-place.
    ep_coupled : bool, default: False
        Whether to use electro-mechanical coupling in material assignment.

    Raises
    ------
    ValueError
        If material creation fails for any part.
    RuntimeError
        If material assignment fails due to part validation errors.

    Examples
    --------
    >>> import ansys.health.heart.models as models
    >>> from ansys.health.heart.settings.material.material_factory import (
    ...     assign_default_mechanics_materials,
    ... )
    >>> model = models.BiVentricle()
    >>> assign_default_mechanics_materials(model, ep_coupled=True)
    >>> all_parts_have_materials = all(part.meca_material is not None for part in model.parts)
    >>> print(all_parts_have_materials)
    True

    Notes
    -----
    This function follows the PyAnsys Heart material assignment hierarchy:

    - Parts with fiber orientation receive myocardium materials with active components
    - Parts without fiber orientation receive passive materials only
    - Active parts retain active material properties, passive parts have active components disabled
    - The function will skip parts that already have valid mechanical materials assigned
    """
    assignments = 0

    for part in model.parts:
        if (
            not isinstance(part.meca_material, MechanicalMaterialModel)
            or part.meca_material is None
        ):
            try:
                if part.fiber:
                    part.meca_material = _default_myocardium_material(ep_coupled=ep_coupled)
                    LOGGER.info(
                        f"Assigned myocardium mechanical material to part '{part.name}' "
                        f"(EP coupled: {ep_coupled})."
                    )

                    # Disable the active module for passive parts
                    if not part.active:
                        part.meca_material.active = None
                        LOGGER.debug(f"Disabled active component for passive part '{part.name}'.")
                else:
                    part.meca_material = _default_passive_material()
                    LOGGER.info(f"Assigned passive mechanical material to part '{part.name}'.")

                assignments += 1

            except Exception as e:
                error_msg = f"Failed to assign mechanical material to part '{part.name}': {e}"
                LOGGER.error(error_msg)
                raise RuntimeError(error_msg) from e

    LOGGER.info(
        f"Successfully assigned default mechanical materials to "
        f"{assignments}/{len(model.parts)} parts "
    )

    return
